#!/usr/bin/env python3

import sys
import math

if len(sys.argv) != 4:
    sys.exit("Usage: {} <gtf> <patients_ReadsPerGene.out.tab> <normals_ReadsPerGene.out.tab>".format(sys.argv[0]))

gtf = sys.argv[1]
patients_file = sys.argv[2]
normals_file = sys.argv[3]

# GTF: gene_id를 key로, (gene_name, chrom, start, end, length_bp)을 value로 갖는 dictionary 저장
gene_info = {}
genefile = open(gtf)
for line in genefile:
    line = line.rstrip()
    if not line or line.startswith('#'):
        continue
    column = line.split('\t')
    if column[2] != 'gene':
        continue
    chrom = column[0]
    start_p = int(column[3])
    end_p = int(column[4])
    length_bp = end_p - start_p + 1

    gene_id = None
    gene_name = None
    for item in column[8].split(';'):
        item = item.strip()
        if item.startswith('gene_id'):
            gene_id = item.split()[1].strip('"')
        elif item.startswith('gene_name'):
            gene_name = item.split()[1].strip('"')

    if gene_id:
        gene_info[gene_id] = (gene_name or gene_id, chrom, start_p, end_p, length_bp)
genefile.close()

# STAR ReadsPerGene.out.tab (column 2 = unstranded counts)
def read_star_counts(path, col=1):
    counts = {}
    f = open(path)
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        c = line.split('\t')
# Header line 생략
        if c[0].startswith('N_'):
            continue
# gene_id를 key로, column 값을 value로 dictionary 생성
        counts[c[0]] = int(c[col])
    f.close()
    return counts

patient_counts = read_star_counts(patients_file)
normal_counts  = read_star_counts(normals_file)

# 출력·TPM 분모에 쓸 gene 집합: 한쪽 sample에서라도 발현(count>0)된 gene 모두 포함 
genes = [g for g in patient_counts
         if g in normal_counts and g in gene_info
         and (patient_counts[g] > 0 or normal_counts[g] > 0)]
genes.sort()

# log, geometric mean 계산 전용 집합: 두 sample 모두 count>0 (geometric mean이 정의되는 gene만)
common = [g for g in genes
          if patient_counts[g] > 0 and normal_counts[g] > 0]

# ---- TPM ----
def calc_tpm(counts):
    rpk = {}
    for g in genes:
# RPK(Reads Per Kilobase) : 각 gene의 length_bp를 kb로 변환 (count 0이면 RPK 0 -> 해당 sample TPM 0)
        length_kb = gene_info[g][4] / 1000.0
# 각 gene_id마다 read counts를 gene length로normalize
        rpk[g] = counts[g] / length_kb
# 모든 값의 합을 1 million으로 나눠 scaled sum을 구한다
    scale = sum(rpk.values()) / 1000000.0
    if scale == 0:
        return {g: 0.0 for g in genes}
# 각 gene_id를 sequencing depth까지 고려한 dictionary를 생성
    return {g: v / scale for g, v in rpk.items()}

tpm_patient = calc_tpm(patient_counts)
tpm_normal  = calc_tpm(normal_counts)

# ---- Median-of-ratios  ----
# Step 1: geometric mean 구해서 ratio 구할 준비
geo_mean = {}
for g in common:
    pc, nc = patient_counts[g], normal_counts[g]
    if pc > 0 and nc > 0:
        geo_mean[g] = math.sqrt(pc * nc)      

# median: 홀수 개면 중앙값, 짝수 개면 중앙에 가까운 두 값의 평균 
# Step 1-2: (count / geometric mean) 비율의 median 구하기 위한 준비
def median(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:        # 홀수 개 → 정가운데 값
        return s[mid]
    else:                 # 짝수 개 → 가운데 두 값의 평균
        return (s[mid - 1] + s[mid]) / 2.0

def size_factor(counts):
# Step 2-3: 각 sample에서 (count / geometric mean) 비율을 구하고, 비율들의 median = size factor
    ratios = [counts[g] / mean for g, mean in geo_mean.items()]
    return median(ratios) if ratios else 1.0

sf_patient = size_factor(patient_counts)
sf_normal  = size_factor(normal_counts)
# Step 4: size factor로 normalized counts 계산
mr_patient = {g: patient_counts[g] / sf_patient for g in genes}
mr_normal  = {g: normal_counts[g]  / sf_normal  for g in genes}

# ---- log2 fold change (patient / normal) ----
def log2fc(a, b):
    # 한쪽이라도 0이면 정의 불가 → NA 처리용 None 반환
    if a <= 0 or b <= 0:
        return None
    return math.log2(a / b)

def fmt_num(x):
    # log/geometric mean이 정의 불가한 값(None)은 NA, 나머지는 소수점 2자리
    return "NA" if x is None else "{:.2f}".format(x)

# ---- print tabular output ----
header = ['gene_name', 'chromosome', 'start_position', 'end_position',
          'TPM_patient', 'TPM_normal', 'log2FC_TPM',
          'Median_of_ratio_patient',  'Median_of_ratio_normal',  'log2FC_MR']

# 모든 행을 먼저 문자열로 만들어 모은다 (열 너비 계산을 위해)
rows = []
for g in genes:
    name, chrom, start_p, end_p, _ = gene_info[g]
    tp, tn = tpm_patient[g], tpm_normal[g]
    mp, mn = mr_patient[g],  mr_normal[g]
    rows.append([
        name, chrom, str(start_p), str(end_p),
        fmt_num(tp), fmt_num(tn), fmt_num(log2fc(tp, tn)),
        fmt_num(mp), fmt_num(mn), fmt_num(log2fc(mp, mn)),
    ])

# 열별 최대 너비를 width에 저장= 헤더와 모든 데이터 셀 길이 중 최댓값 (정렬 유지 목적)
widths = [len(h) for h in header]
for r in rows:
    for i, cell in enumerate(r):
        if len(cell) > widths[i]:
            widths[i] = len(cell)

# gene_name, chromosome(앞 2열)은 왼쪽 정렬, 나머지 숫자 열은 오른쪽 정렬
def fmt_row(cells):
    out = []
    for i, cell in enumerate(cells):
        pad = ' ' * (widths[i] - len(cell))   # 모자란 만큼 공백 생성
        # 앞 2열은 왼쪽 정렬(공백을 뒤에), 나머지는 오른쪽 정렬(공백을 앞에)
        out.append(cell + pad if i < 2 else pad + cell)
    return '  '.join(out)   # 열 사이 구분은 공백 2칸

header_line = fmt_row(header)
print(header_line)
print('-' * len(header_line))   # 구분선
for r in rows:
    print(fmt_row(r))
