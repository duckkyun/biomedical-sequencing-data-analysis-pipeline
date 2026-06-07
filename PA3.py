#!/usr/bin/env python3

import sys


def get_contig_order(vcf_file):
	"""VCF 헤더의 ##contig 라인에서 chromosome 이름들을 등장 순서대로 추출.
		list 형태로 저장하여 이후  sorted()로 처리하기 위한 준비."""
	order = []
	for line in open(vcf_file):
		if not line.startswith('#'):
			break
		if line.startswith('##contig'):
			chrom = line.split('ID=')[1].split(',')[0]
			order.append(chrom)
	return order


def chrom_key(c):
	"""숫자(1-22) → X → Y → M/MT → 기타 순으로 정렬용 튜플 반환.
	   sorted(..., key=chrom_key) 로 호출하면 튜플의 첫 요소로 그룹 분리,
	   두 번째 요소로 그룹 내 정렬이 이루어진다.
	   chr로 표기되는 format과 그렇지 않은 format이 있다고 알려져 있다.
	   따라서 두 경우 모두 cover하도록 했음."""
	
	s = c[3:] if c.startswith('chr') else c

	if s.isdigit():
		return (0, int(s))

	special = {'X': 1, 'Y': 2, 'M': 3, 'MT': 3}
	if s in special:
		return (1, special[s])

	return (2, s)


if len(sys.argv) != 2:
	print('Usage: ./PA3.py <vcf_file>')

else:
	vcf_file = sys.argv[1]

	myfile = open(vcf_file)
	counts = {}
	lengths = {}


# vcf file의 각 line 읽기
	for line in myfile:
		line = line.rstrip()

		if not line or line.startswith('#'):
			continue
		else:
			column = line.split('\t')
			chrom = column[0]
			pos = column[1]
			id_ = column[2]
			ref = column[3]
			alt = column[4]
			qual = float(column[5]) if column[5] != '.' else 0.0
			filter_ = column[6]
			info = column[7]

# INFO에는 여러 정보가 ;를 기준으로 separated되어 있기 때문에 split method를 활용해 list로 변환 후, list를 index별로 읽어 dictionary에 저장. DP는 항상  단일값만 가지기 때문에 이 과제에서 dictionary 내 value가 여러 값인 경우는 고려하지 않고 하나의 string으로 저장되게끔 하였다.
			info_dict = {}
			for item in info.split(';'):
				if '=' in item:
					key, value = item.split('=')
					info_dict[key] = value
# '='가 없는 정보는 다루지 않기 때문에 그대로 저장만 해둔다.
				else:
					info_dict[item] = item

# 처음 보는 chromosome 이면 counts와 lengths 딕셔너리에 추가.
			if chrom not in counts:
				counts[chrom]  = {'substitutions': 0, 'insertions': 0, 'deletions': 0}
				lengths[chrom] = {'substitutions': 0, 'insertions': 0, 'deletions': 0}

# 필터 조건: QUAL >= 20, DP >= 20, ALT 단일값, AF == 1
# .get() 으로 안전하게 조회 (DP/AF 가 없는 라인이면 0 으로 처리되어 자연스럽게 제외됨)
			if qual < 20:
				continue
			if int(info_dict.get('DP', '0')) < 20:
				continue
			if ',' in alt:
				continue
			if float(info_dict.get('AF', '0')) != 1:
				continue

# 처음 보는 chromosome 이면 counts와 lengths 딕셔너리에 추가.
			"""if chrom not in counts:
				counts[chrom]  = {'substitutions': 0, 'insertions': 0, 'deletions': 0}
				lengths[chrom] = {'substitutions': 0, 'insertions': 0, 'deletions': 0}"""

# 변이 종류 분류 후 카운트 & 길이 누적
			if len(ref) == 1 and len(alt) == 1:
				counts[chrom]['substitutions']  += 1
				lengths[chrom]['substitutions'] += 1
			elif len(alt) > len(ref):
				counts[chrom]['insertions']  += 1
				lengths[chrom]['insertions'] += len(alt) - len(ref)
			elif len(ref) > len(alt):
				counts[chrom]['deletions']  += 1
				lengths[chrom]['deletions'] += len(ref) - len(alt)

	myfile.close()


# 출력 ----------------------------------------------------------------
# 헤더의 ##contig 순서를 chrom_key로 정렬
	contig_order = get_contig_order(vcf_file)
# contig 헤더에 없지만 data line에 등장한 chromosome도 포함시키기 위해 union.
	sorted_chroms = sorted(set(contig_order) | set(counts.keys()), key=chrom_key)

# 전체 합계 누적용 카운터
	total_counts  = {'substitutions': 0, 'insertions': 0, 'deletions': 0}
	total_lengths = {'substitutions': 0, 'insertions': 0, 'deletions': 0}

# 표 형식 출력
	fmt = '{:<12} {:>14} {:>11} {:>10} {:>11} {:>11} {:>11}'
	separator = '-' * 86

	print(fmt.format('Chromosome', 'Substitutions', 'Insertions', 'Deletions',
	                 'Sub_Length', 'Ins_Length', 'Del_Length'))
	print(separator)

# chrom_key 함수로 sorting한 순서대로 출력    
	for c in sorted_chroms:
# contig에는 나왔지만 Data lines에는 나오지 않는 Edge Case 고려.
		if c not in counts:
			continue
# 각 Chromosome에 대한 dictionary 값 변수 저장.
		cc = counts[c]
		ll = lengths[c]
		print(fmt.format(c,
			cc['substitutions'], cc['insertions'], cc['deletions'],
			ll['substitutions'], ll['insertions'], ll['deletions']))
# 각 chromosome별로 정보가 저장돼 있기 때문에, sorted_chroms에서 하나씩 불러 모든 total dictionary의 key에 넣어주기만 하면 된다.
		for k in total_counts:
			total_counts[k]  += cc[k]
			total_lengths[k] += ll[k]

	print(separator)
	print(fmt.format('TOTAL',
		total_counts['substitutions'],  total_counts['insertions'],  total_counts['deletions'],
		total_lengths['substitutions'], total_lengths['insertions'], total_lengths['deletions']))
