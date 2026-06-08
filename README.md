# Biomedical sequencing data analysis pipeline

건국대학교(Konkuk University) **Biomedical Data Processing (Spring 2026)** 강좌의
학부 수업 과제용(Programming Assignment) 프로그램이다. NGS sequencing 데이터
분석 파이프라인의 각 단계(QC → mapping → variant calling / DEG / assembly 평가)를
외부 패키지 없이 Python 표준 라이브러리만으로 직접 구현했다.

> 학부 수업 과제용 프로그램 (AI 기반 검색을 통해 직접 작성하였음)

## Environment

- **Python 3.5.2** — Konkuk `bps.konkuk.ac.kr` server (`/usr/bin/python3`)에서 작성·실행
- 표준 라이브러리만 사용 (`sys`, `math`) — 외부 의존성 없음 (`requirements.txt` 불필요)
- f-string 미사용 → **Python 3.5 이상**이면 그대로 실행 가능

## Scripts

| Script | Input | Description |
|---|---|---|
| `PA1.py` | FASTQ | state machine으로 FASTQ를 파싱해 read 수 · 총 길이 · read length 분포(텍스트 + 막대그래프) 출력 |
| `PA2.py` | SAM | FLAG bit로 properly-aligned read pair를 PE / MP로 분류, mean TLEN · orientation fraction · insert type 추정 |
| `PA3.py` | VCF | QUAL ≥ 20 · DP ≥ 20 · AF = 1 필터 후 chromosome별 substitution / insertion / deletion 개수·길이 집계 (염색체 자연 정렬). 추가로 filtering funnel(필터 단계별 변이 잔존), Ti/Tv ratio, indel frameshift 비율 등 callset QC 지표 산출 |
| `PA4.py` | GTF + STAR `ReadsPerGene.out.tab` ×2 | patients vs normals에서 TPM · median-of-ratios normalization · log2 fold change 계산 (DEG) |
| `PA5.py` | FASTA | assembly 통계 — total length · sequence 수 · longest / shortest · N50 · length 분포 |

## Usage

```bash
./PA1.py <fastq_file>
./PA2.py <sam_file>
./PA3.py <vcf_file>
./PA4.py <gtf> <patients_ReadsPerGene.out.tab> <normals_ReadsPerGene.out.tab>
./PA5.py <fasta_file>
```

각 스크립트는 인자 없이 실행하면 usage 메시지를 출력한다.

## Pipeline context

| Track | 단계 | 해당 스크립트 |
|---|---|---|
| ① DNA-seq (WES) variant calling | QC → mapping → variant calling | `PA1` (QC) · `PA2` (mapping) · `PA3` (variant) |
| ② RNA-seq (DEG) | gene count → 정규화 → fold change | `PA4` |
| ③ Genome assembly | de novo assembly 평가 | `PA5` |
