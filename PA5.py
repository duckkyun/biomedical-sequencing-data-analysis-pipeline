#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
	print("Usage: ./PA5.py <fasta_file>")
	sys.exit(1)


fa_file = sys.argv[1]

read_list = []
myfile = open(fa_file)

# multi-line FASTA 대응: 다음 '>' 헤더가 오기 전까지 sequence 길이를 누적했다가 append
current_length = 0
seen_header = False

for line in myfile:
	line = line.rstrip()

	if line.startswith('>'):
		if seen_header:                 # 직전 sequence를 마무리해 저장
			read_list.append(current_length)
		current_length = 0              # 새 sequence 시작 → 누적 초기화
		seen_header = True

	else:
		current_length += len(line)     # 같은 sequence의 여러 줄을 합산

if seen_header:                         # 파일 끝의 마지막 sequence 저장
	read_list.append(current_length)

myfile.close()

if not read_list:                       # FASTA에 sequence가 하나도 없으면 종료
	print("No sequences found in %s" % fa_file)
	sys.exit(1)

counts = len(read_list)         # 1. The number of total sequences (= 헤더 개수)
total_length = sum(read_list) # 2. Total length of sequences
max_length = max(read_list) # 3. The length of the longest sequence
min_length = min(read_list) # 4. The length of the shortest sequence


# 5. N50을 찾자
read_list.sort(reverse=True) # 내림차순으로 배열
pseudo_total = 0
for i in read_list:
	pseudo_total = pseudo_total + i
	if pseudo_total > (total_length) * (0.5):
		N50 = i
		break


# 6. sequence distribution 표현 고민하기
length_dict = {}

for i in read_list:
	b = len(str(i)) # 길이의 자릿수
	length_dict[b] = length_dict.get(b,0) + 1 # 해당 자릿수 value 1 추가
"""
for i in read_list:
	if i in length_dict: # length_dict에 이미 있는 read_length라면, 저장돼 있는 key에 count 하나 추가.
		length_dict[i] = length_dict[i] + 1
	else: # length_dict에 없는 read_length라면, key 하나를 만들고 value도 1 부여.
		length_dict[i] = 1
"""


print("Total length of sequences is %s" % (total_length))
print("The number of total sequences is %s" %  (counts))
print("The length of the longest sequences is %s" %  (max_length))
print("The length of the shortest sequences is %s" %  (min_length))
print("N50 is %s" % (N50))
print('%12s %14s %13s' % ('<range>', '<counts>', '<percentage>'))
for length in sorted(length_dict):
	low = 10 ** (length-1)
	high = (10 ** length) - 1
	cnt = length_dict[length]
	pct = cnt / counts * 100
	print('%6d to %-7d %8d %11.2f%%' % (low, high, cnt, pct))
