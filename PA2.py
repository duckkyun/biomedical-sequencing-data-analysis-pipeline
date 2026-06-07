#!/usr/bin/env python3

import sys

# 실행 파일 + sam file을 각 argument로 받아야만 실행
if len(sys.argv) != 2:
    print('Usage: ./PA2.py <sam_file>')

else:
	sam_file = sys.argv[1]

	myfile = open(sam_file)
	count_all_read = 0
	header_line = 0
	PE = 0
	MP = 0
	RR = 0
	FF = 0
	pe_sum_tlen = 0
	mp_sum_tlen = 0
	strand = ''
	insert_type = ''
	qnames = set()

# sam file의 각 line 읽기
	for line in myfile:
		line = line.rstrip()

# 빈 줄은 건너뛰기 (concat/edit으로 인한 blank line 방어)
		if not line:
			continue
# The header must be prior to the alignments. Header lines start with '@', while alignment lines do not.
		if line.startswith('@'):
			header_line += 1
			continue
# SAM은 TAB-delimited format이기 때문에 '\t'로 각 column을 split한 후 처리
# TLEN > 0인 leftmost alignment line만 사용하여 read pair를 한 번만 count한다. FLAG의 bit 값으로 핵심 정보 파악 후 strand orientation 판정.
		else:
			count_all_read += 1
			column = line.split('\t')
			qname = column[0]
			flag = int(column[1])
			pos = int(column[3])
			mapq = int(column[4])
			rnext = column[6]
			pnext = int(column[7])
			tlen = int(column[8])

			qnames.add(qname)
# tlen이 양수, 즉 leftmost read만 사용하여 RF/FR 판정
			if tlen <= 0:
				continue
# 0x100 : secondary alignment. 하나의 read가 여러 alignment line으로 표현되는 경우 제외. 
# properly-aligned read pairs = secondary/supplementary alignment 없이 제대로 align된 read pairs.
			if (flag & 0x100):
				continue
			if (flag & 0x800):
				continue
# 0x1 : template having multiple segments in sequencing. paired read인지 판정.
			if not (flag & 0x1):
				continue
# 0x2 : each segment properly aligned according to aligner. properly aligned paired read인지 판정.
			if not (flag & 0x2):
				continue
# 0x4 : segment unmapped. align 되지 않은 read인지 판정.
			if (flag & 0x4):
				continue
# 0x8 : next segment in the template unmapped.paired read가 align되지 않았는지 판정.
			if (flag & 0x8):
				continue
# 0x10 : SEQ being reverse complemented. reverse read인지 판정
# 0x20 : SEQ of the next segment in the template being reverse complemented. paired read가 reverse인지 판정
# paired read에는 PE와 MP만 존재한다고 가정. 따라서RR이나 FF시 Error 출력.
			if (flag & 0x10):
				strand = 'reverse'
				if (flag & 0x20):
					RR += 1
					'''print("Error(RR) in", column[0])'''
				else: # RF라면 MP
					MP += 1
					mp_sum_tlen += tlen

			if not (flag & 0x10):
				strand = 'forward'
				if not (flag & 0x20):
					FF += 1
					'''print("Error(FF) in", column[0])'''
				else:
					PE += 1
					pe_sum_tlen += tlen

if MP == 0:
	mp_mean_tlen = 0
else:
	mp_mean_tlen = (mp_sum_tlen) / (MP)

if PE == 0:
	pe_mean_tlen = 0
else:
	pe_mean_tlen = (pe_sum_tlen) / (PE)

if MP + PE == 0:
	mean_tlen = 0
else:
	mean_tlen = (mp_sum_tlen + pe_sum_tlen) / (MP+PE)

if MP > PE:
	insert_type = 'Mate Paired'
if PE > MP:
	insert_type = 'Paired End'
if PE == MP:
	insert_type = 'impossible to estimate.'

all_read_pairs = int(len(qnames))
# PE fraction = PE properly-aligned read pairs / all read pairs
# MP fraction = MP properly-aligned read pairs / all read pairs

pe_fraction = (PE / all_read_pairs)*100
mp_fraction = (MP / all_read_pairs)*100
rr_fraction = (RR / all_read_pairs)*100
ff_fraction = (FF / all_read_pairs)*100
improperly_aligned_fraction = 100 - pe_fraction - mp_fraction - rr_fraction - ff_fraction

fraction_outputs = [
	('Fraction for PE reads orientation', pe_fraction),
	('Fraction for MP reads orientation', mp_fraction)
	]
'''
	('Fraction for Error(RR) orientation', rr_fraction),
	('Fraction for Error(FF) orientation', ff_fraction),
	('Fraction for imporperly-aligned read pairs', improperly_aligned_fraction)
	]
'''

'''print('%-45s : %10i %1s' % ('# of header line', header_line, 'lines'))
print('%-45s : %10i %1s' % ('# of alignment line', count_all_read, 'lines'))
print('%-45s : %10i %1s' % ('# of this file', header_line+count_all_read, 'lines'))'''
print('%-45s : %10i %1s' % ('The number of all read pairs', all_read_pairs, 'pairs'))
print('%-45s : %10i %1s\n%-45s : %10i %1s' % ('# of properly-aligned MP', MP,'pairs','# of properly-aligned PE', PE, 'pairs'))
print('%-45s : %10.*f' % ('Mean TLEN of properly-aligned MP reads', 2, mp_mean_tlen))
print('%-45s : %10.*f' % ('Mean TLEN of properly-aligned PE reads', 2, pe_mean_tlen))
for i, o in fraction_outputs:
	print('%-45s : %10.2f%%' % (i, o))
print('%-45s : [%10s]' % ('The estimated insert type', insert_type))
print('%-45s : %10.*f' % ('Mean TLEN of this file', 2, mean_tlen))
