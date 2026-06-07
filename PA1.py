#!/usr/bin/env python3

import sys

# 실행 파일 / fastq 파일 총 2개를 argument로 받아야만 실행하기
if len(sys.argv) != 2:
    print('Usage: ./PA1.py <fastq_file>')

else:
    fastq_file = sys.argv[1]

# read sequence의 수 / total length / distribution 처리할 dictionary
    read_count = 0
    total_length = 0
    length_dict = {}

# line 단위로 읽기 위한 준비. header-sequence-optional-quality 순서로 등장한다는 가정.
    myfile = open(fastq_file)

# state는 line과 line 사이 관계를 밝히는 object. seq은 저장할 각 read sequence. qual은 seq이 정확히 읽혔는지 검증하고 다음 line으로 넘어가기 위한 object.
    state = 'header'
    seq = ''
    qual = ''

# fastq 파일의 각 line 읽기
    for line in myfile:
        line = line.rstrip()

# state가 'header'일 경우, @로 시작하는지 검증 후 seq/qual 초기화 후 state를 'seq'으로 변경 후 다음 line 읽기
        if state == 'header':
            if line.startswith('@'):
                seq = ''
                qual = ''
                state = 'seq'
                continue
# bad case의 경우 보고 후 프로그램 종료
            else:
                print('Error: invalid FASTQ format. Check Sequence identifier line.')
                break

# state가 'seq'일 경우, seq에 line에 저장된 sequence를 string으로 추가 후 다음 line 읽기. +로 시작하는 optional line일 경우 state를 'qual'로 변경 후 다음 line 읽기
        elif state == 'seq':
            if line.startswith('+'):
                state = 'qual'
                continue
            else:
                seq = seq + line
                continue

# state가 'qual'일 경우, line에 저장된 sequence를 qual에 string으로 추가 후 seq과 길이가 동일한지 검증. 
        elif state == 'qual':
            qual = qual + line
# 길이가 동일하다면 read_length에 seq의 길이 저장 후 read_count를 하나 센 다음, total_length에 read_length의 값 추가.
            if len(qual) == len(seq):
                read_length = len(seq)
                read_count = read_count + 1
                total_length = total_length + read_length
# length_dict에 이미 있는 read_length라면, 저장돼 있는 key에 count만 하나 추가.
                if read_length in length_dict:
                    length_dict[read_length] = length_dict[read_length] + 1
# length_dict에 없는 read_length라면, key를 하나 만들고 value도 1 부여.
                else:
                    length_dict[read_length] = 1
# length_dict 업데이트 후 state를 header로 변경 후 동일한 과정 반복.
                state = 'header'
# bad case의 경우 보고 후 프로그램 종료
            elif len(qual) != len(seq):
                print('Error: invalid FASTQ format. Check Quality score line.')
                break

    myfile.close()
# The number / Total length of read sequences
    print('The number of read sequences is %d' % read_count)
    print('Total length of read sequences is %d' % total_length)
    print()
# fastq 파일의 고유한 length 값들을 리스트 lengths에 저장
    lengths = list(length_dict.keys())
# 길이 순서대로 distribution을 나타낼 것이기 때문에 sort method 사용.
    lengths.sort()

# Distribution
    print('<Text distribution of read lengths>')
    print('%8s%14s%12s' % ('Length', 'Count', 'Abundance'))
    
    for length in lengths:
        abundance = (length_dict[length] / read_count) * 100
        print('%8d%14d%11.2f%%' % (length, length_dict[length], abundance))
    print()

    print('<Graphic distribution of read lengths>')
    print('%8s%14s%12s' % ('Length', 'Distribution', 'Abundance'))

    for length in lengths:
        count = length_dict[length]
        abundance = (count / read_count) * 100
# Abundance(%)를 10% 단위로 시각화
        black = int(abundance // 10)
        

        distribution = chr(0x25A0) * black + chr(0x25A1) * (10 - black)

        print('%8d%14s%11.2f%%' % (length, distribution, abundance))

