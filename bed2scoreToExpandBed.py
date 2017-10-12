#!/usr/bin/env python
import os
import sys
import argparse
from itertools import izip
import bed


parser = argparse.ArgumentParser(description='expanding bed files with scores')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')

args = parser.parse_args()

bedFile = args.i
out = open(args.o, 'w')

normalizationValue = 20000000

totalScore = 0
for bedline in bed.bed(bedFile).read():
	totalScore += float(min(0,bedline.score()))

for bedline in bed.bed(bedFile).read():
	expandedScore = max(0, normalizationValue * (float(bedline.score()) / totalScore))
	newScore = '.'
	fields = bedline.fields()
	fields[4] = newScore
	out.write(int(round(expandedScore)) * (bedline.fields2line(fields) + '\n'))
out.close()


	