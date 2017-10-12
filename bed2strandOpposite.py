#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='changes the strand')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-s', required= False, type=int, default=6, help='column number of strand information')

args = parser.parse_args()

bedFile = args.i
output = args.o
out = open(output, 'w')
strandCol = args.s

d = {'+': '-', '-': '+', '.': '.'}
bed = bed.bed(bedFile)

for bedline in bed.read():
	strand = bedline.strand(strandCol)
	newStrand = d[strand]
	fields = bedline.fields()
	fields[strandCol - 1] = newStrand
	out.write(bedline.fields2line(fields) + '\n')
