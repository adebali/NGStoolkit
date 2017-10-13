#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='count strings')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-fasta', required= True, help='fasta reference')
parser.add_argument('-string', required= True, help='string to count')

args = parser.parse_args()

bedFile = args.i
output = args.o
fasta = args.fasta
string = args.string

bedObject = bed.bed(bedFile)

for bedline in bedObject.read():
	count = bedline.countString(fasta, string)
	fields = bedline.fields()
	fields.append(str(count))
	line = bedline.fields2line(fields)
	out.write(line + '\n')