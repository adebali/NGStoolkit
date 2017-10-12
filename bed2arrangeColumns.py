#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='get insert sizes of a bed file')
parser.add_argument('-i', required=True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-col', nargs="+", required=True, help='rule: $1 $2 $3 . 1 $4')
args = parser.parse_args()

bedFile = args.i
out = args.o
columns = args.col

bed = bed.bed(bedFile)
for bedline in bed.read():
	fields = bedline.fields()
	newFields = ['.'] * len(columns)
	for i in range(columns):
		value = columns[i]
		if '$' in value:
			newFields.append(fields[columns[i].split('$')[1] - 1])
		else:
			newFields.append(value)
	out.write(bedline.fields2line(newFields) + '\n')
out.close()
