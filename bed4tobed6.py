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
parser.add_argument('-s', required=False, default=".", help='score')
parser.add_argument('-n', required=False, default=".", help='name')
args = parser.parse_args()

bedFile = args.i
out = args.o


bed = bed.bed(bedFile)
for bedline in bed.read():
	fields = bedline.fields()
	fields.append(args.n)
	fields.append(args.s)
	fields[5] = fields[3]
	fields[3] = '.'
	out.write(bedline.fields2line(fields) + '\n')
out.close()
