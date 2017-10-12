#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='replace score column with a constant value')
parser.add_argument('-i', required=True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-s', required=False, default=".", help='score')
args = parser.parse_args()

bedFile = args.i
out = args.o


bed = bed.bed(bedFile)
for bedline in bed.read():
	fields = bedline.fields()
	fields[4] = args.s
	out.write(bedline.fields2line(fields) + '\n')
out.close()
