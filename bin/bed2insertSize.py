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
args = parser.parse_args()

bedFile = args.i
out = args.o


bed = bed.bed(bedFile)
for bedline in bed.read():
	out.write(str(bedline.length()))
	out.write('\n')
out.close()
