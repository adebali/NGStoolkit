#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

parser = argparse.ArgumentParser(description='eliminates duplicate fragments')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-c', required=  True, help='columns that you want to be identical to define them as duplicate')

args = parser.parse_args()
bedFile = args.i
outBedFile = args.o
columns = args.c
prevLine = ''
filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
	if not generalUtils.linesAreDuplicate(line, prevLine, columns):
		out.write(line)
		prevLine = line
out.close()
