#!/usr/bin/env python
import os
import sys
import generalUtils
import argparse

parser = argparse.ArgumentParser(description='normalize bed counts by dividing them by the maximum count number')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-t', required=  True, help='count column number')
parser.add_argument('-m', required=  False, help='multiplication factor')

args = parser.parse_args()

bedFile = args.i
outBedFile = args.o
countTabNumber = int(args.t)

if args.m:
	multiplyFactor = float(args.m)
else:
	multiplyFactor = 1

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')
maxCount = 0

for line in filein:
	ll = line.split('\t')
	count = int(ll[countTabNumber - 1])
	maxCount = max(maxCount, count)

filein = open(bedFile, 'r')

for line in filein:
	newLineList = line.strip().split('\t')
	newLineList[countTabNumber - 1] = str(round(float(newLineList[countTabNumber - 1]) * multiplyFactor / maxCount, 4))
	newLine = '\t'.join(newLineList) + '\n'
	out.write(newLine)

out.close()
