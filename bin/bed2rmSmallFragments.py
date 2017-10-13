#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

parser = argparse.ArgumentParser(description='removes small fragments: eg small RNAs')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')
parser.add_argument('-c', required=  True, help='cutoff: remove fragments < cutoff')

args = parser.parse_args()
if args.i:
	bedFile = args.i
	filein = open(bedFile, 'r')
else:
	filein = sys.stdin
if args.o:
	outBedFile = args.o
	out = open(outBedFile, 'w')
else:
	out = sys.stdout

cutoff = int(args.c)

def line2removeDecision(line):
	ll = line.split('\t')
	start = int(ll[1])
	end = int(ll[2])
	length = end - start
	if length < cutoff:
		return True
	else:
		return False
for line in filein:
	removeDecision = line2removeDecision(line)
	if not removeDecision:
		out.write(line)
out.close()
