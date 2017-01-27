#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='convert Fastq to tab separated single line per sequence Fastq file | or reverse process')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('--reverse', required= False, action='store_true', help='output')

args = parser.parse_args()

input = args.i
output = args.o
reverseFlag = args.reverse

def singleFastqLineToFourLines(input, output):
	filein = open(input, 'r')
	out = open(output, 'w')
	for line in filein:
		ll = line.split('\t')
		for e in ll:
			out.write(e.strip() + '\n')
	out.close()

if not reverseFlag:
	myFastq = fastq.fastq(input)
	myFastq.writeSingleLineFastq(output)
else:
	singleFastqLineToFourLines(input, output)