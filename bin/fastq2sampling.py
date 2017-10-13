#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='random sampling of reads')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-c', required= True, help='count')
parser.add_argument('-s', required= False, help='seed')

args = parser.parse_args()

fastqFile = args.i
output = args.o
out = open(output,'w')
count = int(args.c)
defaultSeed= 123

myFastq = fastq.fastq(fastqFile)
sequenceNumber = myFastq.getSequenceCount()

if args.s:
	random.seed(args.s)
else:
	random.seed(defaultSeed)

randomSequenceNumbers = sorted(random.sample(range(1, sequenceNumber), count))
i = 0
with open(fastqFile) as input:
	for line1 in input:
		i += 1
		line2 = input.next()
		line3 = input.next()
		line4 = input.next()
		if i == randomSequenceNumbers[0]:
			# print(i)
			out.write(line1)
			out.write(line2)
			out.write(line3)
			out.write(line4)
			del(randomSequenceNumbers[0])
			if len(randomSequenceNumbers) == 0:
				break
out.close()