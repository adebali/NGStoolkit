#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='random sampling of reads')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-c', required= True, help='count')
parser.add_argument('-s', required= False, help='seed')

args = parser.parse_args()

bedFile = args.i
output = args.o
out = open(output,'w')
count = int(args.c)
defaultSeed= 123

bed = bed.bed(bedFile)
hitNumber = bed.getHitNum()

if args.s:
	random.seed(args.s)
else:
	random.seed(defaultSeed)
print(hitNumber)
print(count)
randomHitNumbers = sorted(random.sample(range(1, hitNumber+1), min(count, hitNumber)))
i = 0
with open(bedFile) as input:
	for line in input:
		i += 1
		if i == randomHitNumbers[0]:
			out.write(line)
			del(randomHitNumbers[0])
			if len(randomHitNumbers) == 0:
				break
out.close()

