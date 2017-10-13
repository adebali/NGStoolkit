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
parser.add_argument('-original', nargs='+', required= True, help='space separated original chromosome names')
parser.add_argument('-new', nargs='+', required= True, help='space separated new chromosome names')

args = parser.parse_args()

bedFile = args.i
output = args.o
out = open(output, 'w')
originalChromosomes = args.original
newChromosomes = args.new

d = {}
i = 0
for chr in originalChromosomes:
	d[chr] = newChromosomes[i]
	i += 1

bed = bed.bed(bedFile)

for bedline in bed.read():
	fields = bedline.fields()
	oldChromosome = bedline.chromosome()
	newChromosome = d[oldChromosome]
	fields[0] = newChromosome
	out.write(bedline.fields2line(fields) + '\n')
