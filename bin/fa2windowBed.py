#!/usr/bin/env python
import os
import argparse

# Requirements
# bedtools makewindows

parser = argparse.ArgumentParser(description='prints sequence header and lengths')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-n', required= True, help='window size')

args = parser.parse_args()
input = args.i
output = args.o
windowSize = args.n
seqDict = {}

filein = open(input, 'r')
header = False

seq = ''

for line in filein:
	if line.startswith('>'):
		if header:
			seqDict[header] = len(seq)
		header = line[1:].strip()
		seq = ''
		seqLength  = 0
	else:
		seq += line.strip()
		# seqLength = len(seq)
		# seqLength += seqLength
else:
	seqDict[header] = len(seq)

temp = open('temp.txt', 'w')

for header in seqDict.keys():
	temp.write(header + '\t' + str(seqDict[header]) + '\n')
temp.close()

os.system('bedtools makewindows -g temp.txt -w ' + str(windowSize) + ' >' + output)
# os.system('rm temp.txt')
