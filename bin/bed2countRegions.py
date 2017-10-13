#!/usr/bin/env python
import os
import sys
import argparse
import tempfile

parser = argparse.ArgumentParser(description='count strings. Assumes bed6 as input.')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')

args = parser.parse_args()

bedFile = open(args.i, 'r')
output = args.o
_temp, temp = tempfile.mkstemp()

previousLine = ''
count = 1
totalCount = 0


for line in bedFile:
	if len(line.split('\t')) != 6:
		raise ValueError("Expected bed6 but got bed" + str(len(line.split('\t'))))
	totalCount += 1
	if line == previousLine:
		count += 1
	else:
		if previousLine == '':
			previousLine = line
		os.write(_temp, previousLine.strip() + '\t' + str(count) + '\n')
		count = 1
	previousLine = line
else:
	if line == previousLine:
		count += 1
	else:
		count = 1
	os.write(_temp, previousLine.strip() + '\t' + str(count) + '\n')

os.close(_temp)
print(temp)

out = open(output, 'w')
for line in open(temp, 'r'):
	ll = line.strip().split('\t')
	count = int(ll[6]) * (float(1000000)/totalCount)
	del(ll[6])
	ll[4] = str(count)
	newLine = '\t'.join(ll) + '\n'
	out.write(newLine)
out.close()
# os.close(_temp)
# os.remove(temp)