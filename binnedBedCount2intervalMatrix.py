#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils
import bed

if not '--test' in sys.argv:
	parser = argparse.ArgumentParser(description='seperate intervals into bins')
	parser.add_argument('-i', required= True, help='input')
	parser.add_argument('-o', required= True, help='output')
	parser.add_argument('-c', required=  True, help='column of interval identifier')
	args = parser.parse_args()

	bedFile = args.i
	outBedFile = args.o
	intervalIdColumn = int(args.c)
else:
	bedFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utils', 'testFiles', 'binnedBedExample_count.bed')
	outBedFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utils', 'testFiles', 'binnedBedExample_count_clustered.bed')
	intervalIdColumn = 7

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

intervalNumber = 0
previousIntervalId = None
newLL = []
for line in filein:
	ll = line.strip().split('\t')
	intervalId = ll[intervalIdColumn - 1]
	if intervalId != previousIntervalId:
		previousIntervalId = intervalId
		if newLL:
			out.write('\t'.join(newLL) + '\n')
		newLL = ll[:intervalIdColumn]
		newLL += [ll[intervalIdColumn]]
	else:
		newLL += [ll[intervalIdColumn]]
else:
	out.write('\t'.join(newLL) + '\n')
out.close()