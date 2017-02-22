#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

parser = argparse.ArgumentParser(description='converts sga to bed')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-g', required=  True, help='genome chromosome limits')
parser.add_argument('-b', required=  True, help='slop b: window size from both sides')

args = parser.parse_args()
infile = args.i
out = args.o
chromosomeSizesFile = args.g
window = int(args.b)

def getChromosomeSizesDict(chromosomeSizesFile):
	myDict = {}
	for line in open(chromosomeSizesFile):
		ll = line.strip().split('\t')
		myDict[ll[0].strip()] = int(ll[1].strip())
	return myDict

def point2interval(chromosome, point, window, chromosomeSizes):
	start = min(max(0, point - window), chromosomeSizes[chromosome])
	end = min(max(0, point + window), chromosomeSizes[chromosome])
	if start == 0 or end == chromosomeSizes[chromosome]:
		return None
	else:
		return [start, end]

def sgaLine2bedLine(line, window, chromosomeSizes):
	ll = line.strip().split('\t')
	point = int(ll[2])
	interval = point2interval(ll[0], point, window, chromosomeSizes)
	if interval == None:
		return None
	start = interval[0]
	end = interval[1]
	newLine = '\t'.join([ll[0], str(start), str(end), ll[3], ll[4]])
	return newLine

chromosomeSizes = getChromosomeSizesDict(chromosomeSizesFile)
#print(chromosomeSizes)


for line in infile:
	newLine = sgaLine2bedLine(line, window, chromosomeSizes)
	if newLine:
		out.write(newLine + '\n')


# generalUtils.lineBasedFileOperation(sgaFile, outBedFile, sgaLine2bedLine, [window, chromosomeSizes])

