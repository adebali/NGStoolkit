#!/usr/bin/env python
import sys
import os
import argparse
import bed
import generalUtils

parser = argparse.ArgumentParser(description='takes bed as input, get the middle point and extend it to both sides')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-g', required= True, help='genomeFile')
parser.add_argument('-w', required= True, help='windowSize')
args = parser.parse_args()
bedFile = args.i
output = args.o
genome = args.g
windowSize = int(args.w)

chrDict = {}

for line in open(genome, 'r'):
    ll = line.split('\t')
    chrDict[ll[0]] = int(ll[1])

def line2newLine(line):
    bedLine = bed.bedline(line)
    chromosome = bedLine.chromosome()
    start = bedLine.start()
    newEnd = min(start + windowSize, chrDict[chromosome])
    return bedLine.newline(start, newEnd)

generalUtils.lineBasedFileOperation(bedFile, output, line2newLine, [])
