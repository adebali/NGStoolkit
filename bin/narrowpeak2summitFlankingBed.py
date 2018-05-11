#!/usr/bin/env python
import sys
import os
import argparse
import bed
import generalUtils

parser = argparse.ArgumentParser(description='takes bed as input, get the middle point and extend it to both sides')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-w', required= True, help='windowSize')
parser.add_argument('-g', required= False, default=False, help='genomeFile')
parser.add_argument('--randomMid', required= False, action='store_true', help='for cases of .5 middle point, randomly select between positions 0 or 1')
args = parser.parse_args()
bedFile = args.i
output = args.o
windowSize = int(args.w)

if args.g:
    chromosomeSizes = {}
    for line in open(args.g, 'r'):
        ll = line.split('\t')
        chromosomeSizes[ll[0]] = int(ll[1])

def getInterval(line, randomness=False):
    bedLine = bed.bedline(line)
    middlePoint = bedLine.start() + int(bedLine.getField(10))
    start = middlePoint - windowSize
    end = middlePoint + windowSize
    if args.g:
        chromosome = bedLine.chromosome()
        chrEnd = chromosomeSizes[chromosome]
        if start > 0 and end < chrEnd:
            return bedLine.newline(start, end)
        return False
    return bedLine.newline(start, end)
    

generalUtils.lineBasedFileOperation(bedFile, output, getInterval, [])
