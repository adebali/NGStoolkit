#!/usr/bin/env python
import sys
import bed
import argparse

parser = argparse.ArgumentParser(description='converts bed file into smaller bins')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-u', required=True, help='upstream length')
parser.add_argument('-d', required=True, help='downstream length')
parser.add_argument('-g', required=True, help='genome chromosome sizes')
parser.add_argument('-s', required=True, help='transcription site: start | end')
parser.add_argument('-strandColumn', required=False, default=6, help='strand column')
args = parser.parse_args()

filein = args.i
out = args.o
genomeChrSizesFile = args.g
upstreamLength = int(args.u)
downstreamLength = int(args.d)
site = args.s
if site != 'start' and site != 'end':
    raise ValueError('site is not given as expected (start or end): ' + site)
strandColumn = int(args.strandColumn)

def getGenomeChrSizes(genomeChrSizesFile):
    myDict = {}
    filein = open(genomeChrSizesFile, 'r')
    for line in filein:
        ll = line.strip().split('\t')
        myDict[ll[0]] = int(ll[1])
    return myDict

for line in filein:
    chrSizes = getGenomeChrSizes(genomeChrSizesFile)
    bedLine = bed.bedline(line)
    chromosome = bedLine.chromosome()
    strand = bedLine.strand(strandColumn)
    if strand == '+':
        if site == 'start':
            start = bedLine.start()
        elif site == 'end':
            start = bedLine.end()
        upstream = start - upstreamLength
        downstream = start + downstreamLength
        newStart = upstream
        newEnd = downstream
    elif strand == '-':
        if site == 'start':
            start = bedLine.end()
        elif site == 'end':
            start = bedLine.start()
        upstream = start + upstreamLength
        downstream = start - downstreamLength
        newStart = downstream
        newEnd = upstream
    else:
        raise ValueError('strand is not as expected: ' + strand)
    if newStart > 0 and newEnd < chrSizes[chromosome]:
        out.write(bedLine.newline(newStart, newEnd) + '\n')