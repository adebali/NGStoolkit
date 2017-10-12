#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='get twoo flanking regions based on the average length of the main intervals')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-g', required= True, help='genome (chromosome) sizes')
# parser.add_argument('--average', required=False, type=bool, default=True, help='Flag if we want to remove the regions at the average length of the edges (default=True)')
parser.add_argument('--fixed', action='store_true', default=False, help='fixed region length flanking. It overrides --average')
parser.add_argument('-l', required=False, type=int, default=1000, help='length of the flanking regions (should be used with --fixed, default=1000)')

args = parser.parse_args()

bedFile = args.i
out = args.o
genomeFile = args.g

chromosomeSizes = {}

for line in open(genomeFile, 'r'):
    ll = line.split('\t')
    chromosomeSizes[ll[0]] = int(ll[1])

lines = []
for line in bedFile:
    lines.append(line)

bedObject = bed.bed(None, lines=lines)

if args.fixed:
    lengthOfRegion = args.l
else:
    averageLength = round(bedObject.getAverageLength())
    lengthOfRegion = averageLength

for line in lines:
    bedLine = bed.bedline(line)
    if bedLine.start() - lengthOfRegion > 0 and bedLine.end() + lengthOfRegion < chromosomeSizes[bedLine.chromosome()]:
        out.write(bedLine.getLine() + '\n')
