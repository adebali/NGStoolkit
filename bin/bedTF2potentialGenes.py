#!/usr/bin/env python
import bed
import sys
import argparse
import copy

parser = argparse.ArgumentParser(description='get two potential gene bodies based on the TFBS')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-l', type=int, default=1000, help='length of the potential gene body')
parser.add_argument('-g', required= True, help='genome (chromosome) sizes')

args = parser.parse_args()

genomeFile = args.g
chromosomeSizes = {}
for line in open(genomeFile, 'r'):
    ll = line.split('\t')
    chromosomeSizes[ll[0]] = int(ll[1])

bedFile = args.i
lines = []
for line in bedFile:
    lines.append(line)

out = args.o

bedObject = bed.bed(None, lines=lines)

flankingLength = args.l

for line in lines:
    bedLine = bed.bedline(line)
    upperLimit = chromosomeSizes[bedLine.chromosome()]
    lowerLimit = 1

    negativeEnd = bedLine.start()
    negativeStart = negativeEnd - flankingLength
    negativeFields =  copy.deepcopy(bedLine.fields())
    negativeFields[1] = negativeStart
    negativeFields[2] = negativeEnd
    negativeFields[5] = '-'
    
    if negativeStart > 0 and negativeEnd < upperLimit:
        out.write(bedLine.fields2line(negativeFields) + '\n')

    positiveStart = bedLine.end()
    positiveEnd = positiveStart + flankingLength
    positiveFields =  copy.deepcopy(bedLine.fields())
    positiveFields[1] = positiveStart
    positiveFields[2] = positiveEnd
    positiveFields[5] = '+'

    if positiveStart > 0 and positiveEnd < upperLimit:
        out.write(bedLine.fields2line(positiveFields) + '\n')
    

