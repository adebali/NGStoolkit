#!/usr/bin/env python
import sys
import os
import argparse
import bed

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-g', required= True, help='genome (chromosome) sizes')

args = parser.parse_args()
bedFile = args.i
genomeFile = args.g

chromosomeSizes = {}

for line in open(genomeFile, 'r'):
    ll = line.split('\t')
    chromosomeSizes[ll[0]] = int(ll[1])

bedFileObject = bed.bed(bedFile)
totalHitNumber = bedFileObject.getHitNum()
scaleFactor = 1000000 / float(totalHitNumber)

chromosomeCount = {}
for bedLine in bedFileObject.read():
    chromosome = bedLine.chromosome()
    chromosomeCount[chromosome] = chromosomeCount.get(chromosome, 0) + 1

normalizedChromosomeCount = {}
for chromosome in chromosomeCount.keys():
    normalizedChromosomeCount[chromosome] = scaleFactor * ((1000 * chromosomeCount[chromosome]) / float(chromosomeSizes[chromosome]))

for chromosome in normalizedChromosomeCount.keys():
    value = normalizedChromosomeCount[chromosome] 
    args.o.write(chromosome + '\t' + str(value) + '\n')