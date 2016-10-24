#!/usr/bin/env python
# DOESN'T WORK!!!!!!!!!!!!
import os
import sys
import fastq
import generalUtils
import argparse

parser = argparse.ArgumentParser(description='writes sequence length distribution')
parser.add_argument('-inputs', nargs='+', help='<Required> Input file list (min 2)', required=True)
parser.add_argument('-o', required= True, help='<Required> Output')
parser.add_argument('-same', nargs='+', required= True, help='<Required> identical (same) column numbers (eg. 1 3 for identical first and third columns)')
parser.add_argument('-merge', nargs='+', required= True, help='<Required> to be merged column no list')


args = parser.parse_args()
inputs = arg.inputs
output = args.o
identicalColumns = args.same
mergingColumns = args.merge

separator = '\t'

out = open(output, 'w')
firstFileIdenticalColumnList = identicalColumnList = []
fileNumber = 1

for file in inputs:
    i = 0
    filein = open(file, 'r')
    for line in filein:
        ll = line.strip().split(separator)
        for identicalColumn in identicalColumns:
            identicalColumnList[i] = ll[int(identicalColumn)] + separator
        if fileNumber != 1 and identicalColumnList[i] != firstFileIdenticalColumnList[i]:
            sys.exit('Identical Columns do not match between files')
        elif fileNumer == 1:
            newLine = identicalColumnList[i]
        for columnNo in mergingColumns:
            value = ll[int(columnNo)]
            newLine += value + separator
        i += 1
    
    if fileNumber == 1:
        firstFileIdenticalColumnList = list(identicalColumnList)

    fileNumber += 1