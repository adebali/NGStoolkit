#!/usr/bin/env python
import os
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='normalize bed counts by read number')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-c', required= True, help='<Required> countTab')
parser.add_argument('-readNum', required= True, help='<Required> total mapped read number of the run')
parser.add_argument('-perNumReads', required=  False, help='per number reads value, default = 1000000')
parser.add_argument('-o', required=  True, help='<Required> output')

args = parser.parse_args()
bedFile =  args.i
outBedFile = args.o
if args.perNumReads:
    perNumReads = int(args.perNumReads) 
else:
    perNumReads = 1000000

readNum = int(args.readNum)

countTabNumber = int(args.c)

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    ll = line.strip().split('\t')
    count = float(ll[countTabNumber - 1])
    newCount = (float(count) / readNum) * float(perNumReads)

    newLine = ''
    for i in range(len(ll)):
        if i == countTabNumber - 1:
            field = str(newCount)
        else:
            field = ll[i]
        newLine += field + '\t'
    newLine = newLine.strip() + '\n'
    out.write(newLine)

out.close()
