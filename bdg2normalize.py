#!/usr/bin/env python
import os
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='fix fragment length if possible')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-o', required=  True, help='<Required> output')
parser.add_argument('-c', required= False, default=4, type=int, help='countTab')
parser.add_argument('-x', required=  False, default=1, type=float, help='Scale Factor')

args = parser.parse_args()
bedFile =  args.i
outBedFile = args.o
countTabNumber = args.c
scaleFactor = args.x

out = open(outBedFile, 'w')

for bedLine in bed.bed(bedFile).read():
    count = float(bedLine.fields()[countTabNumber - 1])
    newCount = str(count * float(scaleFactor))
    newLine = bedLine.newline(bedLine.start(), bedLine.end(), {4: newCount})
    out.write(newLine + '\n')
out.close()
