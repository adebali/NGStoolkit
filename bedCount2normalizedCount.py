#!/usr/bin/env python
import os
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='fix fragment length if possible')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-c', required= True, help='<Required> countTab')
parser.add_argument('-o', required=  True, help='<Required> output')
parser.add_argument('-l', required=  False, help='per sequence length, default = 1000')
parser.add_argument('--bypassLength', required= False, action= "store_true", default= False, help='per sequence length, default = 1000')
parser.add_argument('-m', required=  False, help='multiplyFactor, default = 1')

args = parser.parse_args()
bedFile =  args.i
outBedFile = args.o
if args.l:
    perSequenceLength = int(args.l) 
else:
    perSequenceLength = 1000

if args.m:
    multiplyFactor = float(args.m) 
else:
    multiplyFactor = 1


countTabNumber = int(args.c)
begTab = 2
endTab = 3

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    ll = line.strip().split('\t')
    beg = int(ll[begTab - 1])
    end = int(ll[endTab - 1])
    count = int(ll[countTabNumber - 1])

    if args.bypassLength:
        regionLength = perSequenceLength
    else:
        regionLength = end - beg
    newCount = (float(count) / regionLength) * float(perSequenceLength) * multiplyFactor
    

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
