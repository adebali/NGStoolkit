#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

# Uses chr as chromosome identifier.
parser = argparse.ArgumentParser(description='converts gff to bed')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-r', required= True, help='region')

args = parser.parse_args()
gffFile = args.i
outBedFile = args.o
region = args.r

filein = open(gffFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    if not line.startswith('#'):
        ll = line.split('\t')
        if ll[2] == region:
            start = ll[3]
            end = ll[4]
            strand = ll[6]
            out.write('chr' + '\t' + str(start) + '\t' + str(end) + '\t' + strand + '\n')
out.close()
