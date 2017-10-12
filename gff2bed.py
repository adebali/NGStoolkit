#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

# Uses chr as chromosome identifier.
parser = argparse.ArgumentParser(description='converts gff to bed')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-r', required= False, help='region')

args = parser.parse_args()
gffFile = args.i
outBedFile = args.o
if args.r:
    region = args.r
else:
    region = "NA"

filein = open(gffFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    if not line.startswith('#'):
        ll = line.split('\t')
        if ll[2] == region or region == "NA":
            chromosome = ll[0]
            start = ll[3]
            end = ll[4]
            strand = ll[6]
            name = "NA"
            if len(ll) >= 9:
                description = ll[8]
                if description.startswith("ID"):
                    name = description.split(';')[0].split('=')[1]
            out.write(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + name + '\t.\t' + strand + '\n')
out.close()
