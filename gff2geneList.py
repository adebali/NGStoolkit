#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

# Uses chr as chromosome identifier.
parser = argparse.ArgumentParser(description='retrieves gene names from gff by order')
parser.add_argument('-i', required= True, help='input')

args = parser.parse_args()
gffFile = args.i

filein = open(gffFile, 'r')

for line in filein:
    if not line.startswith('#'):
        ll = line.split('\t')
        if ll[2] == 'gene':
            start = ll[3]
            end = ll[4]
            strand = ll[6]
            metaData = ll[8].split(';')
            gene = 'NaN'
            for info in metaData:
                if 'gene=' in info:
                    gene = info.split('gene=')[1].strip()
            print(gene)    

