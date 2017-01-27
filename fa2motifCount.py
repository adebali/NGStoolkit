#!/usr/bin/env python

import fasta
import argparse
import re

parser = argparse.ArgumentParser(description='gets total motif count')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-m', required= True, help='regex pattern for motif of interest')

args = parser.parse_args()
input = args.i
motif = args.m

totalCount = 0
totalLength = 0
fastaObject = fasta.fasta(args.i)
for seqObject in fastaObject.stream(500000):
    sequence = seqObject['s']
    motifCount = len(re.findall(motif, sequence))
    totalCount += motifCount
    totalLength += len(sequence)

print(totalCount * 1000 / float(totalLength))