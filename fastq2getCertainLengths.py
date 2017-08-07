#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse

parser = argparse.ArgumentParser(description='retrieves sequences with the specified length(s)')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-l', nargs= '+', type=int, required= True, help='sequence length(s) of interest')

args = parser.parse_args()

fastqFile = args.i
out = open(args.o, 'w')
seqLengths = args.l

myFastq = fastq.fastq(fastqFile)
for seqObject in myFastq.stream():
    if len(seqObject.sequence) in seqLengths:
        out.write(seqObject.getAsString())
out.close()
