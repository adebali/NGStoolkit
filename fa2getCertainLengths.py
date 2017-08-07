#!/usr/bin/env python

import fasta
import argparse
import re
import sys

parser = argparse.ArgumentParser(description='returns fasta file having sequence with the desired length(s)')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')
parser.add_argument('-l', nargs='+', type=int, required= True, help='sequence length(s) of interest')

args = parser.parse_args()

fastaObject = fasta.fasta(args.i)
output = args.o
sequenceLengthList = args.l
out = open(output, 'w')
seqDicts = fastaObject.stream2(100*4096)
for seqObject in seqDicts:
    if seqObject.getLength() in sequenceLengthList:
        out.write('>' + seqObject.getHeader() + '\n' + seqObject.getSequence() + '\n')
out.close()
