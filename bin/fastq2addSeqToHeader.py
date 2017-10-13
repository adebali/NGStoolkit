#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse
import random
import hashlib

parser = argparse.ArgumentParser(description='add sequence to the header line')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

input = args.i
out = open(args.o, 'w')


fastqObject = fastq.fastq(input)
for seqObject in fastqObject.stream():
	seqObject.assignHeader(seqObject.header + ' ' + seqObject.sequence)
	out.write(seqObject.getAsString())
out.close()