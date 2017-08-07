#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse
import random
import hashlib

parser = argparse.ArgumentParser(description='add md5sum of a sequence to the description line')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

input = args.i
out = open(args.o, 'w')

def sequence2md5sum(sequence):
	return hashlib.md5(sequence).hexdigest()

fastqObject = fastq.fastq(input)
for seqObject in fastqObject.stream():
	seqObject.assignHeader(seqObject.header + ' ' + sequence2md5sum(seqObject.sequence))
	out.write(seqObject.getAsString())
out.close()