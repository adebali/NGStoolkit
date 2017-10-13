#!/usr/bin/env python
import os
import sys
import fastq
import generalUtils
import argparse

parser = argparse.ArgumentParser(description='writes sequence length distribution')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-l', required= True, help='sequence length of interest')

args = parser.parse_args()

fastqFile = args.i
output = args.o
seqLength = int(args.l)

myFastq = fastq.fastq(fastqFile)
lengthDistro = myFastq.getLengthDistribution(seqLength)
generalUtils.writeDict(lengthDistro, output)
