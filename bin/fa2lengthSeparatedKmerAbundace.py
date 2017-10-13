#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='gets kmer (eg. dimer) distribution for each position')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-k', required= True, help='k of kmer')
parser.add_argument('-r', required= True, help='length range (eg 9-13)')

args = parser.parse_args()
input = args.i
output = args.o
kmer = int(args.k)
lengthRange = args.r

fasta.fasta(input).separateByLengthAndWriteKmerAbundance(kmer, lengthRange, output)