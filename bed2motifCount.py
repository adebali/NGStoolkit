#!/usr/bin/env python

import bed
import argparse

parser = argparse.ArgumentParser(description='remove intervals that are more than N')
parser.add_argument('-bed', required= True, help='input bed file')
parser.add_argument('-fasta', required= True, help='input fasta file')
parser.add_argument('-o', required=  True, help='output')
parser.add_argument('-m', nargs="+", required=  True, help='motifs')
args = parser.parse_args()

motifs = args.m
out= open(args.o, 'w')

bedFile = bed.bed(args.bed)
for bedLine in bedFile.read():
    ll = bedLine.fields()
    sequence = bedLine.getFasta(args.fasta)
    for motif in motifs:
        ll.append(str(sequence.count(motif)))
    out.write("\t".join(ll))
    out.write("\n")
out.close()