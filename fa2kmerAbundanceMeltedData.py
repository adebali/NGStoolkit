#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='gets kmer (eg. dinucleotide) distribution for each position')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-k', required= True, help='k of kmer')
parser.add_argument('-l', required=  False, help='sequence length of interest')
parser.add_argument('-c', nargs = '+', required=  False, help='class dictionary: space-separated comma-separated pairs. Eg: product,CPD replicate,A')
parser.add_argument('--percentage', action='store_true', help = 'Write percentages instead if actual counts')

args = parser.parse_args()

dictionary = {}
if args.c:
    for pair in args.c:
        ll = pair.split(',')
        dictionary[ll[0]] = ll[1]
if args.percentage:
    percentageFlag = True
else:
    percentageFlag = False

Fasta = fasta.fasta(args.i)
if args.l:
    kmerAbundanceDict = Fasta.getKmerAbundance(int(args.k), args.l)
else:
    kmerAbundanceDict = Fasta.getKmerAbundance(int(args.k))

Fasta.writeKmerAbundanceMeltedData(kmerAbundanceDict, args.o, dictionary, percentageFlag)
