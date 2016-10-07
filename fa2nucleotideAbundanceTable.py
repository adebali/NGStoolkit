#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='gets nucleotide distribution for each position')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-l', required=  False, help='sequence length of interest')
parser.add_argument('-n', required=  False, help='nucleotides of interest by order')
parser.add_argument('--percentage', action='store_true', help = 'Write percentages instead if actual counts')

args = parser.parse_args()

defaultNucleotideOrder = 'ATGC'

Fasta = fasta.fasta(args.i)
if args.l:
    nucleotideAbundanceDict = Fasta.getNucleotideAbundance(args.l)
else:
    nucleotideAbundanceDict = Fasta.getNucleotideAbundance()

if args.n:
    nucleotideOrder = args.n
else:
    nucleotideOrder = defaultNucleotideOrder

if args.percentage:
    percentageFlag = True
else:
    percentageFlag = False

Fasta.writeNucleotideAbundanceTable(nucleotideAbundanceDict, args.o, nucleotideOrder, percentageFlag)
