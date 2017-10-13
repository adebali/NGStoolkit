#!/usr/bin/env python
import fasta
import argparse

parser = argparse.ArgumentParser(description='prints sequence header and lengths')
parser.add_argument('-i', required= True, help='input')

args = parser.parse_args()
Fasta = fasta.fasta(args.i)
Fasta.getSequenceLengths()

