#!/usr/bin/env python
import os
import bed
import sys
import argparse
import generalUtils
import fasta
from sequence import DNA

parser = argparse.ArgumentParser(description='convert cufflinks output to bed file with counts')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-s', required= True, help='<Required> strand')
parser.add_argument('-o', required=  True, help='<Required> output')

args = parser.parse_args()

def cuffLinksLine2bedLine(line, strand):
    ll = line.split('\t')
    locus = ll[6]
    locusL = locus.split(':')
    chromosome = locusL[0]
    interval = locusL[1]
    intervalL = interval.split('-')
    start = intervalL[0]
    end = intervalL[1]
    FPKM = ll[9]
    LL = [chromosome, start, end, strand, FPKM]
    return '\t'.join(LL)

generalUtils.lineBasedFileOperation(args.i, args.o, cuffLinksLine2bedLine, [args.s])