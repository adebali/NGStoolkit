#! /usr/bin/env python
from itertools import izip
import sys
import argparse

parser = argparse.ArgumentParser(description='categorizes TS and NTS based on the strand information')
parser.add_argument('-plus', required= True, help='plus strand input')
parser.add_argument('-minus', required= True, help='minus strand input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-s', required= False, type=int, default=6, help='strand column')

args = parser.parse_args()
plusInput =  args.plus
minusInput =  args.minus
out = open(args.o, 'w')
strandColumn = args.s - 1

with open(plusInput) as plusFile, open(minusInput) as minusFile:
    for p, m in izip(plusFile, minusFile):
        pl = p.strip().split("\t")
        ml = m.strip().split("\t")
        if pl[0:-1] != ml[0:-1]:
            raise ValueError("Plus and Minus files do not have identical line starts.\n" + p + "\n" + m)
        pval = pl[-1]
        mval = ml[-1]
        strand = pl[strandColumn]
        if strand == "+":
            TS = mval
            NTS = pval
        elif strand == "-":
            TS = pval
            NTS = mval
        else:
            raise ValueError("strand is neither + nor -")
        newLineList = pl[:-1] + [TS, NTS]
        out.write("\t".join(newLineList) + "\n")
        
