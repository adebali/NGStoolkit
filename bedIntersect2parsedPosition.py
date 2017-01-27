#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed intersect result to position list. Format should be \
chr1 start1 end1 chr1 start2 end2 anythingElse ...')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-s2', required=False, default=5, help='b bed start column no')
parser.add_argument('-e2', required=False, default=6, help='b bed end column no')
args = parser.parse_args()

filein = args.i
out = args.o
s2 = int(args.s2) - 1 
e2 = int(args.e2) - 1 
chr2 = s2 - 1

chr1 = 0
s1 = 1
e1 = 2

def line2position(line):
    ll = line.strip().split('\t')
    chromosomeA = ll[chr1]
    startA = int(ll[s1])
    endA = int(ll[e1])
    midA = (startA + endA) / 2
    startB = int(ll[s2])
    endB = int(ll[e2])
    midB = (startB + endB) / 2
    damgePosition = midB - midA
    return damgePosition

for line in filein:
    if line.strip():
        out.write(str(line2position(line.strip())) + '\n')

