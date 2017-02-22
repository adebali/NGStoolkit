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
args = parser.parse_args()

filein = args.i
out = args.o
# chr2 = int(args.chr2) -1
# s2 = chr2 + 1
# e2 = chr2 + 2

chr1 = 0
s1 = 1
e1 = 2

def line2position(line):
    ll = line.strip().split('\t')
    chr1name = ll[0].strip()
    chr2name = False
    for i in range(1,len(ll)-2):
        if ll[i].strip() == chr1name:
            chr2 = i
            break
    s2 = chr2 + 1
    e2 = chr2 + 2
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

##
#echo -e 'chr1\t10\t100\ta\tb\tc\td\tchr1\t10\t60\nchr1\t10\t100\ta\tb\tc\td\tchr1\t10\t80\n' | ./bedIntersect2parsedPosition.py -chr2 8
# -20
# -10
##