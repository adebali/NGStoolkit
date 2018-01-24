#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedline

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed score to scoreless bed by expanding integered scores')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-placeHolder', type=int, default=0, help='number of placeholder dots (.)')


args = parser.parse_args()
bedFile = args.i
out = args.o

extraWord = ''
if args.placeHolder >=1:
    extraWord = '\t'


for line in bedFile:
    bedLine = bedline(line)
    score = int(bedLine.score())
    for i in range(score):
        out.write(bedLine.getLine() + extraWord + '\t'.join(args.placeHolder * ['.']) + '\n')