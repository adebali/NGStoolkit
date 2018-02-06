#!/usr/bin/env python

import argparse
import os
import sys
import tempfile

PARSER = argparse.ArgumentParser(description='removes neighboring intervals')
PARSER.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdout, help='input')
PARSER.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
PARSER.add_argument('-wnum', required=True, type=int, help='window number')
PARSER.add_argument('-countCol', required=False, default=7, type=int, help='count column')

args = PARSER.parse_args()
wnum = args.wnum
countCol = args.countCol
d = {}
i = 0
lineNum = 0
for line in args.i:
    i += 1
    position = i%wnum
    if position == 0:
        lineNum += 1
        position = wnum
    ll = line.strip().split('\t')
    d[position] = d.get(position, 0) + float(ll[countCol - 1])

for i in range(1, wnum + 1):
    args.o.write(str(i) + '\t' + str(d[i]/float(lineNum)) + '\n')