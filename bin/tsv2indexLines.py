#!/usr/bin/env python

import argparse

PARSER = argparse.ArgumentParser(description='adds a column of index numbers to a tab separated file')
PARSER.add_argument('-i', required=True, help='input')
PARSER.add_argument('-o', required=True, help='output')
PARSER.add_argument('-start', required=True, type=int, default=1, help='index start number')
PARSER.add_argument('-end', required=False, default=False, help='index end number')

ARGS = PARSER.parse_args()
FILE = ARGS.i
OUT_FILE = ARGS.o
INDEX_START = ARGS.start
INDEX_END = ARGS.end
FILE_IN = open(FILE, 'r')
OUT = open(OUT_FILE, 'w')

index = INDEX_START - 1
for line in FILE_IN:
    if INDEX_END:
        if index == INDEX_END + 1:
            index = INDEX_START - 1
    index += 1
    ll = line.strip().split('\t')
    ll.append(str(index))
    out.write('\t'.join(ll))
