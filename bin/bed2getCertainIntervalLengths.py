#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='retrieves certain fragment lengths')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-l', nargs='+', required= True, type=int, help='space-separeted sequence lengths of interest')

args = parser.parse_args()

filein = args.i
out = args.o
lengths = args.l

for line in filein:
    bedLine = bed.bedline(line)
    if bedLine.length() in lengths:
        out.write(line)
out.close()