#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedline

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed to one side bed (start or end)')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-s', default='start', choices=['start', 'end'], help='side (start or end)')
parser.add_argument('-l', type=int, default=1000, help='length of each side')

args = parser.parse_args()
bedFile = args.i
out = args.o
sideLength = int(args.l)
side = args.s

for line in bedFile:
    bedLine = bedline(line)
    if side == 'start':
        referencePoint = bedLine.start()
    else:
        referencePoint = bedLine.end()
    start = referencePoint - sideLength
    end = referencePoint + sideLength
    if start > 0:
        out.write(bedLine.newline(start, end) + '\n')