#!/usr/bin/env python
import bed
import sys
import argparse
from pybedtools import BedTool
from bed import bedline


parser = argparse.ArgumentParser(description='retrieves certain fragment lengths')
parser.add_argument('-b', reqiured=True, help='input')
parser.add_argument('-o', reqiured=True, nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-a', reqiured=True, help='interval file')
parser.add_argument('-fasta', required=False, default=False, help='fasta file')


args = parser.parse_args()

inputFile = args.b
out = args.o
bedFile = args.a

a = BedTool(inputFile)
b = BedTool(bedFile)

a_and_b = a.intersect(b, c=True)
for line in a_and_b:
    bedLine = bedline(str(line))
    TTcount = bedLine.getFasta().count('TT')
    out.write(bedLine.addField([str(TTcount)]))

# out.write(str(a_and_b))