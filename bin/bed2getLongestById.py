#!/usr/bin/env python
import bed
import sys
import argparse
from bed import bed


parser = argparse.ArgumentParser(description='name section must be filled. This script is intended to remove isoforms and keep the longest one. It benefits from the ID (name) value')
parser.add_argument('-i', required=True, help='input')
parser.add_argument('-o', required=False, nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')



args = parser.parse_args()

inputFile = args.i
out = args.o

bedObject = bed(inputFile)

d = {}

def getLongest(bedlineList):
    selectedLine = bedlineList[0]
    for bedline in bedlineList:
        if bedline.length() > selectedLine.length():
            selectedLine = bedline
        elif bedline.length() == selectedLine.length():
            if bedline.score() > selectedLine.score():
                selectedLine = bedline
    return selectedLine

for bedline in bedObject.read():
    name = bedline.name()
    parent = name.split('.')[0]
    if parent not in d.keys():
        d[parent] = []
    d[parent].append(bedline)

for parent in d.keys():
    bedlineList = d[parent]
    longestBedline = getLongest(bedlineList)
    out.write(longestBedline.getLine() + '\n')