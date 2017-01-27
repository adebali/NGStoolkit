#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts given list to count table')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('--noInteger', required= False, action='store_true', help='no integer flag')
args = parser.parse_args()

filein = args.i
out = args.o
keysAreIntegers = not args.noInteger

def list2countTable(theList, keysAreIntegers = True):
    table = sorted(Counter(theList).items())
    return table

theList = []
for line in filein:
    if keysAreIntegers:
        e = int(line.strip())
    else:
        e = line.strip()
    theList.append(e)
table = list2countTable(theList)
for e in table:
    out.write(str(e[0]) + '\t' + str(e[1]) + '\n')

