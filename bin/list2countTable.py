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
    table = {}
    for e in theList:
        table[e] = table.get(e, 0) + 1 
    # table = sorted(Counter(theList).items())
    return table


i = 0
theList = []
for line in filein:
    i += 1
    if i % 1000000 == 0:
        print(i)
    if keysAreIntegers:
        e = int(line.strip())
    else:
        e = line.strip()
    theList.append(e)
table = list2countTable(theList)
for key in sorted(table.keys()):
    out.write(str(key) + '\t' + str(table[key]) + '\n')

