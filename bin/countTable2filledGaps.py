#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='fill gaps with "zero" count for the missing keys')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-s', required= False, help='start')
parser.add_argument('-e', required= False, help='end')

args = parser.parse_args()

filein = args.i
out = args.o

myDict = {}
newDict = {}
for line in filein:
    ll = line.strip().split('\t')
    key = int(ll[0])
    value = int(ll[1])
    myDict[key] = value


if args.s:
    start = int(args.s)
else:
    start = min(myDict.keys())
if args.e:
    end = int(args.e)
else:
    end = max(myDict.keys())
    

for i in range(start, end + 1):
    if i in myDict.keys():
        value = myDict[i]
    else:
        value = 0
    newDict[i] = value

for i in sorted(newDict.keys()):
    out.write(str(i) + '\t' + str(newDict[i]) + '\n')

