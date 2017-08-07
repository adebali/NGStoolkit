#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='add column(s) of same value to the table')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-c', nargs="+", required= True, help='column value(s)')
args = parser.parse_args()

filein = args.i
out = args.o
columns = args.c
separator = "\t"

for line in filein:
    ll = line.strip().split(separator)
    newList = ll + columns
    out.write(separator.join(newList) + "\n")