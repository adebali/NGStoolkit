#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-f', required= False, default='\t', help='delimiter')
parser.add_argument('-code', required= True, help='code')

args = parser.parse_args()

filein = args.i
out = args.o
delimiter = args.f
code = args.code

for line in filein:
    ll = line.strip().split(delimiter)
    exec(code)
