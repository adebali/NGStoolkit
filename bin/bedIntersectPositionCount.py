#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedintersect

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed intersect result to distance percentage.')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-count', type=int, required=True, help='count column')
parser.add_argument('-cat', type=int, nargs="+", required=True, help='category column')
parser.add_argument('-scale', nargs="+", type=float, default=[1.0], help='scale factors for multiplication of the count')
parser.add_argument('-ends', type=str, required=False, default='keep', action='store', choices=['keep', 'double', 'remove'], help='what to do with the ends: double remove keep')
args = parser.parse_args()


filein = args.i
out = args.o
categoryColumns = args.cat



d = {}



scaleFactor = 1
for factor in args.scale:
    scaleFactor = scaleFactor * factor

for line in filein:
    if line.strip() != '':
        ll = line.strip().split('\t')
        position = int(float(ll[args.count - 1].strip()))
        categories = []
        for catCol in categoryColumns:
            categories.append(ll[catCol - 1].strip().replace('|', '_'))
        category_string = '|'.join(categories)
        d[category_string] = d.get(category_string, {})
        d[category_string][position] = d[category_string].get(position, 0) + 1

for category in sorted(d.keys()):
    for position in sorted(d[category].keys()):
        if position == 0 or position == 100:
            if args.ends == 'double':
                count = d[category][position] * 2
                writeFlag = True
            elif args.ends == 'remove':
                count = 0
                writeFlag = False
            elif args.ends == 'keep':
                count = d[category][position]
                writeFlag = True
            else:
                raise ValueError('Unexpected -ends argument. Must be "double", "remove" or "keep"(default), but had ' + args.ends)
        else: 
            count = d[category][position]            
            writeFlag = True

        if writeFlag:
            categories = category.split('|')
            out.write(str(position) + '\t' + str(count * float(scaleFactor)) + '\t' + '\t'.join(categories) + '\n')
