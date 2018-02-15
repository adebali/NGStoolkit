#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedintersect
from operator import itemgetter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed intersect positions to count per position list per gene per strand.')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()


filein = args.i
out = args.o

d = {}

for line in filein:
    ll = line.strip().split('\t')
    if float(ll[4]) > 0:
        d[ll[3]] = d.get(ll[3], {
            'chromosome': ll[0],
            'start': ll[1],
            'end': ll[2],
            'name': ll[3],
            'count': {
                'opposite': [],
                'same': []
            },
            'score': float(ll[4]),
            'strand': ll[5]
        }
        )
        category = ll[7]
        d[ll[3]]['count'][category].append(int(float(ll[6])))

geneNo = 0
for tupleElement in sorted(d.iteritems(), key=lambda (x, y): y['score']):
    geneNo += 1
    # print(tupleElement)
    geneObject = tupleElement[1]
    for category in ['same', 'opposite']:
        for i in range(0,100):
            i_count = geneObject['count'][category].count(i)
            out.write(
                '\t'.join([geneObject['chromosome'], geneObject['start'], geneObject['end'], geneObject['name'], str(geneObject['score']), geneObject['strand'], category, str(geneNo), str(i), str(i_count)]) + '\n'
            )