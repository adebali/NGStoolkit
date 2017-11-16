#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedintersect

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed intersect result to distance percentage or abssolute (binned) position.')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('--flanking', action='store_true', default=False, help='''
    if bed file is for the flanking regions of interest such as TSS upstream and TES downstream.
    Input must go line by line: upstream, downstream, up, down ...
    ''')
parser.add_argument('--fixed', action='store_true', default=False, help='should be used when for flankings we want to get absolute range')
# parser.add_argument('-length', type=int, default=1000, help='if bed file is for the flanking regions of interest such as TSS upstream and TES downstream')
parser.add_argument('-w', type=int, default=1, help='bin (window) size')

args = parser.parse_args()


filein = args.i
out = args.o

for line in filein:
    intersectLine = bedintersect(line)
    name = intersectLine.name()

    if args.flanking:
        if '_upstream' in name:
            region = 'upstream'
        elif '_downstream' in name:
            region = 'downstream'
        else:
            raise ValueError("Expected keyword (up/downstream) not found in name: " + name)
    else:
        region = 'main'

    if region == 'main':
        distancePercentage = intersectLine.getDistancePercentage()
        distance = distancePercentage
    if region != 'main':
        if args.fixed:
            distance = intersectLine.getAbsoluteDistance(region, args.w)
        else:
            distancePercentage = intersectLine.getDistancePercentage()   
            if region == 'upstream':
                distancePercentage -= 100
            elif region == 'downstream':
                distancePercentage += 100
            else:
                raise ValueError('Unexpected region ' + region)
            distance = distancePercentage
    sameStrands = intersectLine.sameStrands()
    if sameStrands:
        strands = 'same'
    else:
        strands = 'opposite'
    allFields = intersectLine.fields()
    name1 = allFields[3]
    fields = allFields[6:]
    fields[3] = name1.replace('_upstream', '').replace('_downstream', '')
    fields.append(str(distance))
    fields.append(strands)
    newLine = intersectLine.fields2line(fields)
    out.write(newLine + '\n')
    