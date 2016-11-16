#!/usr/bin/env python
import os
import bed
import sys
import argparse
import generalUtils

parser = argparse.ArgumentParser(description='fix fragment length if possible')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-c1', required= True, help='<Required> start position tab for first file')
parser.add_argument('-c2', required= True, help='<Required> start position tab for second file')
parser.add_argument('-o', required=  True, help='<Required> output')

args = parser.parse_args()
input =  args.i
output = args.o
s1tab = int(args.c1) - 1
s2tab = int(args.c2) - 1

def bedClosest2distance(line, s1tab, s2tab):
	ll = line.split('\t')
	start1 = int(ll[s1tab])
	end1 = int(ll[s1tab + 1])
	start2 = int(ll[s2tab])
	end2 = int(ll[s2tab + 1])
	pos1 = generalUtils.mean([start1, end1])
	pos2 = generalUtils.mean([start2, end2])
	distance = pos2 - pos1
	return str(distance)

generalUtils.lineBasedFileOperation(input, output, bedClosest2distance, [s1tab, s2tab])