#!/usr/bin/env python
import sys
import argparse

parser = argparse.ArgumentParser(description='converts chromosome sizes (.fai) to bed format')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')

args = parser.parse_args()

for line in args.i:
    ll = line.strip().split('\t')
    bedList = [ll[0], '1', ll[1], '.', '.', '.']
    args.o.write('\t'.join(bedList) + '\n')