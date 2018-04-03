#!/usr/bin/env python

import argparse

PARSER = argparse.ArgumentParser(description='converts bed file into smaller bins')
PARSER.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
PARSER.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
PARSER.add_argument('-s', required=True, help='bin size')
PARSER.add_argument('-n', required=True, help='bin number')
PARSER.add_argument('--noShort', required= False, action='store_true', help='ignore flanking regions shorter than bin size')

ARGS = PARSER.parse_args()
BED_FILE = ARGS.i
OUT_BED_FILE = ARGS.o
BIN_SIZE = int(ARGS.s)
BIN_NUMBER = int(ARGS.n)
FILE_IN = args.i
OUT = args.o
NO_SHORT_FLAG = ARGS.noShort
INTERVAL_NUMBER = 0

for line in FILE_IN:
	INTERVAL_NUMBER += 1
	intervalName = 'interval' + str(INTERVAL_NUMBER)
	ll = line.strip().split('\t')
	chromosome = ll[0]
	newStart = start = int(ll[1])
	end = int(ll[2])
	for i in range(0, BIN_NUMBER):
		newStart = start + i * BIN_SIZE
		newEnd = newStart + BIN_SIZE
		newLine = '\t'.join([chromosome, str(newStart), str(newEnd)] + ll[3:] + [intervalName]) + '\n'
		OUT.write(newLine)
OUT.close()
