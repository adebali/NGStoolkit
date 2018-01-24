#!/usr/bin/env python
import sys
import os
import argparse
import bed

parser = argparse.ArgumentParser(description='print windows of given size')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-w', required= True, type=int, help='windowSize')
parser.add_argument('--noShort', required= False, action='store_true', help='ignore flanking regions shorter than bin size')
parser.add_argument('-n', required= False, type=int, default=0, help='placeholder dot number')

args = parser.parse_args()
bedFile = args.i
windowSize = args.w
placeHolderNum = args.n
if args.noShort:
    noShortFlag = args.noShort
else:
    noShortFlag = False

bed.bed(bedFile).makeWindows(windowSize, noShortFlag, placeHolderNum)