#!/usr/bin/env python
import sys
import os
import argparse
import bed
import generalUtils

parser = argparse.ArgumentParser(description='takes bed as input, get the middle point and extend it to both sides')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-w', required= True, help='windowSize')
parser.add_argument('--randomMid', required= False, action='store_true', help='for cases of .5 middle point, randomly select between positions 0 or 1')
args = parser.parse_args()
bedFile = args.i
output = args.o
windowSize = int(args.w)

def getInterval(line, randomness=False):
    bedLine = bed.bedline(line)
    middlePoint = bedLine.midpoint()
    start = middlePoint - windowSize
    end = middlePoint + windowSize
    return bedLine.newline(start, end)

generalUtils.lineBasedFileOperation(bedFile, output, getInterval, [])
