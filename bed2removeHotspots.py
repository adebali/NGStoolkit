#!/usr/bin/env python
import sys
import os
import argparse

if '--test' in sys.argv:
    testDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utils', 'testFiles')
    bedFile = os.path.join(testDir, 'bedExampleHotspots.bed')
    output = os.path.join(testDir, 'bedExampleHotspots_removedHotspots.bed')
    hotspotCutoff = 2
else:
    parser = argparse.ArgumentParser(description='remove intervals that are more than N')
    parser.add_argument('-i', required= True, help='input')
    parser.add_argument('-o', required= True, help='output')
    parser.add_argument('-n', required=  True, help='the cutoff of number of reads to be considered as a hotspot')
    args = parser.parse_args()
    bedFile = args.i
    output = args.o
    hotspotCutoff = args.n

filein = open(bedFile, 'r')

currentLineIsPreviousLine = False
currentLine = 'currentLine'
previousLine = 'currentLine'
readCount = 0
startFlag = True

out = open(output, 'w')

def printLines(line, number, out):
    for i in range(number):
        out.write(line + '\n')

for line in filein:
    line = line.strip()
    currentLine = '\t'.join(line.split('\t')[0:2])
    currentLineIsPreviousLine = currentLine == previousLine
    if currentLineIsPreviousLine or startFlag:
        readCount += 1
    else:
        if readCount < hotspotCutoff:
            printLines(previousFullLine, readCount, out)
        readCount = 1
    startFlag = False
    previousLine = currentLine
    previousFullLine = line
else:
    if readCount < hotspotCutoff:
        printLines(previousFullLine, readCount, out)