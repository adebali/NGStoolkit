#!/usr/bin/env python
import os
import sys

bedFile = sys.argv[1]
filein = open(bedFile, 'r')
for line in filein:
    ll = line.strip().split('\t')
    strand = ll[3]
    positiveReadNum = float(ll[4])
    negativeReadNum = float(ll[9])
    if strand == '+':
        templateReadNum = negativeReadNum
        codingReadNum = positiveReadNum
    elif strand == '-':
        templateReadNum = positiveReadNum
        codingReadNum = negativeReadNum
    templateMinusCoding = templateReadNum - codingReadNum
    newll = ll[0:4] 
    newll.append(str(templateMinusCoding))
    newLine = '\t'.join(newll)
    print(newLine)
