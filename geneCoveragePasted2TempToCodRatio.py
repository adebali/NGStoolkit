#!/usr/bin/env python
import os
import sys

bedFile = sys.argv[1]
filein = open(bedFile, 'r')
for line in filein:
    ll = line.strip().split('\t')
    strand = ll[3]
    positiveReadNum = int(ll[4])
    negativeReadNum = int(ll[9])
    if strand == '+':
        templateReadNum = negativeReadNum
        codingReadNum = positiveReadNum
    elif strand == '-':
        templateReadNum = positiveReadNum
        codingReadNum = negativeReadNum
    templateToCodingRatio = float(templateReadNum)/max(1,codingReadNum)
    newll = ll[0:4] 
    newll.append(str(templateToCodingRatio))
    newLine = '\t'.join(newll)
    print(newLine)
