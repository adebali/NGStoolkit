#!/usr/bin/env python
import os
import sys

bedFile = sys.argv[1]
outBedFile = sys.argv[2]
multiplyFactor = int(sys.argv[3])

countTabNumber = 7
begTab = 2
endTab = 3


filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    ll = line.split('\t')
    beg = int(ll[begTab - 1])
    end = int(ll[endTab - 1])
    count = int(ll[countTabNumber - 1])

    regionLength = end - beg + 1
    newCount = (float(count) / regionLength) * float(multiplyFactor)

    newLine = ''
    for i in range(len(ll)):
        if i == begTab - 1:
            field = str(beg)
        elif i == endTab - 1:
            field = str(end)
        elif i == countTabNumber - 1:
            field = str(newCount)
        else:
            field = ll[i]

        newLine += field + '\t'

    newLine += '\n'
    out.write(newLine)

out.close()
