#!/usr/bin/env python
import os
import sys

bedFile = sys.argv[1]
outBedFile = sys.argv[2]

countTabNumber = 7

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
	ll = line.split('\t')
	count = float(ll[countTabNumber - 1])
	if count != 0:
		print(line.strip())
out.close()
