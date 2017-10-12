#!/usr/bin/env python
import os
import sys
import argparse
import bed
import math
import numpy 

parser = argparse.ArgumentParser(description='get count ratio: log2(n/(median of all intervals))')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

countTabNumber = 7
bedObject = bed.bed(args.i)
out = open(args.o, 'w')

valueList = []
for bedline in bedObject.read():
	valueList.append(float(bedline.fields()[countTabNumber - 1]))
median = numpy.median(valueList)


bedObject = bed.bed(args.i)
for bedline in bedObject.read():
	fields = bedline.fields()
	count = float(fields[countTabNumber - 1])
	ratio = float(count) / median
	
	try:
		logValue = math.log(ratio)
	except: 
		logValue = "NA"

	fields[countTabNumber - 1] = str(logValue)
	out.write(bedline.fields2line(fields))
	out.write('\n')