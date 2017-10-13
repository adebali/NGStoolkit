#!/usr/bin/env python
import os
import sys
import generalUtils
import argparse
import bed

parser = argparse.ArgumentParser(description='append NULL values to the file')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-n', required=  True, help='append file to N times with a null line')
parser.add_argument('-null', required=  False, default= "NaN", help='null value (default: NaN)')

args = parser.parse_args()

bedFile = args.i
bedObject = bed.bed(bedFile)
hitNum = bedObject.getHitNum()

if args.null:
	nullValue = args.null
else:
	nullValue = "NA"
output = args.o

totalLine = int(args.n)
lineNumberToBeAppended = totalLine - hitNum

columnNumber = bedObject.getColumnNumber()
nullLine = (nullValue + "\t") * columnNumber + "\n"

out = open(output, "w")

for line in open(bedFile):
	out.write(line)
for i in range(lineNumberToBeAppended):
	out.write(nullLine)
out.close()
