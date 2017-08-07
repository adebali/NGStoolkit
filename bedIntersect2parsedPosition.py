#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
from collections import Counter
from bed import bedline

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='converts bed intersect result to position list.')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('--reverse', default=False, action="store_true", help="reverse flag: it will reverse the position of the negative strand")
parser.add_argument('--name', default=False, action="store_true", help="name flag: it will categorize positions based on the length")
parser.add_argument('-rpm', type=int, default=1000000, help="RPM (Reads Per Million mapped reads): value of total mapped reads.")
parser.add_argument('-nameNorm', type=int, default=1000, help="normalize by the count of aligned region: value is per N regions.")
parser.add_argument('-l', type=int, help="length of the window")
args = parser.parse_args()

filein = args.i
out = args.o
reverseFlag = args.reverse
nameFlag = args.name
windowLength = args.l
if nameFlag and not windowLength:
    raise ValueError('-r must be specified with --name')
elif nameFlag:
    windowStart = -1 * (windowLength/2)
    windowEnd = windowLength/2 - 1
    perNameCount = args.nameNorm
    mappedReadNumber = args.rpm



def line2position(line, reverseFlag):
    def reversePosition(position):
        return -1 * (position + 1)
    
    chr1 = 0
    s1 = 1
    e1 = 2

    ll = line.strip().split('\t')
    chr1name = ll[0].strip()
    chr2name = False
    for i in range(1,len(ll)-2):
        if ll[i].strip() == chr1name:
            chr2 = i
            break
    s2 = chr2 + 1
    e2 = chr2 + 2
    chromosomeA = ll[chr1]
    startA = int(ll[s1])
    endA = int(ll[e1])
    midA = (startA + endA) / 2
    startB = int(ll[s2])
    endB = int(ll[e2])
    midB = (startB + endB) / 2
    damagePosition = midB - midA
    if reverseFlag:
        bedLine = bedline(line)
        strand = bedLine.strand()
        if strand == "-":
            damagePosition = reversePosition(damagePosition)
        else:
            if strand != "+" and strand != ".":
                raise ValueError("strand is an unexpected character:" + strand)
    return damagePosition

def fillTheGaps(dictionary):
    for key in dictionary.keys():
        for i in range(windowStart, windowEnd + 1):
            dictionary[key][i] = dictionary[key].get(i, 0)        
    return dictionary

def normalize(dictionary, nameCounts):
    normalizedD = {}
    for key in dictionary.keys():
        normalizedD[key] = {}
        for i in range(windowStart, windowEnd + 1):
            count = dictionary[key][i]
            normalizedCount = ((float(count) * float(1000000)) / float(mappedReadNumber)) * float(perNameCount)/float(nameCounts[name])
            normalizedD[key][i] = normalizedCount
    return normalizedD

if not nameFlag:
    for line in filein:
        if line.strip():
            out.write(str(line2position(line.strip(), reverseFlag)) + '\n')
elif nameFlag:
    myDict = {}
    nameCounts = {}
    for line in filein:
        if line.strip():
            bedLine = bedline(line)
            name = bedLine.name()
            position = line2position(line.strip(), reverseFlag)
            myDict[name] = myDict.get(name, {})
            myDict[name][position] = myDict[name].get(position, 0) + 1
            nameCounts[name] = nameCounts.get(name, 0) + 1
    filledGaps = fillTheGaps(myDict)
    normalizedD = normalize(filledGaps, nameCounts)
    for name in normalizedD.keys():
        for i in range(windowStart, windowEnd + 1):
            count = normalizedD[name][i]
            lineValues = [str(i), str(count), name]
            out.write("\t".join(lineValues) + "\n")
out.close()

##
#echo -e 'chr1\t10\t100\ta\tb\tc\td\tchr1\t10\t60\nchr1\t10\t100\ta\tb\tc\td\tchr1\t10\t80\n' | ./bedIntersect2parsedPosition.py 
# -20
# -10
##

##echo -e 'chr1\t10\t100\ta\tb\t-\td\tchr1\t10\t60\nchr1\t10\t100\ta\tb\t+\td\tchr1\t10\t80\n' | ./bedIntersect2parsedPosition.py --reverse
#19
#-10