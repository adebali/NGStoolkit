
from collections import defaultdict

upstream = 100
downstream = 1000
genomeLength = 4639675


def overlaps(beg1,end1,beg2,end2):
    if (beg2 > beg1 and beg2 < end1) or (end2 > beg1 and end2 < end1):
        return True
    return False
prevBeg = 0
prevEnd = 0
lineNo = 0
flagDict = defaultdict()
allLines = []
filein = open("TSS.txt", "r")
prevNewLine = ''
for line in filein:
    ll = line.strip().split("\t")
    position= int(ll[0])
    strand = ll[1]
    if strand == "+":
        beg = position - upstream
        end = position + downstream
    elif strand == "-":
        beg = position - downstream
        end = position + upstream
    else:
        raise ValueError("Unexpected strand: " + strand)
    if beg > 0 and end < genomeLength:
        newList = ["chr", str(beg), str(end), ".", ".", strand]
        newLine = "\t".join(newList)
        if not overlaps(prevBeg, prevEnd, beg, end):
            flagDict[lineNo + 1] = True
        else:
            flagDict[lineNo + 1] = flagDict[lineNo] = False
        # print(flagDict)
        prevBeg = beg
        prevEnd = end
        prevLine = newLine
    lineNo += 1
    allLines.append("\t".join(newList))

for key in sorted(dict(flagDict).keys()):
    if flagDict[key]:
        print(allLines[key - 1])