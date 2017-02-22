import os
import sys

filein = open("/nas/longleaf/home/adebali/ogun/ENCODE/Nucleosome_Gm12878_rep1.sga.bed.1K.intersect.bed")
myDict = {}

scoreColumn = 8

scoreCutoff = [4,8]

previousStartName = ''
theList = []
i = 0

for line in filein:
    i += 1
# for i in range(1000):
    # line = filein.next().strip()
    ll = line.strip().split('\t')
    chromosome = ll[0]
    start = ll[1]
    score = int(ll[scoreColumn - 1])
    startName = chromosome + "_" + start
    if startName != previousStartName:
        if theList != []:
            sortedList = sorted(theList, key=lambda k: k['score'])
            for e in sortedList:
                # if i%10 == 0:
                if i%1 == 0:
                    print(e['line'].strip())
            # maxValueLineDict = sortedList[-1]
            # currentScore = maxValueLineDict['score']
            # if currentScore >= scoreCutoff[0] and currentScore <= scoreCutoff[1]:
            #     print(maxValueLineDict['line'].strip())
        theList = []
    theList.append({'line':line, 'score': score})
    previousStartName = startName
