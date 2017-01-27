import os
import sys
import generalUtils
from fastq import fastq

groupDictionary = generalUtils.table2dictionary('dataDir/samples_readNumbers.csv', 'group', ',')
# print(groupDictionary)
# sys.exit()
out = open('dataDir/samples_minReadCount.csv', 'w')

newGroupDictionary = {}

maxValue = 1000000000000
groupMinDict = {}

for group in sorted(groupDictionary.keys()):
    groupMinimumReadNum = maxValue
    groupMinDict[group] = groupMinimumReadNum
    for sampleDict in groupDictionary[group]:
        readCount = int(sampleDict['sequenceCount'])
        groupMinDict[group] = min(groupMinDict[group], readCount)

for group in sorted(groupDictionary.keys()):
    newGroupDictionary[group] = []
    for sampleDict in groupDictionary[group]:
        sampleDict['minReadCount'] = groupMinDict[group]
        newGroupDictionary[group].append(sampleDict)

out.write(generalUtils.dictionary2table(newGroupDictionary))
out.close()






