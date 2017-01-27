import os
import sys
import generalUtils
from fastq import fastq

groupDictionary = generalUtils.table2dictionary('dataDir/samples.csv', 'group', ',')
# print(groupDictionary)
# sys.exit()
out = open('dataDir/samples_readNumbers.csv', 'w')

newGroupDictionary = {}


for group in sorted(groupDictionary.keys()):
    newGroupDictionary[group] = []
    for sampleDict in groupDictionary[group]:
        fastqFile = 'dataDir/' + sampleDict['sample']
        fastqObject = fastq(fastqFile)
        sampleDict['sequenceCount'] = str(fastqObject.getSequenceCount())
        # sampleDict['sequenceCount'] = str(1)
        newGroupDictionary[group].append(sampleDict)

out.write(generalUtils.dictionary2table(newGroupDictionary))
out.close()






