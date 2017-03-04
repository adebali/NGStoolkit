#!/usr/bin/env python

import os
import sys
import generalUtils
from glob import glob

sampleDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples.csv'), 'group')

def getStrand(fileName):
    if 'Plus' in fileName:
        return 'Plus'
    if 'Minus' in fileName:
        return 'Minus'
    else:
        return None

extension = ".toBa.bam"
# extension = ".toBa.bam.bai"
groups = ['1','2']
for group in groups:
    groupList = sampleDict[group]
    for sample in groupList:
        cell = sample["cell"]
        product = sample["product"]
        replicate = sample["rep"]
        treatment = sample["treatment"]
        wildcard = "dataDir/0106/" + sample["sample"].replace(".fastq", "") + "*" + extension
        print(wildcard)
        filename = glob(wildcard)[0]
        strand = getStrand(filename)
        newName = "_".join([cell, product, treatment, replicate, strand]) + extension
        print(filename)
        code = "cp " + filename + " dataDir/0106/" + newName
        print(code)
        if '--run' in sys.argv:
            os.system(code)