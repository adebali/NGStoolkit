import os
import generalUtils
import bed
from glob import glob
import sys

def minOrPls(treatment, fileName):
    if 'Minus' in fileName:
        extension = '_Min'
    elif 'Plus' in fileName:
        extension = '_Pls'
    else:
        raise ValueError('no minus pr plus')
    return treatment + extension


def globFile(file):
    files = glob(file)
    if files:
        return files[0]
    else:
        return None
        
SAMPLE_STAT_FILE = 'samples.csv'
samples = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'sample')
myDict = {}
for sample in samples.keys():
    files = sorted([globFile('dataDir/0131/' + sample.split('.')[0] + '.cu.bo.hg19.coToBa.coToBe.unSo.seSt_Minus.bed'), \
    globFile('dataDir/0131/' + sample.split('.')[0] + '.cu.bo.hg19.coToBa.coToBe.unSo.seSt_Plus.bed')])
    for file in files:
        if file:            
            treatment_title = samples[sample][0]['treatment_title']
            if treatment_title.endswith("_A") or treatment_title.endswith("_B"):
                newName = minOrPls(treatment_title, file).replace("_A_", "_M_").replace("_B_", "_M_")
                myDict[newName] = myDict.get(newName,0) + bed.bed(file).getHitNum()
for key in myDict.keys():
    print(key + '\t' + str(myDict[key]))