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
        
reference = "hg19"
reference = "hg19nuc"

SAMPLE_STAT_FILE = 'samples.csv'
samples = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'sample')
for sample in samples.keys():
    files = sorted([globFile('dataDir/0106/' + sample.split('.')[0] + '.1.cu.bo.' + reference + '.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.seSt_Minus.bed'), \
    globFile('dataDir/0106/' + sample.split('.')[0] + '.1.cu.bo.' + reference + '.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.seSt_Plus.bed')])
    for file in files:
        if file:
            print(minOrPls(samples[sample][0]['treatment_title'], file) + '\t' + str(bed.bed(file).getHitNum()))
