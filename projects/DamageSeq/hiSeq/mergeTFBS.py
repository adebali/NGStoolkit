from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples_minReadCount.csv'), 'sample')


fileList = glob('dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.inWiTFBS.paInRe.liOfNuToCo.txt')
myDict = {}

for file in fileList:
    print(file)
    sampleName = os.path.basename(file).split('.')[0] + '.1.fastq'
    treatment = sampleDict[sampleName][0]['treatment_title']
    lines = open(file).read().splitlines()
    values = []
    for line in lines:
        ll = line.strip().split('\t')
        values.append(ll[1])
    
    if 'Plus' in file:
        myDict[treatment + '_Pls'] = values
    elif 'Minus' in file:
        myDict[treatment + '_Min'] = values

with open('dataDir/mergedTFBS.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, values in myDict.items():
        writer.writerow([key] + values)
