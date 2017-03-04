from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples.csv'), 'sample')
def fileList2csv(fileList, output):
    myDict = {}
    for file in fileList:
        # print(file)
        sampleName = os.path.basename(file).split('.')[0] + '.1.fastq'
        treatment = sampleDict[sampleName][0]['treatment_title']
        values = open(file).read().splitlines() 
        print(sampleName + ' ' + treatment)        
        if 'Plus' in file:
            myDict[treatment + '_Pls'] = values
        elif 'Minus' in file:
            myDict[treatment + '_Min'] = values

    # with open('dataDir/mergedNucleosomeData_rep' + str(i) + '.csv', 'wb') as csv_file:
    with open(output, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in myDict.items():
            writer.writerow([key] + values)


# fileList = glob('dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.inWiTFBS.inToPo.txt')
keys = ['TFBS', 'DNa', 'CTCF', 'BoTFBSUp', 'BoTFBS', 'STAT3', 'NuPeaks']

cell = 'NHF1'
cell = 'GM12878'

referenceGenome = 'hg19'
# referenceGenome = 'hg19nuc'

for key in keys:
    fileList = glob('dataDir/0106/*1.cu.bo.' + referenceGenome + '.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.seSt_*.inWi' + key + '.inToPo.txt')
    fileList2csv(fileList, 'dataDir/merged' + key + '_' + cell + '.csv')

