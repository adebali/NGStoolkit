from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('samples.csv'), 'sample')
def fileList2csv(fileList, output):
    myDict = {}
    for file in fileList:
        # print(file)
        sampleName = os.path.basename(file).split('.')[0]
        treatment = sampleDict[sampleName][0]['treatment_title']
        lines = open(file).read().splitlines() 
        
        print(sampleName + ' ' + treatment)        
        if 'Plus' in file:
            myDict[treatment + '_Pls' + '_Pos'] = lines[0].split('\t')
            myDict[treatment + '_Pls' + '_Neg'] = lines[1].split('\t')
        elif 'Minus' in file:
            myDict[treatment + '_Min' + '_Pos'] = lines[0].split('\t')
            myDict[treatment + '_Min' + '_Neg'] = lines[1].split('\t')

    # with open('dataDir/mergedNucleosomeData_rep' + str(i) + '.csv', 'wb') as csv_file:
    with open(output, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in myDict.items():
            writer.writerow([key] + values)


# fileList = glob('dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.inWiTFBS.inToPo.txt')

keys = ['TSS', 'TES']
cell = 'GM12878'
for key in keys:
    fileList = glob('dataDir/0308/*unSo.seSt_*.co' + key + '.biCoToPo.txt')
    fileList2csv(fileList, 'dataDir/merged' + key + '_' + cell + '.csv')

