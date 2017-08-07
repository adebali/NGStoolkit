from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('samples.csv'), 'sample')
def fileList2csv(fileList, output):
    myDict = {}
    mergedDict = {}

    for file in fileList:

        sampleName = os.path.basename(file).split('.')[0] + '.fastq'
        
        treatment = sampleDict[sampleName][0]['treatment_title']
        if treatment.endswith("_A") or treatment.endswith("_B"):
            lines = open(file).read().splitlines() 
            
            totalNum = len(lines[0].split('\t'))
            initialList = [0] * totalNum

            tl = treatment.split("_")
            replicate = tl[-1]
            newName = "_".join(tl[:-1]) + "_M"

            print(sampleName + ' ' + treatment)
            l0 = lines[0].strip().split('\t')
            l1 = lines[1].strip().split('\t')
            if 'Plus' in file:
                myDict[treatment + '_Pls' + '_Pos'] = l0
                myDict[treatment + '_Pls' + '_Neg'] = l1
            
                mergedDict[newName + '_Pls' + '_Pos'] =  [int(x) + int(y) for x, y in zip(mergedDict.get(newName + '_Pls' + '_Pos', initialList), l0)]
                mergedDict[newName + '_Pls' + '_Neg'] =  [int(x) + int(y) for x, y in zip(mergedDict.get(newName + '_Pls' + '_Neg', initialList), l1)]
            elif 'Minus' in file:
                myDict[treatment + '_Min' + '_Pos'] = l0
                myDict[treatment + '_Min' + '_Neg'] = l1

                mergedDict[newName + '_Min' + '_Pos'] =  [int(x) + int(y) for x, y in zip(mergedDict.get(newName + '_Min' + '_Pos', initialList), l0)]
                mergedDict[newName + '_Min' + '_Neg'] =  [int(x) + int(y) for x, y in zip(mergedDict.get(newName + '_Min' + '_Neg', initialList), l1)]

        # treatment = sampleDict[sampleName][0]['treatment_title']
        # lines = open(file).read().splitlines() 
        # print(file)
        # print(sampleName + ' ' + treatment)
        # if 'Plus' in file:
        #     myDict[treatment + '_Pls' + '_Pos'] = lines[0].split('\t')
        #     myDict[treatment + '_Pls' + '_Neg'] = lines[1].split('\t')
        # elif 'Minus' in file:
        #     myDict[treatment + '_Min' + '_Pos'] = lines[0].split('\t')
        #     myDict[treatment + '_Min' + '_Neg'] = lines[1].split('\t')

    # with open('dataDir/mergedNucleosomeData_rep' + str(i) + '.csv', 'wb') as csv_file:
    with open(output, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in myDict.items():
            writer.writerow([key] + values)

    with open(output.replace(".csv", "_mergedRep.csv"), 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in mergedDict.items():
            writer.writerow([key] + values)

# fileList = glob('dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.inWiTFBS.inToPo.txt')

keys = ['TSS', 'TES']
cell = 'NHF1'
for key in keys:
    fileList = glob('dataDir/0131/*.cu.bo.hg19.coToBa.coToBe.unSo.seSt_*.co' + key + '.biCoToPo.txt')
    fileList2csv(fileList, 'dataDir/merged' + key + '_' + cell + '.csv')

