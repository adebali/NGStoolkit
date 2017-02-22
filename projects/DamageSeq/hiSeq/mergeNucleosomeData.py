from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples.csv'), 'sample')

def fileList2csv(fileList, output):
    myDict = {}
    for file in fileList:
        print(file)
        sampleName = os.path.basename(file).split('.')[0] + '.1.fastq'
        treatment = sampleDict[sampleName][0]['treatment_title']
        values = open(file).read().splitlines() 
        if 'Plus' in file:
            myDict[treatment + '_Pls'] = values
        elif 'Minus' in file:
            myDict[treatment + '_Min'] = values

    # with open('dataDir/mergedNucleosomeData_rep' + str(i) + '.csv', 'wb') as csv_file:
    with open(output, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in myDict.items():
            writer.writerow([key] + values)
# for i in range(2,10):
for i in range(1,2):
    # fileList = glob('dataDir/0106/*cu.bow.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.geClDaOnAlrep' + str(i) + '.adExLo.nuCoPl.txt')
    # fileList = glob('dataDir/0106/*.1.cu.bo.hg19nuc.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.inWiNu.inToPo.txt')
    fileList = glob('dataDir/0106/*.1.cu.bo.hg19nuc.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.seSt_*.inWiNuPeaks.inToPo.txt')
    # fileList = glob('dataDir/0106/*.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.seSt_*.inWiNuAdar.inToPo.txt')
    # fileList2csv(fileList, 'dataDir/mergedNucleosomeData_rep.intersect_rep3.csv')
    fileList2csv(fileList, 'dataDir/mergedNucleosomeData_1K.csv')
    # fileList2csv(fileList, 'dataDir/mergedNucleosomeData_rep.intersect_adar.csv')

