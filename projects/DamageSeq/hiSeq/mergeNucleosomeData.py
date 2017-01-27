from glob import glob
import os
import generalUtils
import csv

sampleDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples_minReadCount.csv'), 'sample')

for i in range(2,10):

    fileList = glob('dataDir/0106/*cu.bow.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.seSt_*.geClDaOnAlrep' + str(i) + '.adExLo.nuCoPl.txt')
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

    with open('dataDir/mergedNucleosomeData_rep' + str(i) + '.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, values in myDict.items():
            writer.writerow([key] + values)