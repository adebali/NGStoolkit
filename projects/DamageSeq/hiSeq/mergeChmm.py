from glob import glob
import os
import generalUtils
import csv

def group2mergedChmm(group):
    filesToBeMerged = []
    headersFile = os.path.join('dataDir', group + '.chmm_headers.txt') 
    headerOut = open(headersFile, 'w')
    headers = ['chromosome', 'start', 'end', 'state']
    for sampleDict in sorted(groupDict[group], key=lambda k: int(k['no'])):
        sample = sampleDict['sample']
        treatment = sampleDict['treatment_title']
        headers.append(treatment)
        filesToBeMerged.append(os.path.realpath(os.path.join('dataDir/0106', sample.replace('.fastq', '') + '.cu.bow.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.saFrBe.soBe.coCh.noCh.bed')))
    code = 'mergeFileColumns.py -input ' + ' '.join(filesToBeMerged) + ' -o ' + os.path.join('dataDir', group + '.chmm.txt') + ' -same 1 2 3 4 -merge 5'
    print(code)
    os.system(code)
    headerOut.write('\t'.join(headers))
    headerOut.close()

groupDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples_minReadCount.csv'), 'group')
# print(groupDict.keys())
# print(sorted(groupDict[group], key=lambda k: k['no']))
group2mergedChmm('1')

# for group in groupDict.keys():
    # group2mergedChmm(group)