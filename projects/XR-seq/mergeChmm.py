from glob import glob
import os
import generalUtils
import csv

def group2mergedChmm(groups):
    filesToBeMerged = []
    headers = ['chromosome', 'start', 'end', 'state', 'color']
    strGroups = []
    for e in groups:
        strGroups.append(str(e))
    print(strGroups)
    groupName = '_'.join(strGroups)
    headersFile = os.path.join('dataDir', groupName + '.chmm_headers.txt') 
    headerOut = open(headersFile, 'w')
    for group in groups:
        group = str(group)
        for sampleDict in sorted(groupDict[group], key=lambda k: int(k['no'])):
            sample = sampleDict['sample']
            treatment = sampleDict['treatment_title']
            headers.append(treatment)
            filesToBeMerged.append(os.path.realpath(os.path.join('dataDir/0131', sample.replace('.fastq', '') + '.cu.bo.hg19.coToBa.coToBe.unSo.coCh.bed')))
    code = 'mergeFileColumns.py -input ' + ' '.join(filesToBeMerged) + ' -o ' + os.path.join('dataDir', groupName + '.chmm.txt') + ' -same 1 2 3 4 5 -merge 6'
    # print(code)
    os.system(code)

    print(headers)
    headerOut.write('\t'.join(headers))
    headerOut.close()

groupDict = generalUtils.table2dictionary(generalUtils.file('dataDir/samples.csv'), 'group')
# print(groupDict.keys())
# print(sorted(groupDict[group], key=lambda k: k['no']))
group2mergedChmm([1])
group2mergedChmm([2])

# for group in groupDict.keys():
    # group2mergedChmm(group)