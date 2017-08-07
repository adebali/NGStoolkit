from glob import glob
import os
import generalUtils
import csv
groupDict = generalUtils.table2dictionary('samples.csv', 'group')

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
            filesToBeMerged.append(os.path.realpath(os.path.join('dataDir/0308', sample + '.unSo.coCh.bed')))
    code = 'mergeFileColumns.py -input ' + ' '.join(filesToBeMerged) + ' -o ' + os.path.join('dataDir', groupName + '.chmm.txt') + ' -same 1 2 3 4 5 -merge 6'
    print(code)
    os.system(code)

    print(headers)
    headerOut.write('\t'.join(headers))
    headerOut.close()


group2mergedChmm([1])
group2mergedChmm([2])
