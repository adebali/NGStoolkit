from glob import glob
import os
import generalUtils
import csv
groupDict = generalUtils.table2dictionary('samples.csv', 'group')

def group2mergedChmm(groups):
    filesToBeMerged = []
    headers = ['chromosome', 'start', 'end', 'strand', 'gene']
    strGroups = []
    for e in groups:
        strGroups.append(str(e))
    print(strGroups)
    groupName = '_'.join(strGroups)
    headersFile = os.path.join('dataDir', groupName + '.tempTable_headers.txt') 
    headerOut = open(headersFile, 'w')
    for group in groups:
        group = str(group)
        for sampleDict in sorted(groupDict[group], key=lambda k: int(k['no'])):
            sample = sampleDict['sample'].replace(".1.fastq", "")
            treatment = sampleDict['treatment_title']
            headers.append(treatment + "_TS")
            headers.append(treatment + "_NTS")
            filesToBeMerged.append(os.path.realpath(os.path.join('dataDir/0125/', sample + '.1.cu.bo.mm10.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.geDaSi.seSt_Plus.geMa.caTSvNTS.txt')))
    tableFile = os.path.join('dataDir', groupName + '.tempTable.txt')
    code = 'mergeFileColumns.py -input ' + ' '.join(filesToBeMerged) + ' -o ' + tableFile + ' -same 1 2 3 4 5 -merge 6 7'
    print(code)
    os.system(code)

    print(headers)
    headerOut.write('\t'.join(headers))
    headerOut.close()

    os.system("cat " + headersFile + " " + tableFile + " >" + groupName + ".genes.txt")
    os.system("rm " + headersFile)
    os.system("rm " + tableFile)


group2mergedChmm([2])
