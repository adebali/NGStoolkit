from glob import glob
import os
import generalUtils
import csv
import subprocess
numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,17,18,19,20,21,22,23,24,25,26,27,28,31,34,37,38,41,44,47,48,57,58,59,60]

sampleDict = generalUtils.table2dictionary(generalUtils.file('samples.csv'), 'no')

# for i in range(1,61):
#     fileName = 'dataDir/0106/' + sampleDict[str(i)][0]['sample'].replace('.fastq', '') + '.cu.bo.hg19.coToBa.coToBe.coBeToSiFr.geInSi.siLiIn.txt'
#     if os.path.isfile(fileName):
#         result = str(i) + ' ' + open(fileName).readline().strip()
#     else:
#         result = str(i) + ' NA NA'
#     print(result)

# for i in range(1,61):
#     file1 = 'dataDir/raw/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '1.fastq') 
#     file2 = 'dataDir/raw/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '2.fastq') 
#     if os.path.isfile(file1):
#         md5_1 = subprocess.check_output('md5sum ' + file1, shell=True)
#         md5_2 = subprocess.check_output('md5sum ' + file2, shell=True)
#     else:
#         md5_1 = md5_2 = "NA"
#     print(str(i) + '\t' + md5_1 + '\t' + md5_2)

# for i in numbers:
#     f = 'dataDir/0106/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '')  + '1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed'
#     md5 = subprocess.check_output('md5sum ' + f, shell=True)
#     print(str(i) + '\t' + f + '\t' + md5)

os.system("mkdir -p /pine/scr/a/d/adebali/GEO/oadebali@gmail.com")

for i in numbers:
    code = 'cp ' + 'dataDir/raw/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '1.fastq') + ' mkdir -p /pine/scr/a/d/adebali/GEO/oadebali@gmail.com/'
    os.system(code)
    code = 'cp ' + 'dataDir/raw/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '2.fastq') + ' mkdir -p /pine/scr/a/d/adebali/GEO/oadebali@gmail.com/'
    os.system(code)
    code = 'cp ' + 'dataDir/0106/' + sampleDict[str(i)][0]['sample'].replace('1.fastq', '')  + '1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed' + ' mkdir -p /pine/scr/a/d/adebali/GEO/oadebali@gmail.com/'
    os.system(code)
