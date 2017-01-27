

import os
import sys
from glob import glob
import generalUtils

samples = generalUtils.table2dictionary(generalUtils.file('dataDir/samples.csv'), 'no')
headers = ['chromosome', 'start', 'end', 'state']
columnNo = 5
headers = headers[0:columnNo - 1]
files = []
for no in sorted(samples.keys()):
    sample = samples[no][0]
    treatment = sample['treatment_title']
    headers.append(treatment)
    files.append(glob('dataDir/0110/' + sample['sample'].split('.')[0] + '*coWiHi.bed')[0])




columns = []
for i in range(columnNo):
    columns.append(str(i + 1))

cutString = ','.join(columns)
for i in range(len(files) - 1):
    cutString += ',' + str(columnNo + columnNo*(i+1))

headerText = '\t'.join(headers)
headersOut = open('dataDir/0110/mergedChmm_headers.csv', 'w')
headersOut.write(headerText)
headersOut.close()

commandList = ['paste'] + files + ['|'] + ['cut -f ' + cutString] + ['>'] + ['dataDir/0110/mergedChmm.csv']
print(' '.join(commandList))
os.system(' '.join(commandList))

