import os
import sys
import generalUtils
import slurm

parameters = {
    "--job-name=": "mergeSamples",
    "-n ": 1,
    "--mem=": 8000,
    "--time=": "2-00:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}

extension = '.cu.bo_TAIR10.coToBa.coToBe.unSo.bed'

sampleDictionary = generalUtils.table2dictionary('samples.csv', 'group')
# i = 0
# for group in ['2','3']:
#     i += 1
#     l = sampleDictionary[group]
#     samples = []
#     for e in l:
#         samples.append('dataDir/0707/' + e['sample'].replace('.fastq', '') + extension)
#     code = 'cat ' + ' '.join(samples) + ' | sort -k1,1 -k2,2n -k3,3n >' + 'dataDir/0707/experiment' + str(i) + extension 
#     job = slurm.Slurm(code)
#     job.assignParams(parameters)
#     job.printScript()
#     job.run()

for group in ['4']:
    i = 12
    l = sampleDictionary[group]
    samples = []
    for e in l:
        samples.append('dataDir/0707/' + e['sample'].replace('.fastq', '') + extension)
    code = 'cat ' + ' '.join(samples) + ' | sort -k1,1 -k2,2n -k3,3n >' + 'dataDir/0707/experiment' + str(i) + extension 
    job = slurm.Slurm(code)
    job.assignParams(parameters)
    job.printScript()
    job.run()
