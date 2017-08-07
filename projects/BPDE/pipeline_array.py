#!/usr/bin/env python
import sys
import slurm
import argparse
import generalUtils
import os

parser = argparse.ArgumentParser()
parser.add_argument('-g', required= False, help='group no')
parser.add_argument('-s', required= False, help='comma separated sample list')
parser.add_argument('-e', required= False, help='comma separated samples to be excluded')
parser.add_argument('--report', action='store_true',help='get report only')
args = parser.parse_args()

reportFlag = args.report

SAMPLE_STAT_FILE = './samples.csv'

def group2sample(group):
    if group == '0':
        return generalUtils.sorted_nicely(generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'no').keys())
    groupDict = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'group')
    samples = []
    for sample in groupDict[group]:
        samples.append(sample["no"])
    return samples

if args.g:
    group = args.g
    samples = group2sample(group)
elif args.s:
    samples = args.s.split(',')
else:
    raise ValueError("Either -g or -s must be given.")

if args.e:
    exclude = args.e.split(',')
    samples = list(set(samples) - set(exclude))

parameters = {
    "--job-name=": "Damage_Seq_Pipeline",
    "-n ": 8,
    # "--mem=": 128000,
    "--mem=": 32000,
    "--time=": "2-00:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--array=": ",".join(samples),
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}
print(samples)
if not reportFlag:
    job = slurm.Slurm('python pipeline.py -n $SLURM_ARRAY_TASK_ID')
    job.assignParams(parameters)
    job.printScript()
    job.run()
else:
    os.system('rm report.txt')
    os.system('python pipeline.py -n 1 --outputCheck | cut -d\' \' -f 1 | xargs | sed -e \'s/ /,/g\' >report.txt')
    for no in samples:
        code = 'python pipeline.py -n ' + str(no) + ' --outputCheck | cut -d\' \' -f 3 | xargs | sed -e \'s/ /,/g\' >>report.txt'
        print(code)
        os.system(code)
        