#!/usr/bin/env python
import sys
import slurm
import argparse
import generalUtils

parser = argparse.ArgumentParser()
parser.add_argument('-g', required= False, help='group no')
parser.add_argument('-s', required= False, help='comma separated sample list')
args = parser.parse_args()

SAMPLE_STAT_FILE = 'dataDir/samples.csv'

def group2sample(group):
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

parameters = {
    "--job-name": "Damage_Seq_Pipeline",
    # "-n": 1,
    "-n": 4,
    "--mem": 32000,
    "--time": "24:00:00",
    "--output": "./log/%A_%a.out",
    "--error": "./log/%A_%a.err",
    "--array": ",".join(samples)
}

job = slurm.Slurm('python pipeline.py -n $SLURM_ARRAY_TASK_ID')
job.assignParams(parameters)
job.printScript()
job.run()