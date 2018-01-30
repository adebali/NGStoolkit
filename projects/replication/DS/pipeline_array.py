#!/usr/bin/env python
import sys
import slurm
import argparse
import generalUtils
import os

parser = argparse.ArgumentParser()
parser.add_argument('-g', nargs='+', type=int, required= False, help='group no')
parser.add_argument('-s', required= False, type=str, nargs='+', help='space separated sample list')
parser.add_argument('-e', required= False, type=str, nargs='+', help='space separated samples to be excluded')
parser.add_argument('-samples', required= False, default='samples.csv', help='sample list (default:samples.csv)')
parser.add_argument('-script', required= False, default='pipeline.py', help='script command (default:pipeline.py)')
parser.add_argument('--report', action='store_true',help='get report only')
parser.add_argument('-reportDir', required= False, default='dataDir', help='report directory')

args = parser.parse_args()
reportFlag = args.report

SAMPLE_STAT_FILE = args.samples
if not os.path.isfile(SAMPLE_STAT_FILE):
    raise ValueError('No such a file: ' + SAMPLE_STAT_FILE)

def group2sample(group):
    if group == '0':
        return generalUtils.sorted_nicely(generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'no').keys())
    groupDict = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'group')
    samples = []
    for sample in groupDict[group]:
        samples.append(sample["no"])
    return samples

if args.g:
    groups = args.g
    samples = []
    for group in groups:
        samples += group2sample(str(group))
elif args.s:
    samples = args.s
else:
    raise ValueError("Either -g or -s must be given.")

if args.e:
    exclude = args.e
    samples = list(set(samples) - set(exclude))

os.system('mkdir -p ./log')


pipelineParameters = {
    "--job-name=": "DS_replication",
    "-n ": 1,
    # "-n ": 8,
    "--mem=": 16000,
    "--time=": "1-00:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--array=": ",".join(samples),
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}

catParameters = {
    "--job-name=": "DS_rep_cat",
    "-n ": 1,
    "--mem=": 4000,
    "--time=": "02:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--array=": ",".join(samples),
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}

print(samples)
if not reportFlag:
    job = slurm.Slurm('python ' + args.script + ' run -n $SLURM_ARRAY_TASK_ID')
    job.assignParams(pipelineParameters)
    job.printScript()
    jobId = job.run()

    # catJob = slurm.Slurm('python ' + args.script + ' cat')
    # catJob.assignParams(catParameters)
    # catJob.setDependencies([jobId])
    # catJob.printScript()
    # catJobId = catJob.run()

else:
    reportDirectory = args.reportDir
    os.system('rm ' + os.spathjoin(reportDirectory, 'report.txt'))
    os.system('python pipeline.py run -n 1 --outputCheck | cut -d\' \' -f 1 | xargs | sed -e \'s/ /,/g\' >' + os.spathjoin(reportDirectory, 'report.txt'))
    for no in samples:
        code = 'python pipeline.py run -n ' + str(no) + ' --outputCheck | cut -d\' \' -f 3 | xargs | sed -e \'s/ /,/g\' >>' + os.path.join(reportDirectory, 'report.txt')
        print(code)
        os.system(code)