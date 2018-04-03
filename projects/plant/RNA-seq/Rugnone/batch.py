import os
import json
import slurm
import argparse
from sampleHelpers import samples

parser = argparse.ArgumentParser(description='Batch SLURM submission of the general master pipeline', prog="batch.py")
parser.add_argument('-json', default='sample.json', required= False, help='file for sample configuration in JSON format')
parser.add_argument('-e', required= True, type=int, help='Experiment number')
parser.add_argument('-s', required= False, type=int, nargs='+', help='space separated sample list')
parser.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')

args = parser.parse_args()
sampleClass = samples(args.json)
if args.s:
    sampleNoList = args.s
else:
    sampleNoList = sampleClass.experiment2sampleNoList(args.e)



pipelineParameters = {
    "--job-name=": "db_pipeline",
    "-n ": 8,
    # "-n ": 1,
    "--mem=": 16000,
    "--time=": "1-00:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--array=": ','.join(str(x) for x in sampleNoList),
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}

job = slurm.Slurm('python pipeline.py run -e ' + str(args.e) + ' -n $SLURM_ARRAY_TASK_ID')
job.assignParams(pipelineParameters)
job.printScript()
if not args.mock:
    jobId = job.run()
