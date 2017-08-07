#!/usr/bin/env python
import sys
import slurm
import argparse
import generalUtils
import os

samples = []
s = range(10)
for e in s:
    samples.append(str(e+1))

parameters = {
    "--job-name=": "Barcode",
    "-n ": 8,
    # "--mem=": 128000,
    "--mem=": 4000,
    "--time=": "1-00:00:00",
    "--output=": "./log/%A_%a.out",
    "--error=": "./log/%A_%a.err",
    "--array=": ",".join(samples),
    "--mail-type=": "END,FAIL",      # notifications for job done & fail
    "--mail-user=": "oadebali@gmail.com" # send-to address
}

job = slurm.Slurm('python barcode_simulation.py $SLURM_ARRAY_TASK_ID')
job.assignParams(parameters)
job.printScript()
job.run()

    