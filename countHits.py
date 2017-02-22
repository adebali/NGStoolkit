#!/usr/bin/env python

import os
import sys
import subprocess
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('-files', required= True, help='files wildcard')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()
files = sorted(glob(args.files))
out = args.o

def getExtension(fileName):
    return fileName.strip().split('.')[-1]

def countHits(fileName):
    extension = getExtension(fileName)
    ignore = False
    if extension == 'bed' or extension == 'bedpe' or extension == 'sam':
        pattern = "^"
    elif extension == 'fastq':
        pattern = "^+"
    elif extension == "fa":
        pattern = "^>"
    else:
        ignore = True

    if not ignore:
        count = subprocess.check_output('grep -c "' + pattern + '" ' + fileName, shell=True)
    else:
        count = "NA"
    return count

def getInitial(fileName):
    return fileName.strip().split('.')[0]

for f in files:
    out.write(f + '\t' + countHits(f) + '\n')
