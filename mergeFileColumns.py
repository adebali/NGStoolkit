#!/usr/bin/env python
# DOESN'T WORK!!!!!!!!!!!!
import os
import sys
import fastq
import generalUtils
import argparse
from subprocess import Popen, PIPE, STDOUT

parser = argparse.ArgumentParser(description='writes sequence length distribution')
parser.add_argument('-inputs', nargs='+', help='<Required> Input file list (min 2)', required=True)
parser.add_argument('-o', required= True, help='<Required> Output')
parser.add_argument('-same', nargs='+', required= True, help='<Required> identical (same) column numbers (eg. 1 3 for identical first and third columns)')
parser.add_argument('-merge', nargs='+', required= True, help='<Required> to be merged column no list')

args = parser.parse_args()
inputs = args.inputs
output = args.o
identicalColumns = args.same
mergingColumns = args.merge

code = ['paste <(cat ' + os.path.realpath(inputs[0]) + ' | cut -f ' + ','.join(identicalColumns) + ') ']
for input in inputs:
    code .append('<(cat ' + os.path.realpath(input) + ' | cut -f ' + ','.join(mergingColumns) + ') ')
code.append(' > ' + os.path.realpath(output))
out = open('temp.sh','w')
out.write(' '.join(code))
out.close()

os.system('bash temp.sh && rm temp.sh')