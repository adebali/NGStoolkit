#!/usr/bin/env python
import os
import sys
import bed
import argparse
import random
import tempfile

parser = argparse.ArgumentParser(description='add NTS reads to TS')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-transcribed', required= True, help='transcribed region list bed file')

args = parser.parse_args()

bedFile = args.i
output = args.o
transcribedFile = args.transcribed
intermediateTemporaryFile = tempfile.NamedTemporaryFile()
intermediateFile = intermediateTemporaryFile.name

codeList = [
    'bedtools',
    'intersect',
    '-a', transcribedFile,
    '-b', bedFile,
    '-s',
    '-wb',
    '-F', '0.49',
    '>', intermediateFile
]

os.system(' '.join(codeList))

temporaryOutputFile = tempfile.NamedTemporaryFile()
unsortedOutputFile = temporaryOutputFile.name

os.system('cat ' + bedFile + ' ' + intermediateFile + ' > ' + unsortedOutputFile)
os.system('sort -k1,1 -k2,2n -k3,3n ' + unsortedOutputFile + ' > ' + output)

intermediateTemporaryFile.close()
temporaryOutputFile.close()