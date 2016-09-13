#!/usr/bin/python

import os
import sys
from glob import glob
import argparse
import datetime

parser = argparse.ArgumentParser(description='bsub execution looping through files.')
parser.add_argument('-code', required= True, help='code that you are interested in: place #IN as input and #OUT as output')
parser.add_argument('-files', required= True, help='Files that you are interested in. Example: *.fastq')
parser.add_argument('-extension', required= True, help='extension that you want to add to the outputfile')
parser.add_argument('-log', help='optional log file to save the commands')
parser.add_argument('--mock', action='store_true', help = 'Do not run anything, only print commands')


args = parser.parse_args()
fileList = glob(args.files)

outputDirectory = '/'.join(fileList[0].split('/')[:-1])

if args.log:
	out = open(outputDirectory + '/' + args.log, 'w')
	out.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + '\n')

for fileName in fileList:
	code = args.code
	code = "bsub " + code.replace("#IN", fileName).replace("#OUT", fileName + "." + args.extension)
	print(code)
	if args.log:
		out.write(code + '\n')
	if not args.mock:
		os.system(code)
		print('The job is submitted')
	else:
		print('The job is mocked')


if args.log:
	out.close()
