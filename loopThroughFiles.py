<<<<<<< HEAD
#!/usr/bin/env python
=======
#!/usr/bin/python
>>>>>>> 587350ecbb57997bed8cbdecda2c4aae48b1e76c

import os
import sys
from glob import glob
import argparse
import datetime

<<<<<<< HEAD
totalCodeString = ''
for e in sys.argv:
	totalCodeString += e

if '#OUT' in totalCodeString:
	extensionRequirement = True
else:
	extensionRequirement = False

parser = argparse.ArgumentParser(description='bsub execution looping through files.')
parser.add_argument('-code', required= True, help='code that you are interested in: place #IN as input and #OUT as output')
parser.add_argument('-files', required= True, help='Files that you are interested in. Example: *.fastq')
parser.add_argument('-ignore', required= False, help='Files to be ignored. Example: Specific_file*.fastq')
parser.add_argument('-ignoreList', required= False, help='Files (separated by comma) to be ignored. Example: file1.fastq,file2.fastq')
parser.add_argument('-extension', required= extensionRequirement, help='extension that you want to add to the outputfile')
=======
parser = argparse.ArgumentParser(description='bsub execution looping through files.')
parser.add_argument('-code', required= True, help='code that you are interested in: place #IN as input and #OUT as output')
parser.add_argument('-files', required= True, help='Files that you are interested in. Example: *.fastq')
parser.add_argument('-extension', required= True, help='extension that you want to add to the outputfile')
>>>>>>> 587350ecbb57997bed8cbdecda2c4aae48b1e76c
parser.add_argument('-log', help='optional log file to save the commands')
parser.add_argument('--mock', action='store_true', help = 'Do not run anything, only print commands')


args = parser.parse_args()
<<<<<<< HEAD
fileListInput = set(glob(args.files))

toBeIgnoredFiles = set()

if args.ignore:
	toBeIgnoredFiles = toBeIgnoredFiles.union(set(glob(args.ignore)))
if args.ignoreList:
	toBeIgnoredFiles = toBeIgnoredFiles.union(set(args.ignoreList.split(',')))

fileList = list(fileListInput - toBeIgnoredFiles)
print(fileListInput)
=======
fileList = glob(args.files)

>>>>>>> 587350ecbb57997bed8cbdecda2c4aae48b1e76c
outputDirectory = '/'.join(fileList[0].split('/')[:-1])

if args.log:
	out = open(outputDirectory + '/' + args.log, 'w')
	out.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + '\n')

for fileName in fileList:
	code = args.code
<<<<<<< HEAD

	if extensionRequirement:
		outputFile = fileName + '.' + args.extension
	else:
		outputFile = fileName + '.out'

	logFile = outputFile + '.log'
	code = code.replace("#IN", fileName).replace("#OUT", outputFile)
	if code.split(' ')[0] == 'bsub':
		code = code.replace('bsub ', 'bsub -o ' + logFile + ' ')
	# print(code)
=======
	code = "bsub " + code.replace("#IN", fileName).replace("#OUT", fileName + "." + args.extension)
	print(code)
>>>>>>> 587350ecbb57997bed8cbdecda2c4aae48b1e76c
	if args.log:
		out.write(code + '\n')
	if not args.mock:
		os.system(code)
<<<<<<< HEAD
		# print('The job is submitted')
	# else:
	# 	print('The job is mocked')
=======
		print('The job is submitted')
	else:
		print('The job is mocked')
>>>>>>> 587350ecbb57997bed8cbdecda2c4aae48b1e76c


if args.log:
	out.close()
