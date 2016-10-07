import os
import sys

def run(codeList, runFlag=True):
	code = list2gappedString(codeList)
	allStringList = list2allStringList(codeList)
	if runFlag:
		print('running the following code:')
		print(code)
		failedHere = os.system(code)
		if failedHere:
			sys.exit('we cannot execute the code: ' + code)
	else:
		a = 1
		print('skipping the following code:')
		print(code)

def file(input):
	if os.path.isfile(input):
		return input
	else:
		sys.exit('no such a file: ' + input)

def runGet(code):
	import subprocess
	return subprocess.check_output(code, shell=True)

def list2allStringList(theList):
	newList = []
	for e in theList:
		newList.append(str(e))
	return newList

def list2gappedString(theList):
	string = ''
	for e in theList:
		string += str(e) + ' '
	return string[:-1]

def in2out(input, oldExtension, newExtension):
	if input.endswith(oldExtension):
		output = replaceLast(input, oldExtension, newExtension)
	else:
		output = input + oldExtension
	return output

def out2log(output):
		return output + '_log.txt'

def replaceLast(source_string, replace_what, replace_with):
	head, sep, tail = source_string.rpartition(replace_what)
	return head + replace_with + tail

def writeDict(theDict, output):
    out = open(output, 'w')
    for key in theDict.keys():
        out.write(str(key) + '\t' + str(theDict[key]) + '\n')
    return 1

def linesAreDuplicate(line1, line2, columns):
	if not line1 and not line2:
		return True
	if not line1 or not line2:
		return False
	separator = '\t'
	columnList = columns.split(',')
	ll1 = line1.strip().split(separator)
	ll2 = line2.strip().split(separator)

	for field in columnList:
		if ll1[int(field) - 1].strip() != ll2[int(field) - 1].strip():
			return False
	return True
