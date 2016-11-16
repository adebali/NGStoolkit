import os
import sys
import re

def sorted_nicely(l):
	""" Sorts the given iterable in the way that is expected.

	Required arguments:
	l -- The iterable to be sorted.

	"""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key = alphanum_key)

def run(codeList, runFlag=True):
	""" System call of the given list
	
	Arguments:
	codeList -- list of the code
	Optional:
	runFlag -- Set False for mocking.
	
	"""
	code = list2gappedString(codeList)
	allStringList = list2allStringList(codeList)
	if runFlag:
		# print('running the following code:')
		print('-->\t' + code)
		failedHere = os.system(code)
		if failedHere:
			sys.exit('we cannot execute the code: ' + code)
	else:
		a = 1
		# print('skipping the following code:')
		print('--|\t' + code)

def file(input):
	if os.path.isfile(input):
		return input
	else:
		sys.exit('no such a file: ' + input)

def runGet(code):
	import subprocess
	if not '--mock' in sys.argv:
		return subprocess.check_output(code, shell=True)
	else:
		print(code)
		return False

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

def in2tempOut(input, oldExtension, newExtension, tempDir):
	if input.endswith(oldExtension):
		output = replaceLast(input, oldExtension, newExtension)
	else:
		output = input + oldExtension
	temporaryOutput = os.path.join(tempDir, os.path.basename(output))
	return temporaryOutput

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

def line2switchColumnsBasedOnValue(line, valueColNo, value, col1no, col2no, separator = '\t'):
	ll = line.strip().split(separator)
	newll = list(ll)
	if ll[valueColNo] == value:
		 newll[col1no] = ll[col2no]
		 newll[col2no] = ll[col1no]
	return separator.join(newll)

def lineBasedFileOperation(input, output, function, arguments):
	filein = open(input, 'r')
	out = open(output, 'w')
	for line in filein:
		currentArguments = list([line] + arguments)
		newLine = function(*currentArguments)
		if newLine:
			out.write(newLine.strip() + '\n')
	out.close()

def dna2reverseComplement(seq):
	seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
				'a':'t', 't':'a', 'g':'c', 'c':'g',
				'n':'n', 'N':'N', 'x':'x', 'X':'X'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def subseq(seq, start, end):
	return seq[start - 1 : end]

def motifCount(sequence, motif):
	count = start = 0
	while True:
		start = sequence.find(motif, start) + 1
		if start > 0:
			count+=1
		else:
			return count

def mean(numbers):
	return float(sum(numbers)) / max(len(numbers), 1)

def list2concatString(l, string):
	newList = []
	for e in l:
		newList.append(e + string)
	return newList