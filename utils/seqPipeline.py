# Requirements that should be in path:
# samtools
# bedtools

import os
import sys
import fastq
import bed
import fasta
import generalUtils


class SeqPipeline:
	def __init__(self, input, runMode=True):
		self.input = input
		self.latestOutput = input
		self.executedModules = []
		self.in2out = generalUtils.in2out
		self.funIn2out = generalUtils.funIn2out
		self.in2tempOut = generalUtils.in2tempOut
		self.out2log = generalUtils.out2log
		self.runGet = generalUtils.runGet
		self.runMode = runMode
		if '--mock' in sys.argv:
			self.runMode = False

	def run(self, codeList, runFlag=True):
		code = generalUtils.list2gappedString(codeList)
		code = code.replace('"', '\\"')
		code = code.replace('(', '\\(')
		code = code.replace(')', '\\)')
		allStringList = generalUtils.list2allStringList(codeList)
		if runFlag:
			print('-->\t' + code)
			if self.runMode:
				failedHere = os.system(code)
				if failedHere:
					raise ValueError('we cannot execute the code: ' + code)
			# else:
			# 	print("gave up running the code, because the command is not given in the default 'RUN' mode.")
		else:
			print('X\t' + code)

	def runM(self, codeDict, runFlag=True):
		# codeDict = {'codeList': ['code', 'input', 'output'],
		# 			'inputs': inputs,
		# 			'outputF': ['.in', '.out']
		# }
		codeList = codeDict['codeList']
		inputs = codeDict['inputs']
		def input2output(input):
			if 'outputF' in codeDict.keys():
				outputF = codeDict['outputF']
				return self.in2out(input, outputF[0], outputF[1])
			else:
				outputFunction = codeDict['outputFunction']
				return outputFunction(input)

		outputs = []
		for input in inputs:
			output = input2output(input)
			outputs.append(output)
			argPlacedCodeList = [input if x=='#IN' else x for x in codeList]
			argPlacedCodeList = [output if x=='#OUT' else x for x in argPlacedCodeList]
			self.run(argPlacedCodeList, runFlag)
		return outputs

	def internalRun(self, function, arguments, runFlag=True):
		functionName = function.__name__
		if runFlag:
			print('i->\t' + functionName)
			if self.runMode:
				return function(*arguments)
		else:
			print('iX\t' + functionName)
			return False

	def checkInput(self, pairedEnd=False):
		expectedExtension = '.fastq'
		if pairedEnd:
			self.input1 = self.input
			self.input2 = self.in2out(self.input, '.1' + expectedExtension, '.2' + expectedExtension)
			expectedExtension = '.1' + expectedExtension
		if self.input.endswith(expectedExtension):
			return True
		else:
			sys.exit('The format is not in the expected extension ' + expectedExtension)

	def checkPairedEndInput(self):
		return self.checkInput(True)

	def checkDependency_(self, flags):
		for flag in flags:
			if getattr(self, flag) != True:
				sys.exit('dependecy ' + flag + ' is not completed, exiting...')
		return True

	def bed2nonRedundantBed(self, runFlag=True):
		self.checkDependency_(['sorted'])
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.NR.bed')
		singleCodeList = [
			'bed2nonRedundantBed.py',
			input,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def sortBam(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bam', '.sorted.bam')
		log = self.out2log(output)
		temp = '/tmp/' + input.split('/')[-1] + '.temp'
		singleCodeList = [
			'samtools',
			'sort',
			'-T', temp, # temporary file
			'-o', output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.sorted = True
		self.latestOutput = output
		return self

	def sam2bam(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.sam', '.bam')
		log = self.out2log(output)
		singleCodeList = [
			'samtools',
			'view',
			'-bS',
			'-o',
			output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bam2bed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bam', '.bed')
		log = self.out2log(output)
		singleCodeList = [
			'bedtools',
			'bamtobed',
			'-i', input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		self.latestBed = output
		return self

	def bed2firstReadOnlyBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.fro.bed')
		log = self.out2log(output)
		singleCodeList = [
			'bed2firstReadOnlyBed.py',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def sam2bed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.sam', '.bed')
		log = self.out2log(output)
		singleCodeList = [
			'sam2bed',
			input,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def sortBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.sorted.bed')
		self.sorted = True
		log = self.out2log(output)
		singleCodeList = [
			'sortBed',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self
