import os
import sys
import subprocess
import generalUtils
from seqPipeline import SeqPipeline
from glob import glob
import bed

referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'

class rnaSeq(SeqPipeline):
	'''rna-seq Paired End Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkInput()
		self.referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'
		self.referenceFasta = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.fa')
		self.e_coli_K12_windowsBedFile = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed')
		self.referencesGenesBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.genes.bed')
		self.fiftyWindowBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.win50.bed')


	# Class Specific Methods
	##########################################################

	def trimNonqualifiedNucleotides(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'trimmed.fastq')
		log = self.out2log(output)
		singleCodeList = [
			'fastq_quality_trimmer',
			'-t', 20,
			'-Q', 33,
			'-i', input,
			'-o', output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def removePolyAtail(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'noAtail.fastq')
		log = self.out2log(output)
		singleCodeList = [
			'cutadapt',
			'-a', 'A{100}',
			'--minimum-length', 12,
			'-o', output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def align(self, runFlag=True):
		# make sure that index file source (FASTA) starts with >chr followed by the sequence
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'bowtie_ecoli.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t', self.referenceBowtieIndex,
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAQ
			'-S', # Output in SAM format
			'--seed 123', # Randomization parameter in bowtie,
			input,
			output,
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
		bedObject = bed.bed(output)
		self.totalReadNumber = bedObject.getHitNum()
		self.averageReadLength = bedObject.getAverageLength()
		return self

	def splitBedByStrand(self, runFlag=True):
		input = self.latestOutput
		outputs = [	self.in2out(input, '.bed', '.Plus.bed'),
					self.in2out(input, '.bed', '.Minus.bed')
		]
		
		for output in outputs:
			singleCodeListPlus = [
				'grep',
				'-P', '"\+$"',
				input,
				'>', output
			]
			self.run(singleCodeListPlus, runFlag)

		self.latestOutputs = outputs
		return self

	def windowCoverage(self, runFlag=True):
		inputs = self.latestOutputs
		outputMinus = self.in2out(self.latestOutputMinus, '.bed', '.covCnt.bed')
		outputs = [outputPlus, outputMinus]
		e_coli_K12_windowsBedFile = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed'

		for input in inputs:
			output = self.in2out(input, '.bed', '.covCnt')
			singleCodeList = [
					'bedtools',
					'coverage',
					'-counts',
					'-a', e_coli_K12_windowsBedFile,
					'-b', 
					'>', output
				]
			self.run(singleCodeList, runFlag)
		self.latestOutput = outputs[0]
		self.latestOutputPlus = outputPlus
		self.latestOutputMinus = outputMinus
		return self

	def bedCount2percentageByMax(self, runFlag=True):
		inputs = [self.latestOutputPlus, self.latestOutputMinus]
		outputPlus = self.in2out(self.latestOutputPlus, '.bed', '.perc.bed')
		outputMinus = self.in2out(self.latestOutputMinus, '.bed', '.perc.bed')
		outputs = [outputPlus, outputMinus]
		for i in range(len(inputs)):
			singleCodeList = [
					'bedCount2percentageByMax.py',
					'-i', inputs[i],
					'-o', outputs[i],
					'-t', 4,
					'-m', 100
				]
			self.run(singleCodeList, runFlag)
		self.latestOutput = outputs[0]
		self.latestOutputPlus = outputPlus
		self.latestOutputMinus = outputMinus
		return self

	def normalizeByGeneLengthAndAverageReadLength(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		multiplyFactor = 13 / self.averageReadLength
		for input in inputs:
			output = self.in2out(input, '.bed', '.nGene.bed')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedCount.py',
				'-i', input,
				# '-c', 5, # Count tab number geneCoverage
				'-c', 4, # Count tab number
				'-l', 1000, # Value will be per 1000 nucleotides
				'-m', multiplyFactor,
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def normalizeByReadNumber(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.nRead.bed')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedByReadNumber.py',
				'-i', input,
				'-c', 4,
				'-readNum', self.totalReadNumber, # Count tab number
				'-perNumReads', 1000000, # Value will be per million reads
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def genesCoverageCount(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.genesCov.bed')
			outputs.append(output)
			singleCodeList = [
				'bedtools',
				'coverage',
				'-counts',
				'-a', self.referencesGenesBed,
				'-b', input,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def windowCoverageCount(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.50win.bed')
			outputs.append(output)
			singleCodeList = [
				'bedtools',
				'coverage',
				'-counts',
				'-a', self.fiftyWindowBed,
				'-b', input,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def lateralPaste(self, runFlag=True):
		inputs = self.latestOutputs
		output = self.in2out(inputs[0], '.bed', '.mrg.bed')
		singleCodeList = [
			'paste',
			inputs[0],
			inputs[1],
			# '| cut -f 1-5,10',
			'| cut -f 1-4,8',
			'>', output
		]
		self.latestOutput = output
		self.run(singleCodeList, runFlag)
		return self

	def sortByCount(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.cSorted.bed')
			outputs.append(output)
			singleCodeList = [
				'sort',
				'-r',
				'-k5,5n',
				input,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def makeTemplateAndCodingStrands(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.TempCode.bed')
		valueColNo = 3
		value = '+'
		col1no = 4
		col2no = 5
		separator = '\t'
		subArguments = [valueColNo, value, col1no, col2no, separator]
		arguments = [
			input,
			output,
			generalUtils.line2switchColumnsBasedOnValue,
			subArguments
		]
		self.latestOutput = output
		self.internalRun(generalUtils.lineBasedFileOperation, arguments, runFlag)
		return self


input = sys.argv[1]

pipeline = rnaSeq(input)

pipeline\
	.trimNonqualifiedNucleotides(False)\
	.removePolyAtail(False)\
	.align(False)\
	.sam2bam(False)\
	.sortBam(False)\
	.bam2bed(False)\
	.splitBedByStrand(False)\
	.windowCoverageCount(False)\
	.normalizeByGeneLengthAndAverageReadLength(True)\
	.normalizeByReadNumber(True)\
	.lateralPaste(True)\
	#.genesCoverageCount(False)\
	#.makeTemplateAndCodingStrands(True)


	# .sortByCount(False)
	# .windowCoverage(True)\
	# .bedCount2percentageByMax(True)
