import os
import sys
import subprocess
import generalUtils
from seqPipeline import SeqPipeline
from glob import glob
import bed


class rnaSeq(SeqPipeline):
	'''rna-seq Paired End Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkInput()
		self.basename = os.path.basename(self.input).replace('.fastq','')
		self.dataDir = os.path.dirname(self.input)
		self.ecoliReferenceRoot = '/nas02/home/a/d/adebali/ncbi/ecoli'
		self.referenceVersion = 'NC_000913.2'
		self.referenceRoot = os.path.join(self.ecoliReferenceRoot, self.referenceVersion)
		# self.referenceBowtieIndex = os.path.join(self.referenceRoot, self.referenceVersion)
		self.referenceBowtieIndex = os.path.join(self.referenceRoot, self.referenceVersion) + '_chr'
		self.referenceTophatIndex = self.referenceBowtieIndex
		# self.referenceFasta = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.fa'))
		self.referenceGenesBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.genes.bed'))
		# self.referenceGenesBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.genes.bed'))
		self.TSSbed = generalUtils.file(os.path.join(self.dataDir, 'TSSmap.bed'))
		self.referenceGFF3 = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.gff3'))
		# self.referenceGFF = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.gff'))
		self.referenceGFF = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.gff'))
		# self.referenceWindowsBedFile = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.bed'))
		# self.fiftyWindowBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win50.bed'))
		# self.referenceGenomeSize = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chrom.sizes'))
		self.referenceGenomeSize = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.sizes'))
		self.referenceGenomeSize = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.sizes'))
		# self.twohundredWindowPlusBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.Plus.bed'))
		# self.twohundredWindowMinusBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.Minus.bed'))
		# self.fiftyWindowBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.2/NC_000913.2.win50.bed')



	# Class Specific Methods
	##########################################################

	# def reademptionAlign(self, runFlag=True):
	# 	input = self.input
	# 	output = 

	def trimNonqualifiedNucleotides(self, runFlag=True):
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
		# make sure that index file source (FASTA) starts with >chr followed by the sequence
		input = self.latestOutput
		tophatOutputDir = os.path.join(self.dataDir, 'tophat_' + self.basename)
		output = os.path.normpath(os.path.join(tophatOutputDir, '..', self.basename + '.bam'))
		log = self.out2log(output)
		singleCodeList = [
			'mkdir -p ',
			tophatOutputDir,
			'&&',
			'tophat',
			'--library-type', 'fr-firststrand',
			'-o', tophatOutputDir,
			self.referenceTophatIndex,
			input,
			'&&',
			'cp', os.path.join(tophatOutputDir, 'accepted_hits.bam'), output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def splitBamByStrand(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		outputPlus = self.in2out(input, '.bam', '.Plus.bam')
		outputMinus = self.in2out(input, '.bam', '.Minus.bam')
		outputs = [outputPlus, outputMinus]
		codeList = [
			'samtools',
			'view',
			'-bF', '0x10',
			input,
			'>', outputPlus
		]
		self.run(codeList, runFlag)

		codeList = [
			'samtools',
			'view',
			'-bf', '0x10',
			input,
			'>', outputMinus
		]
		self.run(codeList, runFlag)

		self.latestOutputs = outputs
		return self

	def bam2bai(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bam','.bam.bai'],
			'codeList': [
				'samtools',
				'index',
				'#IN',
				'#OUT'
			]
		}
		self.runM(codeDict, runFlag)
		return self

	# def sortBam(self, runFlag=True):
	# 	input = self.latestOutput
	# 	output = self.in2out(input, '.bam', '.sorted.bam')
	# 	log = self.out2log(output)
	# 	temp = '/tmp/' + input.split('/')[-1] + '.temp'
	# 	singleCodeList = [
	# 		'samtools',
	# 		'sort',
	# 		'-T', temp, # temporary file
	# 		'-o', output,
	# 		input,
	# 		'>', log
	# 	]
	# 	self.run(singleCodeList, runFlag)
	# 	self.sorted = True
	# 	self.latestOutput = output
	# 	self.sortedBam = output
	# 	return self

	def bam2bed(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		outputF = ['.bam', '.bed']
		codeList = [
			'bedtools',
			'bamtobed',
			'-i', '#IN',
			'>', '#OUT'
		]
		codeDict = {'inputs':inputs, 'outputF': outputF, 'codeList': codeList}
		self.latestOutputs = self.runM(codeDict, runFlag)

		# bedObject = bed.bed(output)
		#self.totalReadNumber = bedObject.getHitNum()
		#self.averageReadLength = bedObject.getAverageLength()
		return self

	def splitBedByStrand(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		outputs = [	self.in2out(input, '.bed', '.Plus.bed'),
					self.in2out(input, '.bed', '.Minus.bed')
		]
		strands = ['+', '-']
		i = 0
		for output in outputs:
			singleCodeListPlus = [
				'grep',
				'-P', '\'\\' + strands[i] + '$\'',
				input,
				'>', output
			]
			self.run(singleCodeListPlus, runFlag)
			i += 1

		self.latestOutputs = outputs
		return self

	def bed2bam(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bed','.bam'],
			'codeList': [
				'bedToBam',
				'-i', '#IN',
				'-g', self.referenceGenomeSize,
				'>', '#OUT'
			]
		}
		self.latestOutputs = self.runM(codeDict, runFlag)
		return self
	
	def sortBam(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		temp = '/tmp/' + os.path.basename(inputs[0]) + '.temp'
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bam', '.sorted.bam'],
			'codeList': [
				'samtools',
				'sort',
				'-T', temp, # temporary file
				'-o', '#OUT',
				'#IN'
			]
		}
		self.latestOutputs = self.runM(codeDict, runFlag)
		return self

	def windowCoverage(self, runFlag=True):
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		outputs = []
		# multiplyFactor = 13 / self.averageReadLength
		multiplyFactor = 1
		for input in inputs:
			output = self.in2out(input, '.bed', '.nGene.bed')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedCount.py',
				'-i', input,
				'-c', 5, # Count tab number geneCoverage
				# '-c', 4, # Count tab number
				'-l', 1000, # Value will be per 1000 nucleotides
				'-m', multiplyFactor,
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def normalizeByReadNumber(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.nRead.bed')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedByReadNumber.py',
				'-i', input,
				'-c', 5,
				'-readNum', self.totalReadNumber, # Count tab number
				'-perNumReads', 10000000, # Value will be per million reads
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def intervalCoverageCount(self, runFlag=True):
		if runFlag == "Skip": return self
		coverageType = 'genes'
		coverageReference = self.referenceGenesBed
		inputs = self.latestOutputs
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bed', '.' + coverageType + 'Cov.bed'],
			'codeList': [
				'bedtools',
				'coverage',
				'-counts',
				'-a', coverageReference,
				'-b', '#IN',
				'>', '#OUT'
			]
		}
		self.latestOutputs = self.runM(codeDict, runFlag)
		return self

	def tssCoverageCount(self, runFlag=True):
		if runFlag == "Skip": return self
		coverageType = 'Tss'
		coverageReference = self.TSSbed
		inputs = self.latestOutputs
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bed', '.' + coverageType + 'Cov.bed'],
			'codeList': [
				'bedtools',
				'coverage',
				'-counts',
				'-s',
				'-a', coverageReference,
				'-b', '#IN',
				'>', '#OUT'
			]
		}
		self.latestOutputs = self.runM(codeDict, runFlag)
		return self

	def windowCoverageCount(self, runFlag=True):
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
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


	def bam2assembledTranscripts(self, runFlag=True):
		if runFlag == "Skip": return self
		# def addCufflinks():
		# 	self.run(['module add cufflinks'], runFlag)

		#addCufflinks()

		def in2out(input):
			return os.path.join(os.path.dirname(input),os.path.basename(input)).replace('.','_') + '_cufflinks'
		
		inputs = self.latestOutputs
		for input in inputs:
			outputDir = in2out(input)
			codeList = ['mkdir -p', outputDir]
			self.run(codeList, runFlag)

		codeDict = {
			'inputs': inputs,
			'outputFunction': in2out,
			'codeList' : [
				'source /nas02/apps/Modules/default/init/bash',
				'&& module unload samtools',
				'&& module load cufflinks',
				'cufflinks',
				'#IN',
				'--library-type', 'fr-firststrand',
				'-r', self.referenceGFF,
				'-o', '#OUT'
			]
		}
		self.latestOutputs = self.runM(codeDict, runFlag)

	def normalizeBedToBedGraph(self, runFlag=False):
		if runFlag == "Skip": return self
		inputs = self.latestOutputs
		codeDict = {
			'inputs': self.latestOutputs,
			'outputF': ['bed', 'n1Mread.bg'],
			'codeList': [
				'bedtools', 'genomecov',
				'-i', '#IN',
				'-bg',
				'-scale', 1000000,
				'-g', self.referenceGenomeSize,
				'>', '#OUT'
			]
		}
		self.latestBedGraphs = self.runM(codeDict, runFlag)
		return self

	def bedGraph2bigWig(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestBedGraphs
		codeDict = {
			'inputs': inputs,
			'outputF': ['.bg', '.bw'],
			'codeList': [
				'bedGraphToBigWig',
				'#IN',
				self.referenceGenomeSize,
				'#OUT'
			]
		}
		self.latestBigWigFiles = self.runM(codeDict, runFlag)
		return self

	def sendBigWigFiles(self, runFlag=True):
		if runFlag == "Skip": return self
		inputs = self.latestBigWigFiles
		codeDict = {
			'inputs': inputs,
			'outputF': ['NaN','Nan'],
			'codeList': [
				'scp',
				'#IN',
				'ftp:ftp/ecoli'
			]
		}
		self.runM(codeDict, runFlag)
		return self


input = sys.argv[1]

pipeline = rnaSeq(input)

pipeline\
	.trimNonqualifiedNucleotides(False)\
	.removePolyAtail(False)\
	.align(False)\
	.splitBamByStrand(False)\
	.sortBam(False)\
		.bam2bai(False)\
	.bam2bed(False)\
		.normalizeBedToBedGraph(False)\
		.bedGraph2bigWig(False)\
		.sendBigWigFiles(False)\
	.intervalCoverageCount(False)\
	.normalizeByGeneLengthAndAverageReadLength(False)

	# .tssCoverageCount(True)\

	# .bam2assembledTranscripts(True)
	
	# .splitBedByStrand(False)\
	# .bed2bam(False)\
	# .normalizeByReadNumber(False)\
	# .lateralPaste(False)\
	# .windowCoverageCount(False)\
	#.makeTemplateAndCodingStrands(True)

# pipeline.bam2assembledTranscripts(True)

	# .sortByCount(False)
	# .windowCoverage(True)\
	# .bedCount2percentageByMax(True)
