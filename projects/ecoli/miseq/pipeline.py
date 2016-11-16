import os
import sys
import subprocess
import generalUtils
import fasta
import bed
from seqPipeline import SeqPipeline
from glob import glob


def name2treatment(name):	
	treatments = {
		"WT": "WT",
		"UVRD": "UvrD",
		"MFD": "Mfd",
		"PHR": "Phr"
	}
	for key in treatments.keys():
		if name.startswith(key):
			return treatments[key]
	raise ValueError('no treatment found')

tempDirectory = ' /netscr/adebali'

def fileName2adapter(fileName):
	illumina_small_RNA_adaptor = '/nas02/home/a/d/adebali/ogun/seq/illumina_small_RNA_adaptor.fa'
	## Assumes index sequence is placed between two underscores (_) in the file name such as UVRD10_TAGCTT_L007_R1_001.fastq
	adapter = ''
	indexSeq = fileName.split('_')[1]
	adapters = fasta.fasta(illumina_small_RNA_adaptor)
	adapterDict = adapters.read()
	adapterFound = False
	for header in adapterDict.keys():
		sequence = adapterDict[header]
		if sequence[33:33+len(indexSeq)] == indexSeq:
			if adapterFound:
				sys.exit('More than one adapters matching the criteria. Exiting...')
			adapter = sequence
			adapterFound = True
	if adapterFound:
		return adapter
	else:
		sys.exit('Adapter is not found: ' + indexSeq)


class XRseqPipeline(SeqPipeline):
	'''XRseq Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkInput()
		self.ecoliReferenceRoot = '/nas02/home/a/d/adebali/ncbi/ecoli'
		self.referenceVersion = 'NC_000913.2'
		self.referenceRoot = os.path.join(self.ecoliReferenceRoot, self.referenceVersion)
		self.referenceBowtieIndex = os.path.join(self.referenceRoot, self.referenceVersion)
		self.referenceFasta = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.fa'))
		self.referenceGenesBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.genes.bed'))
		self.referenceGFF = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.gff3'))
		# self.referenceWindowsBedFile = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.bed'))
		# self.fiftyWindowBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win50.bed'))
		self.referenceGenomeSize = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chrom.sizes'))
		# self.twohundredWindowPlusBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.Plus.bed'))
		# self.twohundredWindowMinusBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.win200.Minus.bed'))
		
		self.baseFileName = os.path.basename(self.input)
		self.adapter = fileName2adapter(self.baseFileName)
		self.treatment = name2treatment(self.baseFileName)
		# self.latestOutputs = []


	# Class Specific Methods
	##########################################################

	def cutadapt(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'cutadapt.fastq')
		log = self.out2log(output)
		singleCodeList = [
			'cutadapt',
			'-a', self.adapter, # The adapter is located at the 3' end
			'-m', 13, # Minimum length
			'-M', 13, # Maximum length
			'-o', output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.executedModules.append(self.cutadapt.__name__)
		self.latestOutput = output
		return self

	def fastq2lengthDistribution(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.fastq', '.fastq.lengthDist.txt')
		singleCodeList = [
			'fastq2lengthDistribution.py',
			'-i', input,
			'-o', output,
			'-l', 50
		]
		self.run(singleCodeList, runFlag)
		return self

	def fastqSampling(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.fastq', '.2900K.fastq')
		singleCodeList = [
			'fastq2sampling.py',
			'-i', input,
			'-o', output,
			'-c', 2900000
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def bowtie(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'bowtie_ecoli.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t', self.referenceBowtieIndex,
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAC
			'--phred33-quals', # Depends on the sequencing platform
			'-S', # Output in SAM format
			'--all', 
			'--strata',
			'--best',
			'--seed 123', # Randomization parameter in bowtie,
			input,
			output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def splitBamByStrand(self, runFlag=True):
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

	def bed2bam(self, runFlag=True):
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

	def bam2bai(self, runFlag=True):
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
		self.runM(codeDict)
		return self


	# def bed2nonRedundantBed(self, runFlag=True):
	# 	self.checkDependency_(['sorted'])
	# 	input = self.latestOutput
	# 	output = self.in2out(input, '.bed', '.NR.bed')
	# 	singleCodeList = [
	# 		'bed2nonRedundantBed.py',
	# 		'-i', input,
	# 		'-o', output,
	# 		'-c', 3
	# 	]
	# 	self.run(singleCodeList, runFlag)
	# 	self.latestOutput = output
	# 	return self

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
		self.sortedBam = output
		return self

	def randomSelection(self, runFlag=True):
		input = self.latestOutput
		count = 2400000
		output = self.in2out(input, '.bed', '.2M.bed')
		codeList = [
			'bed2sampling.py',
			'-i', input,
			'-c', count,
			'-o', output
		]
		self.run(codeList, runFlag)
		self.latestOutput = output
		return self


	def splitBedByStrand(self, runFlag=True):
		input = self.latestOutput
		outputPlus = self.in2out(input, '.bed', '.Plus.bed')
		outputMinus = self.in2out(input, '.bed', '.Minus.bed')
		singleCodeListPlus = [
			'grep',
			'-P', '\'\+$\'',
			input,
			'>', outputPlus
		]
		self.run(singleCodeListPlus, runFlag)
		self.latestOutput = outputPlus
		self.latestOutputPlus = outputPlus

		singleCodeListMinus = [
			'grep',
			'-P', '\'\-$\'',
			input,
			'>', outputMinus
		]
		self.run(singleCodeListMinus, runFlag)
		self.latestOutputMinus = outputMinus
		self.firstPlusBed = self.latestOutputPlus
		self.firstMinusBed = self.latestOutputMinus
		self.latestOutputs = self.firstBeds = [outputPlus, outputMinus]
		self.treatments = []
		self.treatments.append(self.treatment + '-Plus')
		self.treatments.append(self.treatment + '-Minus')
		return self

	def windowCoverage(self, runFlag=True):
		inputs = [self.latestOutputPlus, self.latestOutputMinus]
		outputPlus = self.in2out(self.latestOutputPlus, '.bed', '.covCnt.bed')
		outputMinus = self.in2out(self.latestOutputMinus, '.bed', '.covCnt.bed')
		outputs = [outputPlus, outputMinus]

		for i in range(len(inputs)):
			singleCodeList = [
					'bedtools',
					'coverage',
					'-counts',
					'-a', self.referenceWindowsBedFile,
					'-b', inputs[i], # Read file, to get the read counts overlapping with the primary genomic locations
					'>', outputs[i]
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
					'-t', 4
				]
			self.run(singleCodeList, runFlag)
		self.latestOutput = outputs[0]
		self.latestOutputPlus = outputPlus
		self.latestOutputMinus = outputMinus
		return self

	def bed2fasta(self, runFlag=False):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.bed.fa')
		referenceFasta = self.referenceFasta
		singleCodeList = [
			'bedtools',
			'getfasta',
			'-fi', referenceFasta,
			'-bed', input,
			'-fo', output,
			'-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
		]
		self.run(singleCodeList, runFlag)
		firstFasta = fasta.fasta(output)
		# self.totalReadNumber = self.internalRun(firstFasta.getSequenceCount, [], runFlag)
		self.latestFasta = output
		self.latestOutput = output
		return self

	def fasta2bed_GetTT(self, runFlag=False):
		input = self.latestOutput
		output = self.in2out(input, '.fa', '.TT.bed')
		singleCodeList = [
			'fa2bedByCriteria.py',
			'-i', input,
			'-o', output,
			'-c', 'sequence[7:9]=="TT"' # Get sequences having TT dimer at the positions 8 and 9.
		]
		self.run(singleCodeList, runFlag)
		bedObject = bed.bed(output)
		self.totalReadNumber = self.internalRun(bedObject.getHitNum, [], runFlag)
		self.latestBed = output
		self.latestOutput = output
		return self


	def normalizeBedToBedGraph(self, runFlag=False):
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


	def fa2nucleotideAbundanceTable(self, runFlag=False):
		input = self.latestFasta
		output = self.in2out(input, '.fa', '.fa.nucAbu.csv')
		nucleotideOrder = 'TCGA'
		singleCodeList = [
			'fa2nucleotideAbundanceTable.py',
			'-i', input,
			'-o', output,
			'-n', nucleotideOrder
		]
		self.run(singleCodeList, runFlag)
		self.nucleotideAbundance = output
		return self

	# def normalizeByMotifNumber(self, runFlag=False):
	# 	input = self.latestOutput
	# 	output = self.in2out(input, '.bed', '.nTT.bed')
	# 	singleCodeList = [
	# 		'bedCount2normalizedByMotifNumber.py'
	# 		'-i', input,
	# 		'-f', self.referenceFasta,
	# 		'-m', 'TT',
	# 		'-c', 4,
	# 		'-s', 3,
	# 		'perNmotif', 13,
	# 		'-o', output
	# 	]
	# 	self.run(singleCodeList, runFlag)
	# 	self.latestOutput = output
	# 	return self

	def fa2dimerAbundanceTable(self, runFlag=False):
		input = self.latestFasta
		output = self.in2out(input, '.fa', '.fa.dimerAbu.csv')
		singleCodeList = [
			'fa2kmerAbundanceTable.py',
			'-i', input,
			'-o', output,
			'-k', 2,
			'-l', 13,
			'--percentage'
		]
		self.run(singleCodeList, runFlag)
		return self

	def plotNucleotideAbundance(self, runFlag=False):
		input = self.nucleotideAbundance
		nucleotideOrder = 'TCGA'
		singleCodeList = [
			'plotNucleotideAbundance.r',
			input,
			self.treatment
		]
		self.run(singleCodeList, runFlag)
		return self

	def genesCoverageCount(self, runFlag=True):
		codeDict = {
			'inputs' : self.latestOutputs,
			'outputF': ['.bed', '.geneCov.bed'],
			'codeList' : [
				'bedtools',
				'coverage',
				'-counts',
				'-a', self.referenceGenesBed,
				'-b', '#IN',
				'>', '#OUT'
			]
		}
		self.bedCountFiles = self.runM(codeDict, runFlag)
		return self

	def windowCoverageCount(self, runFlag=True):
		inputs = self.firstBeds
		covereageFiles = [self.twohundredWindowPlusBed, self.twohundredWindowMinusBed]
		outputs = []
		i = 0
		for input in inputs:
			output = self.in2out(input, '.bed', '.50win.bed')
			outputs.append(output)
			coverageFile = coverageFiles[i]
			singleCodeList = [
				'bedtools',
				'coverage',
				'-counts',
				'-a', covereageFile,
				'-b', input,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self
	
	def normalizeByGeneLength(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.nGene.bed')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedCount.py',
				'-i', input,
				'-c', 5, # Count tab number
				# '-c', 4, # Count tab number
				'-l', 1000, # Value will be per 1000 nucleotides
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def normalizeByTT(self, runFlag=True):
		inputs = self.bedCountFiles
		strands = ['+', '-']
		outputs = []
		i = 0
		for input in inputs:
			strand = strands[i]
			i += 1
			output = self.in2out(input, '.bed', '.nTT.bed')
			TTcountOutput = self.in2out(input, '.bed', '.bed.TTcount.txt')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedByMotifNumber.py',
				'-i', input,
				'-f', self.referenceFasta,
				'-m', 'TT',
				'-s', strand,
				'-perNmotif', 100,
				'-c', 5, # Count tab number
				'--conPrint',
				'-o', output,
				'>', TTcountOutput
			]
			self.run(singleCodeList, runFlag)
		self.bedCountFiles = outputs
		return self


	def normalizeByReadNumber(self, runFlag=True):
		codeDict = {
			'inputs': self.bedCountFiles,
			'codeList': [
				'bedCount2normalizedByReadNumber.py',
				'-i', '#IN',
				'-c', 5,
				'-readNum', self.totalReadNumber, # Count tab number
				'-perNumReads', 1000000, # Value will be per million reads
				'-o', '#OUT'
			],
			'outputF': ['.bed', '.nRead.bed']
		}
		self.bedCountFiles = self.runM(codeDict, runFlag)
		return self

	def bed2bedgraph(self, runFlag=True):
		inputs = [self.firstPlusBed, self.firstMinusBed]
		outputs = []

		for input in inputs:
			output = self.in2out(input, '.bed', '.bdg')
			outputs.append(output)
			singleCodeList = [
				'bedtools',
				'genomecov',
				'-i', input,
				'-g', self.referenceGenomeSize,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def addHeaderToBedGraph(self, runFlag=True):
		inputs = self.latestOutputs
		treatments = ['Plus', 'Minus']
		outputs = []
		i = 0
		for input in inputs:
			output = self.in2out(input, '.bdg', '.head.bdg')
			outputs.append(output)
			treatment = self.treatment + '_' + treatments[i]
			i += 1
			trackLine = 'track type=bedGraph name="' + treatment + '"'
			singleCodeList = [
				'echo', trackLine, '>', output,
				'&&',
				'tail -n +2', input, '>>', output
				]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self



	def bed2wig(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.wig')
			outputs.append(output)
			singleCodeList = [
				'igvtools',
				'count',
				input,
				output,
				self.referenceGenomeSize
				]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def addHeaderToWig(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		i = 0
		for input in inputs:
			output = self.in2out(input, '.wig', '.head.wig')
			outputs.append(output)
			singleCodeList = [
				'echo "track type=wiggle_0 name=' + self.treatments[i] + '"', '>', output,
				'&&',
				'tail -n +2', input, '>>', output
			]
			self.run(singleCodeList, runFlag)
			i += 1
		self.latestOutputs = outputs
		return self

	def lateralConcatanate(self, runFlag=True):
		inputs = self.latestOutputs
		output = self.in2out(inputs[0], '.bed', '.mrg.bed')
		singleCodeList = [
			'paste',
			inputs[0],
			inputs[1],
			'| cut -f 1-4,8',
			#'| cut -f 1-5,10', for gene coverage
			'>', output
		]
		self.latestOutput = output
		self.run(singleCodeList, runFlag)
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

	def templateMinusCoding(self,runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.tmc.bed')
		singleCodeList = [
			'geneCoveragePasted2TempMinusCod.py',
			input,
			'>', output
		]
		self.latestOutput = output
		self.run(singleCodeList, runFlag)
		return self

	def bedGraph2bigWig(self, runFlag=True):
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

	def bed2bigBed(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.bb')
			outputs.append(output)
			singleCodeList = [
				'bedToBigBed',
				input,
				self.referenceGenomeSize,
				output
				]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self


	def unixSortBed(self, runFlag=True):
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.sorted.bed')
			outputs.append(output)
			self.sorted = True		
			log = self.out2log(output)
			singleCodeList = [
				'sort',
				input,
				'-k1,1',
				'-k2,2n',
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self
	

class XRseqPipeline_seqLengths(XRseqPipeline):
	'''XRseq Pipeline class'''

	def cutadapt(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'cutadaptAll.fastq')
		log = self.out2log(output)
		singleCodeList = [
			'cutadapt',
			'-a', self.adapter, # The adapter is located at the 3' end
			'-m', 9,
			'-M', 15,
			'-o', output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.executedModules.append(self.cutadapt.__name__)
		self.latestOutput = output
		return self


	def separateByLengthAndWriteKmerAbundance(self, runFlag=True):
		input = self.latestFasta
		output = self.in2out(input, '.fa', '.nucAbu.fa')
		log = self.out2log(output)
		singleCodeList = [
			'fa2lengthSeparatedKmerAbundace.py',
			'-i', input,
			'-o', output,
			'-k', 1,
			'-r', '9-15'
		]
		self.run(singleCodeList, runFlag)
		return self


input = sys.argv[1]

pipeline = XRseqPipeline(input)

pipeline\
	.cutadapt(False)\
		.fastq2lengthDistribution(False)\
	.fastqSampling(False)\
	.bowtie(False)\
	.sam2bam(False)\
	.sortBam(True)\
	.bam2bed(True)\
	.bed2fasta(True)\
		.fa2nucleotideAbundanceTable(True)\
			.plotNucleotideAbundance(True)\
		.fa2dimerAbundanceTable(True)\
	.fasta2bed_GetTT(True)\
	.splitBedByStrand(True)\
		.normalizeBedToBedGraph(True)\
			.bedGraph2bigWig(True)\
			.sendBigWigFiles(True)\
		.genesCoverageCount(True)\
			.normalizeByReadNumber(True)\
	.bed2bam(True)\
		.bam2bai(True)\


# 1224167
		# .normalizeByTT(False)\

		# .bed2bedgraph(False)\
		# .addHeaderToBedGraph(False)\

		# .bed2wig(True)\
		# .addHeaderToWig(True)\
	# .lateralConcatanate(True)\
	# .normalizeByGeneLength(False)\
	#.windowCoverageCount(False)\
	#.makeTemplateAndCodingStrands(True)
	##.templateMinusCoding(True)
	# .windowCoverage(False)\
	# .bedCount2percentageByMax(False)




# pipeline_spearateLengths = XRseqPipeline_seqLengths(input)

# pipeline_spearateLengths\
# 	.cutadapt(False)\
# 	.bowtie(False)\
# 	.sam2bam(False)\
# 	.sortBam(False)\
# 	.bam2bed(False)\
# 		.bed2fasta(False)\
# 			.separateByLengthAndWriteKmerAbundance(False)
			