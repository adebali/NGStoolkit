import os
import sys
import subprocess
import generalUtils
import fasta
import bed
from seqPipeline import SeqPipeline
from glob import glob


treatments = {
	"WT04_TGACC_L001_R1_001.fastq": "WT",
	"UVRD10_TAGCTT_L007_R1_001.fastq": "uvrD",
	"MFD09_GATCAG_L007_R1_001.fastq": "Mfd",
	"PHR11_GGCTAC_L007_R1_001.fastq": "phr"
}

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
		self.referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'
		self.referenceFasta = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.fa')
		self.referencesGenesBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.genes.bed')
		self.e_coli_K12_windowsBedFile = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed')
		self.referenceGenomeLength = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.genome')
		self.genomeSize = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.chrom.sizes')
		self.fiftyWindowBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.win50.bed')
		self.twohundredWindowPlusBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.win200.Plus.bed')
		self.twohundredWindowMinusBed = generalUtils.file('/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.win200.Minus.bed')
		
		self.baseFileName = os.path.basename(self.input)
		self.adapter = fileName2adapter(self.baseFileName)
		self.treatment = treatments[self.baseFileName]
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
			# '-m 4', # Do not report the reads that are mapped on to more than 4 genomic locations
			'--seed 123', # Randomization parameter in bowtie,
			input,
			output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
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
					'-a', self.e_coli_K12_windowsBedFile,
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
		self.totalReadNumber = firstFasta.getSequenceCount()
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
		self.totalReadNumber = bedObject.getHitNum()
		self.latestBed = output
		self.latestOutput = output
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
		inputs = self.firstBeds
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
		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bed', '.nTT.bed')
			TTcountOutput = self.in2out(input, '.bed', '.bed.TTcount.txt')
			outputs.append(output)
			singleCodeList = [
				'bedCount2normalizedByMotifNumber.py',
				'-i', input,
				'-f', self.referenceFasta,
				'-m', 'TT',
				'-s', 4,
				'-perNmotif', 100,
				'-c', 5, # Count tab number
				'--conPrint',
				'-o', output,
				'>', TTcountOutput
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
				'-c', 5,
				'-readNum', self.totalReadNumber, # Count tab number
				'-perNumReads', 1000000, # Value will be per million reads
				'-o', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
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
				'-g', self.referenceGenomeLength,
				'>', output
			]
			self.run(singleCodeList, runFlag)
		self.latestOutputs = outputs
		return self

	def bedgraph2bigwig(self, runFlag=True):

		inputs = self.latestOutputs
		outputs = []
		for input in inputs:
			output = self.in2out(input, '.bedgraph', '.bw')
			outputs.append(output)
			singleCodeList = [
				'bedtools',
				'genomecov',
				'-i', input,
				'-g', self.referenceGenomeLength,
				'>', output
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
				self.genomeSize
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
	.bowtie(False)\
	.sam2bam(False)\
	.sortBam(False)\
	.bam2bed(False)\
	.bed2fasta(False)\
		.fa2nucleotideAbundanceTable(False)\
			.plotNucleotideAbundance(False)\
		.fa2dimerAbundanceTable(False)\
	.fasta2bed_GetTT(False)\
	.splitBedByStrand(False)\
	.genesCoverageCount(False)\
	.normalizeByTT(True)\
	.normalizeByReadNumber(True)\
		.bed2wig(True)\
		.addHeaderToWig(True)\
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
			