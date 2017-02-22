import os
import sys
import subprocess
import generalUtils
from seqPipeline import SeqPipeline
from glob import glob

# Constants
kHomeDir = os.path.join(os.path.expanduser("~"), "ogun")
kChromHmmDir = os.path.join(kHomeDir, "chromHMM")

chromHMMfiles = glob(os.path.join(kHomeDir, 'chromHMM', '*_chromHMM'))


class XRseqPipeline(SeqPipeline):
	'''XRseq Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkInput()
		self.referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'
		self.referenceFasta = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.fa')
		self.e_coli_K12_windowsBedFile = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed')



	# Class Specific Methods
	##########################################################

	def cutadapt(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'cutadapt.fastq')
		log = self.out2log(output)
		adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'
		singleCodeList = [
			'cutadapt',
			'-a', adapter, # The adapter is located at the 3' end
			'-o', output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.executedModules.append(self.cutadapt.__name__)
		self.latestOutput = output
		return self

	def bowtie(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, 'fastq', 'bowtie_hg19.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t', self.referenceBowtieIndex,
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAC
			'--phred33-quals', # Depends on the sequencing platform
			'-S', # Output in SAM format
			'-n', 2, # No more than 2 mismatches
			'-e', 70, # The sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed E 
			'-m 4', # Do not report the reads that are mapped on to more than 4 genomic locations
			'--seed 123', # Randomization parameter in bowtie,
			input,
			output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def sam2bam(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.sam', '.bam')
		log = self.out2log(output)
		singleCodeList = [
			'samtools',
			'view',
			'-bf', '0x2', #	each segment properly aligned according to the aligner
			'-o',
			output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2nonRedundantBed(self, runFlag=True):
		self.checkDependency_(['sorted'])
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.NR.bed')
		singleCodeList = [
			'bed2nonRedundantBed.py',
			'-i', input,
			'-o', output,
			'-c', 3
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		self.finalBed = output
		return self

	def splitBedByStrand(self, runFlag=True):
		input = self.latestOutput
		outputPlus = self.in2out(input, '.bed', '.Plus.bed')
		outputMinus = self.in2out(input, '.bed', '.Minus.bed')
		singleCodeListPlus = [
			'grep',
			'-P', '"\+$"',
			input,
			'>', outputPlus
		]
		self.run(singleCodeListPlus, runFlag)
		self.latestOutput = outputPlus
		self.latestOutputPlus = outputPlus

		singleCodeListMinus = [
			'grep',
			'-P', '"\-$"',
			input,
			'>', outputMinus
		]
		self.run(singleCodeListMinus, runFlag)
		self.latestOutputMinus = outputMinus
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
		input = self.latestBed
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
		self.latestFasta = output
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
		return self

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

	def coverageChromHMM(self, runFlag=False):
		input = self.finalBed
		for chromHMMfile in chromHMMfiles:
			regionName = chromHMMfile.split('/')[-1].split('_NHFL_')[0]
			output = self.in2out(input, 'bed', 'coverage_' + regionName + '.bed')
			singleCodeList = [
				'bedtools',
				'coverage',
				'-counts',
				'-a', chromHMMfile, # Primary genomic locations
				'-b', input, # Read file, to get the read counts overlapping with the primary genomic locations
				'>', output
			]
			self.run(singleCodeList, runFlag)

			# Normalization between Chromatin States
			coverageInput = output
			normalizedOutput = self.in2out(coverageInput, '.bed', '.norm.bed')
			multiplyFactor = 1000
			singleCodeList = [
				'bedCount2normalizedCount.py',
				'-i', coverageInput,
				'-o', normalizedOutput,
				'-l', multiplyFactor,
				'-c', 7
			]
			self.run(singleCodeList, runFlag)
			#########################################

		return self

	def fa2dimerAbundanceTable(self, runFlag=False):
		input = self.latestOutput
		output = self.in2out(input, '.sam', '.bed')
		singleCodeList = [
			'sam2bed',
			input,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


input = sys.argv[1]

pipeline = XRseqPipeline(input)

pipeline\
	.cutadapt(False)\
	.bowtie(False)\
	.sam2bed(True)\
	.sortBed(True)\
	.bed2nonRedundantBed(True)
	# .coverageChromHMM(False)

