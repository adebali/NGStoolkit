import os
import sys
import subprocess
from seqPipeline import SeqPipeline
from glob import glob

# Constants
kHomeDir = os.path.join(os.path.expanduser("~"), "ogun")
kChromHmmDir = os.path.join(kHomeDir, "chromHMM")

chromHMMfiles = glob(os.path.join(kHomeDir, 'chromHMM', '*_chromHMM'))

# General utils

class DamageSeqSingleReadPipeline(SeqPipeline):
	'''DamageSeq Single End Pipeline class'''

	# Methods
	##########################################################

	def cutadapt(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.fastq', '.cutadapt.fastq')
		log = self.out2log(output)
		adapter = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
		singleCodeList = [
			'cutadapt',
			'--discard', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
			'-g', adapter, # The adapter is located at the 5' end
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
		output = self.in2out(input, '.fastq', '.bowtie_hg19.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t /proj/seq/data/bowtie/hg19/hg19', # hg19 Human Reference Genome, Feb. 2009, Genome Reference Consortium GRCh37
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAC
			'--phred33-quals', # Depends on the sequencing platform
			'-S', # Output in SAM format
			'-m 4', # Do not report the reads that are mapped on to more than 4 genomic locations
			'--seed 123', # Randomization parameter in bowtie,
			input,
			output,
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
			'-c', '1,2,3' # columns
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def slopBed(self, runFlag=True):
		input = self.latestOutput
		slopL = 3
		slopR = -49
		output = self.in2out(input, 'bed', 'slopL' + str(slopL) + 'R' + str(slopR) + '.bed')
		genome = kHomeDir + '/seq/hg19/hg19.bed'
		log = self.out2log(output)
		singleCodeList = [
			'bedtools',
			'slop',
			'-i', input,
			'-g', genome, # Chomosomal lengths, needed in case region shift goes beyond the chromosomal limits.
			'-l', slopL, # The number of base pairs to subtract from the start coordinate.
			'-r', slopR, # The number of base pairs to add to the end coordinate.
			'-s', # Take strand into account
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		self.finalBed = output
		return self

	def bed2fasta(self, runFlag=False):
		input = self.finalBed
		output = self.in2out(input, '.bed', '.bed.fa')
		referenceFasta = '/proj/seq/data/hg19/sam_index/hg19.fasta'
		singleCodeList = [
			'bedtools',
			'getfasta',
			'-fi', referenceFasta,
			'-bed', input,
			'-fo', output,
			'-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
		]
		self.run(singleCodeList, runFlag)
		self.fasta = output
		return self

	def fa2nucleotideAbundanceTable(self, runFlag=False):
		input = self.fasta
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
				coverageInput,
				normalizedOutput,
				multiplyFactor
			]
			self.run(singleCodeList, runFlag)
			#########################################

		return self

input = sys.argv[1]
pipeline = DamageSeqSingleReadPipeline(input)

pipeline\
	.cutadapt(False)\
	.bowtie(False)\
	.sam2bam(False)\
	.sortBam(False)\
	.bam2bed(False)\
	.bed2nonRedundantBed(False)\
	.slopBed(False)\
		.bed2fasta(False)\
			.fa2nucleotideAbundanceTable(True)\
		.coverageChromHMM(False)
