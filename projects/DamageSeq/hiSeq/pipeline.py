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

class DamageSeqPairedEndPipeline(SeqPipeline):
	'''DamageSeq Paired End Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkPairedEndInput()

	# Methods
	##########################################################

	def cutadapt(self, runFlag=True):
		input = self.latestOutput
		output1 = self.in2out(self.input1, '1.fastq', 'cutadapt.1.fastq')
		output2 = self.in2out(self.input2, '2.fastq', 'cutadapt.2.fastq')
		log = self.out2log(output1)
		adapter = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'

		singleCodeList = [
			'cutadapt',
			'--discard-trimmed', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
			'-g', adapter, # The adapter is located at the 5' end
			'-o', output1,
			'-p', output2,
			self.input1,
			self.input2,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.executedModules.append(self.cutadapt.__name__)
		self.latestOutput = output1
		self.latestOutput1 = output1
		self.latestOutput2 = output2
		return self

	def bowtie(self, runFlag=True):
		input1 = self.latestOutput1
		input2 = self.latestOutput2
		output = self.in2out(self.latestOutput, '.fastq', '.bowtie_hg19.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t', '/proj/seq/data/bowtie/hg19/hg19', # hg19 Human Reference Genome, Feb. 2009, Genome Reference Consortium GRCh37
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAC
			'--phred33-quals', # Depends on the sequencing platform
			'-S', # Output in SAM format
			'-m', 4, # Do not report the reads that are mapped on to more than 4 genomic locations
			'-X', 1000,
			'--seed', 123, # Randomization parameter in bowtie,
			'-1', input1,
			'-2', input2,
			output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def slopBed(self, runFlag=True):
		input = self.latestOutput
		slopB = 3
		output = self.in2out(input, 'bed', 'slopB' + str(slopB) + '.bed')
		genome = kHomeDir + '/seq/hg19/hg19.bed'
		log = self.out2log(output)
		singleCodeList = [
			'bedtools',
			'slop',
			'-i', input,
			'-g', genome, # Chomosomal lengths, needed in case region shift goes beyond the chromosomal limits.
			'-b', slopB,
			'-s', # Take strand into account
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
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

	def fa2dimerAbundanceTable(self, runFlag=False):
		input = self.fasta
		output = self.in2out(input, '.fa', '.fa.dimerAbu.csv')
		singleCodeList = [
			'fa2kmerAbundanceTable.py',
			'-i', input,
			'-o', output,
			'-k', 2,
			'--percentage'
		]
		self.run(singleCodeList, runFlag)
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

	def bam2bed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bam', '.bed')
		log = self.out2log(output)
		singleCodeList = [
			'bedtools',
			'bamtobed',
			'-bedpe',
			'-mate1',
			'-i', input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2sortedBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.sorted.bed')
		log = self.out2log(output)
		singleCodeList = [
			'sortBed',
			'-i', input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		self.sorted=True
		return self

	def bed2nonRedundantBed(self, runFlag=True):
		self.checkDependency_(['sorted'])
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.NR.bed')
		singleCodeList = [
			'bed2nonRedundantBed.py',
			'-i', input,
			'-o', output,
			'-c', '1,2,3,4,5,6,7,9' # columns
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2uniquelySortedBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.uSorted.bed')
		singleCodeList = [
			'sort',
			'-u',
			'-k1,1',
			'-k2,2',
			'-k3,3',
			'-k4,4',
			'-k5,5',
			'-k6,6',
			#'-k9,9',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def bedpe2bed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.sf.bed')
		singleCodeList = [
			'bedpe2bed.py',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2fixedRangeBed(self, runFlag=True):
		input = self.latestOutput
		side = 'left'
		sideInput = side[0]
		fixedRange = 10
		output = self.in2out(input, '.bed', '.fr' + sideInput.upper() + str(fixedRange) + '.bed')
		singleCodeList = [
			'bed2fixedRangeBed.py',
			'-i', input,
			'-s', sideInput,
			'-l', fixedRange,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		self.finalBed = output
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

#bsub -M 32 'grep -Pv "\t0" iterativeCov.4.bed | sort -u -k1,1 -k2,2 -k3,3 >escapingUniqueDamages.bed'

input = sys.argv[1]
pipeline = DamageSeqPairedEndPipeline(input)

pipeline\
	.cutadapt(False)\
	.bowtie(False)\
	.sam2bam(False)\
	.bam2bed(False)\
	.bed2uniquelySortedBed(False)\
	.bedpe2bed(False)\
	.slopBed(False)\
	.bed2fixedRangeBed(False)\
		.bed2fasta(False)\
			.fa2nucleotideAbundanceTable(False)\
			.fa2dimerAbundanceTable(True)\
		.coverageChromHMM(False)
