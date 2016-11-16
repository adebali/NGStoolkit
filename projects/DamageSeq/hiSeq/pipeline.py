import os
import sys
import subprocess
from seqPipeline import SeqPipeline
from glob import glob
import bed
import generalUtils

# Constants
kHomeDir = os.path.join(os.path.expanduser("~"), "ogun")
kChromHmmDir = os.path.join(kHomeDir, "chromHMM")

chromHMMfiles = generalUtils.sorted_nicely(glob(os.path.join(kHomeDir, 'chromHMM', '*_chromHMM')))

treatments = {
	'NCPD-A6N.1.fastq':  'NHF1_CPD_Control_A',
	'NC-B17N.1.fastq': 'NHF1_CPD_Control_B',	
	'NCPD0HA1N.1.fastq': 'NHF1_CPD_10J_0h_A',
	'NC0HB13N.1.fastq': 'NHF1_CPD_10J_0h_B',
	'NC1HB14N.1.fastq': 'NHF1_CPD_10J_1h_A',
	'NC1HA2N.1.fastq': 'NHF1_CPD_10J_1h_B',	
	'NC8HB15N.1.fastq': 'NHF1_CPD_10J_8h_A',
	'NC8HA3N.1.fastq': 'NHF1_CPD_10J_8h_B',	
	'NC1DB16N.1.fastq': 'NHF1_CPD_10J_24h_A',
	'NC1DA4N.1.fastq': 'NHF1_CPD_10J_24h_B',	
	'NCPD2DA5N.1.fastq': 'NHF1_CPD_10J_48h_A',
	'NC2DB18N.1.fastq': 'NHF1_CPD_10J_48h_B',
	'NCPD2DB42.1.fastq': 'NHF1_CPD_10J_48h_B2',
	##
	'N6-4-A33N.1.fastq': 'NHF1_6-4_Control_A',
	'N6-4-C35N.1.fastq': 'NHF1_6-4_Control_B',
	'N60HA7N.1.fastq': 'NHF1_6-4_20J_0h_A',	
	'N6-0HB19N.1.fastq': 'NHF1_6-4_20J_0h_B',
	'N620MA8N.1.fastq': 'NHF1_6-4_20J_20min_A',
	'N6-20MB9N.1.fastq': 'NHF1_6-4_20J_20min_B',
	'N62HA11N.1.fastq': 'NHF1_6-4_20J_2h_A',
	'N6-1HB20N.1.fastq': 'NHF1_6-4_20J_2h_B',
	'N61HA10N.1.fastq': 'NHF1_6-4_20J_1h_A',
	'N6-2HB12N.1.fastq': 'NHF1_6-4_20J_1h_B',
	'N64HC34N.1.fastq': 'NHF1_6-4_20J_4h_A',
	'N6-4HB21N.1.fastq': 'NHF1_6-4_20J_4h_B',
	##
	'G-6-4A22.1.fastq': 'GM12787_6-4_20J_nakedDNA_A',
	'G-6-4B24.1.fastq': 'GM12787_6-4_20J_nakedDNA_B',
	'G6-4A23.1.fastq': 'GM12787_6-4_20J_cell_A',
	'G6-4B25.1.fastq': 'GM12787_6-4_20J_cell_B',
	##
	'G-CPDA29.1.fastq': 'GM12787_CPD_20J_nakedDNA_A',
	'G-CPDB31.1.fastq': 'GM12787_CPD_20J_nakedDNA_B',
	'GCPDA36.1.fastq': 'GM12787_CPD_20J_cell_A',
	'GCPDB32.1.fastq': 'GM12787_CPD_20J_cell_B',
	##
	'N6-4HC-29.1.fastq': 'NHF1_6-4_20J_4h_no_CPD_Photolyase',
	'N6-4HCP31.1.fastq': 'NHF1_6-4_20J_4h_with_CPD_Photolyase',
}


# Sample ID	Pool	Dose/damage			Time			Rep (A/B)
# NCPD0hA1N	UVDT1	10J/CPD				0h				A
# NCPD2dA5N	UVDT1	10J/CPD				48h(2days)		A
# NCPD-A6N	UVDT1	0J(No UV)/CPD		-				A
# N6-0hB19N	UVDT1	20J/6-4				0h				B
# N6-20mB9N	UVDT1	20J/6-4				20min			B
# N6-1hB20N	UVDT1	20J/6-4				1h				B
# N6-2hB12N	UVDT1	20J/6-4				2h				B
# N6-4hB21N	UVDT1	20J/6-4				4h				B
# #########
# NC0HB13N	UVDT2	10J/CPD				0h				B
# NC1HB14N	UVDT2	10J/CPD				1h				B
# NC8HB15N	UVDT2	10J/CPD				8h				B
# NC1DB16N	UVDT2	10J/CPD				24h(1day)		B
# NC2DB18N	UVDT2	10J/CPD				48h(2days)		B
# NC-B17N	UVDT2	0J(No UV)/CPD		-				B
# N6-4-C35N	UVDT2	0J(No UV)/6-4		-				C (B)
# NC1HA2N	UVDT2	10J/CPD				1h				A
# NC8HA3N	UVDT3	10J/CPD				8h				A
# NC1DA4N	UVDT3	10J/CPD				24h(1day)		A
# N60HA7N	UVDT3	20J/6-4				0h				A
# N620MA8N	UVDT3	20J/6-4				20min			A
# N61HA10N	UVDT3	20J/6-4				1h				A
# N62HA11N	UVDT3	20J/6-4				2h				A
# N64HC34N	UVDT3	20J/6-4				4h				C (A)
# N6-4-A33N	UVDT3	0J(No UV)/6-4		-				A
# G-6-4A22	HLUVB	20J/6-4				0h (naked DNA)	A
# G6-4A23	HLUVB	20J/6-4				0h (cell)		A
# G-6-4B24	HLUVB	20J/6-4				0h (naked DNA)	B
# G6-4B25	HLUVB	20J/6-4				0h (cell)		B
# G-CPDA29	HLUVB	20J/CPD				0h (naked DNA)	A
# GCPDA36	HLUVB	20J/CPD				0h (cell)		A
# G-CPDB31	HLUVB	20J/CPD				0h (naked DNA)	B
# GCPDB32	HLUVB	20J/CPD				0h (cell)		B

# General utils

class DamageSeqPairedEndPipeline(SeqPipeline):
	'''DamageSeq Paired End Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkPairedEndInput()
		self.treatment = treatments[os.path.basename(self.input)]
		self.genomeSize = '~/ogun/seq/hg19/hg19.bed'
		self.dnaseBed = '/nas02/home/a/d/adebali/ucsc/wgEncodeUwDnaseNHDFAd.fdr01peaks.hg19.bed'

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
			'-p', 8,
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
		slopB = 6
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
		input = self.latestOutput
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
			'-n', nucleotideOrder,
			'--percentage'
		]
		self.nucleotideAbundance = output
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
		self.dinucleotideAbundance = output
		self.run(singleCodeList, runFlag)
		return self

	def plotDinucleotideAbundance(self, runFlag=False):
		input = self.dinucleotideAbundance
		singleCodeList = [
			'plotNucleotideFreqLine.r',
			input,
			self.treatment
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
			#'-Sb'
			'-o',
			output,
			input,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self
#-bf 0X2 bam bedpe

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

	# def bed2sortedBed(self, runFlag=True):
	# 	input = self.latestOutput
	# 	output = self.in2out(input, '.bed', '.sorted.bed')
	# 	log = self.out2log(output)
	# 	singleCodeList = [
	# 		'sortBed',
	# 		'-i', input,
	# 		'>', output
	# 	]
	# 	self.run(singleCodeList, runFlag)
	# 	self.latestOutput = output
	# 	self.sorted=True
	# 	return self

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
			'-k2,2n',
			# '-k3,3n',
			# '-k4,4',
			# '-k5,5',
			# '-k6,6',
			#'-k9,9',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2sortedBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.srt.bed')
		singleCodeList = [
			'sort',
			'-k1,1',
			'-k2,2n',
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
		# self.finalBed = output
		return self

	def bed2bedGraph(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.bdg')
		singleCodeList = [
			'bedtools',
			'genomecov',
			'-i', input,
			'-g', self.genomeSize,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2wig(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.wig')
		singleCodeList = [
			'igvtools',
			'count',
			input,
			output,
			self.genomeSize
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def bed2bigBed(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bed', '.bb')
		singleCodeList = [
			'bedToBigBed',
			input,
			self.genomeSize,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def bedGraph2bigWig(self, runFlag=True):
		input = self.latestOutput
		output = self.in2out(input, '.bdg', '.bw')
		singleCodeList = [
			'wigToBigWig',
			'-clip',
			input,
			self.genomeSize,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def addHeaderToBedGraph(self, runFlag=True):
		input= self.latestOutput
		output = self.in2out(input, '.bdg', '.head.bdg')
		singleCodeList = [
			'echo "track type=bedGraph name=' + self.treatment + '"', '>', output,
			'&&',
			'tail -n +2', input, '>>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def coverageChromHMM(self, runFlag=False):
		input = self.latestOutput
		self.totalMappedReads = self.internalRun(bed.bed(input).getHitNum, [], runFlag)
		outputs = []
		os.system("mkdir -p " + os.path.join(os.path.dirname(input), "chromHMM"))
		for chromHMMfile in chromHMMfiles:
			regionName = chromHMMfile.split('/')[-1].split('_NHFL_')[0]
			os.system("mkdir -p " + os.path.join(os.path.dirname(input), "chromHMM", self.treatment))
			outputDir = os.path.join(os.path.dirname(input), "chromHMM", self.treatment)
			output = os.path.join(outputDir, regionName + '.bed')
			singleCodeList = [
				'bedtools',
				'coverage',
				'-counts',
				'-a', chromHMMfile, # Primary genomic locations
				'-b', input, # Read file, to get the read counts overlapping with the primary genomic locations
				'>', output
			]
			self.run(singleCodeList, runFlag)
			outputs.append(output)
		self.coverageChromHMMoutputs = outputs
		return self

	def normalizeCoverageChromHMM(self, runFlag=False):
		# Normalization between Chromatin States
		multiplicationFactor = float(10000000)/self.totalMappedReads if self.totalMappedReads else 1 # Normalizes counts by total mapped read number
		inputs = self.coverageChromHMMoutputs
		codeDict = {
			"inputs": inputs,
			"outputF": [".bed", ".n1K.bed"],
			"codeList": [
				'bedCount2normalizedCount.py',
				'-i', '#IN',
				'-c', 7,
				'-o', '#OUT',
				'-l', 1000,
				'-m', multiplicationFactor
			]
		}
		self.coverageChromHMMoutputsN1K = self.runM(codeDict, runFlag)
		return self

	def appendNull(self, runFlag=False):
		inputs = self.coverageChromHMMoutputsN1K
		nullValue = "NA"
		maxEntry = 115306
		codeDict = {
			"inputs": inputs,
			"outputF": ['.bed', '.appNA.bed'],
			"codeList": [
				'bedCount2appendedVoid.py',
				'-i', '#IN',
				'-o', '#OUT',
				'-n', maxEntry,
				'-null', 'NaN'
			]
		}
		self.chromHMMappendedByVoidOutputs = self.runM(codeDict, runFlag)
		return self

	def getCountColumn(self, runFlag=False):
		inputs = self.chromHMMappendedByVoidOutputs
		columnNo = 7
		codeDict = {
			"inputs": inputs,
			"outputF": ['.bed', '.txt'],
			"codeList": [
				'cat', '#IN',
				'|',
				'cut', '-f', columnNo,
				'>', '#OUT'
			]
		}
		self.chromHMMcountLists = self.runM(codeDict, runFlag)
		return self

	def mergeChromHMM(self, runFlag=False):
		inputs = self.chromHMMcountLists
		output = os.path.join(os.path.dirname(self.latestOutput), self.treatment + "_mergedChromHMM.txt")
		codeList = ['paste'] + inputs + ['>' + output]
		self.run(codeList, runFlag)
		self.mergedChromHMM = output
		return self

	def plotChromHMM(self, runFlag=False):
		input = self.mergedChromHMM
		codeList = [
			"plotChromHMM.r",
			input,
			self.treatment
		]
		self.run(codeList, runFlag)
		return self

	def dnaseClosest(self, runFlag=False):
		input = self.latestOutput
		output = self.in2out(input,'.bed','.dnaseClosest.bed')
		self.dnaseClosestBed = output
		codeList = [
			'bedtools',
			'closest',
			'-a', self.dnaseBed,
			'-b', input,
			'-k', 1,
			'>', output
		]
		self.run(codeList, runFlag)
		return self

	def retrieveOnlyTheClosest(self, runFlag=False):
		input = self.dnaseClosestBed
		output = self.in2out(input,'.bed','.k1.bed')
		self.dnaseClosestBedK1 = output
		codeList = [
			'sort',
			'-u',
			'-k1,1',
			'-k2,1n',
			'-k3,3n',
			input,
			'>', output
		]
		self.run(codeList, runFlag)
		return self

	def bedClosest2distance(self, runFlag=False):
		input = self.dnaseClosestBedK1
		output = self.in2out(input,'.bed','.bed.distance.txt')
		codeList = [
			'bedClosest2distance.py',
			'-i', input,
			'-c1', 2,
			'-c2', 12,
			'-o', output
		]
		self.run(codeList, runFlag)
		return self

#bsub -M 32 'grep -Pv "\t0" iterativeCov.4.bed | sort -u -k1,1 -k2,2 -k3,3 >escapingUniqueDamages.bed'

input = sys.argv[1]
pipeline = DamageSeqPairedEndPipeline(input)

pipeline\
	.cutadapt(True)\
	.bowtie(True)\
	.sam2bam(True)\
	.bam2bed(True)\
	.bed2uniquelySortedBed(True)\
	.bedpe2bed(True)\
	.slopBed(True)\
	.bed2fixedRangeBed(True)\
	.bed2sortedBed(True)\
		.bed2fasta(True)\
			.fa2nucleotideAbundanceTable(True)\
				.plotNucleotideAbundance(True)\
			.fa2dimerAbundanceTable(True)\
				.plotDinucleotideAbundance(True)\
		.coverageChromHMM(True)\
			.normalizeCoverageChromHMM(True)\
			.appendNull(True)\
			.getCountColumn(True)\
			.mergeChromHMM(True)\
			.plotChromHMM(True)\
		.dnaseClosest(False)\
		.retrieveOnlyTheClosest(False)\
		.bedClosest2distance(False)\
	.bed2bigBed(True)
	# .bed2bedGraph(True)\
	# .addHeaderToBedGraph(True)\
	# .bedGraph2bigWig(True)
	# .bed2wig(True)\
