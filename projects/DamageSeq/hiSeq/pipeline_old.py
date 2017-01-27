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

class DamageSeqPairedEndPipeline(SeqPipeline):
	'''DamageSeq Paired End Pipeline class'''

	def __init__(self, input, reference='nucleosome', runDir=False):
		SeqPipeline.__init__(self, input)
		fileBaseName = os.path.basename(self.input)
		self.checkPairedEndInput()
		samples = generalUtils.table2dictionary(generalUtils.file('dataDir/samples_minReadCount.csv'), 'sample')
		sampleDictionary = samples[fileBaseName][0]
		self.treatment = sampleDictionary['treatment_title']
		self.minimumReadCount = round(int(sampleDictionary['minReadCount']) - 500000)
		self.group = sampleDictionary['group']
		self.sampleMinPyrCount = sampleDictionary['minPyrHitNo']
		if '-group' in sys.argv:
			if self.group != int(sys.argv[sys.arg.index('-group') + 1]):
				self.runMode = False
		self.dnaseBed = '/nas02/home/a/d/adebali/ucsc/wgEncodeUwDnaseNHDFAd.fdr01peaks.hg19.bed'
		self.fName = generalUtils.getFunctionName

		
		if reference == 'nucleosome':
			self.extra_name = 'ncl'
			self.reference = reference
			self.bowtie_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19'
			self.fasta_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19.fa'
			self.chromosome_limits = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19.chrSizes.bed'
		else:
			raise ValueError('a wrong reference is stated')


		if runDir != False:
			self.outputDirectory = os.path.join(os.path.dirname(self.input), runDir)
			os.system('mkdir -p ' + self.outputDirectory)


	# Methods
	##########################################################

	def cutadapt_fastq2fastq(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output1basename = os.path.basename(self.in2out(self.input1, '1.fastq', self.extra_name + '.cut.1.fastq'))
		output2basename = os.path.basename(self.in2out(self.input2, '2.fastq', self.extra_name + 'cut.2.fastq'))
		output1 = os.path.join(self.outputDirectory, output1basename)
		output2 = os.path.join(self.outputDirectory, output2basename)
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
		self.latestOutput = output1
		self.latestOutput1 = output1
		self.latestOutput2 = output2
		return self

	def bowtie_fastq2sam(self, runFlag=True):
		if runFlag == "Skip": return self
		input1 = self.latestOutput1
		input2 = self.latestOutput2
		output = self.in2out(self.latestOutput, '.fastq', '.bow.sam')
		log = self.out2log(output)
		singleCodeList = [
			'bowtie',
			'-t', self.bowtie_reference,
			'-q', # FASTAQ input (default)
			'--nomaqround', # Do NOT round MAC
			'--phred33-quals', # Depends on the sequencing platform
			'-S', # Output in SAM format
			'-m', 4, # Do not report the reads that are mapped on to more than 4 genomic locations
			'-X', 1000,
			'--seed', 123, # Randomization parameter in bowtie,
			'-p', 4,
			'-1', input1,
			'-2', input2,
			output,
			'>', log
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def slopBed_bed2bed(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		slopB = 6
		# output = self.in2out(input, 'bed', 'slopB' + str(slopB) + '.bed')
		output = self.funIn2out(generalUtils.getFunctionName(), input, 'b6')

		genome = self.chromosome_limits
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

	def convertBedToFasta_bed2fa(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.latestOutput
		# output = self.in2out(input, '.bed', '.bed.fa')
		output = self.funIn2out(generalUtils.getFunctionName(), input)

		referenceFasta = self.fasta_reference
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
		self.latestOutput = output
		return self

	def getPyrimidineDimers_fa2bed(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)
		# output = self.in2out(input, '.fa', '.Pyr.bed')
		singleCodeList = [
			'fa2bedByChoosingReadMotifs.py',
			'-i', input,
			'-o', output,
			'-r', '\'.{4}[T|C][T|C].{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
		]
		self.run(singleCodeList, runFlag)
		self.finalBed = output

		self.latestBed = output
		self.latestOutput = output
		return self

	def sampleFromBed_bed2bed(self, runFlag=False):
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)
		singleCodeList = [
			'bed2sampling.py',
			'-i', input,
			'-o', output,
			'-c', self.sampleMinPyrCount
		]
		self.run(singleCodeList, runFlag)
		self.latestBed = output
		self.latestOutput = output
		return self

	def getNucleotideAbundanceTable_fa2csv(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.fasta
		output = self.funIn2out(generalUtils.getFunctionName(), input)
#		output = self.in2out(input, '.fa', '.fa.nucAbu.csv')
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

	def plotNucleotideAbundance_csv2pdf(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.nucleotideAbundance
		nucleotideOrder = 'TCGA'
		singleCodeList = [
			'plotNucleotideAbundance.r',
			input,
			self.treatment
		]
		self.run(singleCodeList, runFlag)
		return self

	def getDimerAbundanceTable_fa2csv(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.fasta
		output = self.funIn2out(generalUtils.getFunctionName(), input)
		#output = self.in2out(input, '.fa', '.fa.dimerAbu.csv')
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

	def plotDinucleotideAbundance_csv2pdf(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.dinucleotideAbundance
		singleCodeList = [
			'plotNucleotideFreqLine.r',
			input,
			self.treatment
		]
		self.run(singleCodeList, runFlag)
		return self



	def convertToBam_sam2bam(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
#		output = self.in2out(input, '.sam', '.bam')
		output = self.funIn2out(generalUtils.getFunctionName(), input)
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

	def convertToBed_bam2bedpe(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
#		output = self.in2out(input, '.bam', '.bed')
		output = self.funIn2out(generalUtils.getFunctionName(), input)
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
	#	if runFlag == "Skip": return self
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

	def removeRedundancy_bedpe2bedpe(self, runFlag=True):
		if runFlag == "Skip": return self
		self.checkDependency_(['sorted'])
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)
#		output = self.in2out(input, '.bed', '.NR.bed')
		singleCodeList = [
			'bed2nonRedundantBed.py',
			'-i', input,
			'-o', output,
			'-c', '1,2,3,4,5,6,7,9' # columns
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def uniqueSort_bedpe2bedpe(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bed', '.uSorted.bed')
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

	def sortBed_bed2bed(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bed', '.srt.bed')
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


	def convertBedpeToSingleFrame_bedpe2bed(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		#output = self.in2out(input, '.bed', '.sf.bed')
		output = self.funIn2out(generalUtils.getFunctionName(), input)
		singleCodeList = [
			'bedpe2bed.py',
			input,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def convertToFixedRange_bed2bed(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		side = 'left'
		sideInput = side[0]
		fixedRange = 10
#		output = self.in2out(input, '.bed', '.fr' + sideInput.upper() + str(fixedRange) + '.bed')
		output = self.funIn2out(generalUtils.getFunctionName(), input, str(fixedRange))
		
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

	def convertToBedGraph_bed2bdg(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bed', '.bdg')
		singleCodeList = [
			'bedtools',
			'genomecov',
			'-i', input,
			'-g', self.chromosome_limits,
			'>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def convertToWig_bed2wig(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
#		output = self.in2out(input, '.bed', '.wig')
		output = self.funIn2out(generalUtils.getFunctionName(), input)
		
		singleCodeList = [
			'igvtools',
			'count',
			input,
			output,
			self.chromosome_limits
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def convertToBigBed_bed2bb(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bed', '.bb')
		singleCodeList = [
			'bedToBigBed',
			input,
			self.chromosome_limits,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self


	def convertToBigWig_bdg2bw(self, runFlag=True):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bdg', '.bw')
		singleCodeList = [
			'wigToBigWig',
			'-clip',
			input,
			self.chromosome_limits,
			output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def addHeader_bdg2bdg(self, runFlag=True):
		if runFlag == "Skip": return self
		input= self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input, '.bdg', '.head.bdg')
		singleCodeList = [
			'echo "track type=bedGraph name=' + self.treatment + '"', '>', output,
			'&&',
			'tail -n +2', input, '>>', output
		]
		self.run(singleCodeList, runFlag)
		self.latestOutput = output
		return self

	def coverageChromHMM(self, runFlag=False):
		if runFlag == "Skip": return self
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

	def normalizeCoverageChromHMM_bed2bed(self, runFlag=False):
		if runFlag == "Skip": return self
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

	def appendNull_bed2bed(self, runFlag=False):
		if runFlag == "Skip": return self
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

	def getCountColumn_bed2txt(self, runFlag=False):
		if runFlag == "Skip": return self
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
		if runFlag == "Skip": return self
		inputs = self.chromHMMcountLists
		output = os.path.join(os.path.dirname(self.latestOutput), self.treatment + "_mergedChromHMM.txt")
		codeList = ['paste'] + inputs + ['>' + output]
		self.run(codeList, runFlag)
		self.mergedChromHMM = output
		return self

	def plotChromHMM(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.mergedChromHMM
		codeList = [
			"plotChromHMM.r",
			input,
			self.treatment
		]
		self.run(codeList, runFlag)
		return self

	def dnaseClosest_bed2bed(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.latestOutput
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input,'.bed','.dnaseClosest.bed')
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

	def retrieveOnlyTheClosest_bed2bed(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.dnaseClosestBed
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input,'.bed','.k1.bed')
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

	def getDistanceFromClosest_bed2txt(self, runFlag=False):
		if runFlag == "Skip": return self
		input = self.dnaseClosestBedK1
		output = self.funIn2out(generalUtils.getFunctionName(), input)		
		#output = self.in2out(input,'.bed','.bed.distance.txt')
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
pipeline = DamageSeqPairedEndPipeline(input, 'nucleosome', '1229')

(pipeline
	.cutadapt_fastq2fastq(False)	
	.bowtie_fastq2sam(False)
	.convertToBam_sam2bam(False)
	.convertToBed_bam2bedpe(False)
	.uniqueSort_bedpe2bedpe(False)
	.convertBedpeToSingleFrame_bedpe2bed(False)
	.slopBed_bed2bed(False)
	.convertToFixedRange_bed2bed(False)
	.sortBed_bed2bed(False)
	.convertBedToFasta_bed2fa(False)
		
		.getNucleotideAbundanceTable_fa2csv(False)
		.plotNucleotideAbundance_csv2pdf(False)
		
		.getDimerAbundanceTable_fa2csv(False)
		.plotDinucleotideAbundance_csv2pdf(False)
	
	.getPyrimidineDimers_fa2bed(False)
	.sampleFromBed_bed2bed(True)
		
		.coverageChromHMM(False)
			.normalizeCoverageChromHMM_bed2bed(False)
			.appendNull_bed2bed(False)
			.getCountColumn_bed2txt(False)
			.mergeChromHMM(False)
			.plotChromHMM(False)
		
		.dnaseClosest_bed2bed(False)
			.retrieveOnlyTheClosest_bed2bed(False)
			.getDistanceFromClosest_bed2txt(False)
	
	.convertToBigBed_bed2bb(False)
)
	# .bed2bedGraph(True)\
	# .addHeaderToBedGraph(True)\
	# .bedGraph2bigWig(True)
	# .bed2wig(True)\
