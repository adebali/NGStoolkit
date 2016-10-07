import os
import sys
import subprocess
import generalUtils
from seqPipeline import SeqPipeline
from glob import glob

referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'

class rnaSeq(SeqPipeline):
	'''rna-seq Paired End Pipeline class'''

	def __init__(self, input):
		SeqPipeline.__init__(self, input)
		self.checkInput()
		self.referenceBowtieIndex = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.0'
		self.referenceFasta = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.fa')
		self.e_coli_K12_windowsBedFile = generalUtils.file('/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed')

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
		e_coli_K12_windowsBedFile = '/nas02/home/a/d/adebali/ncbi/GCF_000005845.2/GCF_000005845.2.chr.w200.bed'

		for i in range(len(inputs)):
			singleCodeList = [
					'bedtools',
					'coverage',
					'-counts',
					'-a', e_coli_K12_windowsBedFile,
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
					'-t', 4,
					'-m', 100
				]
			self.run(singleCodeList, runFlag)
		self.latestOutput = outputs[0]
		self.latestOutputPlus = outputPlus
		self.latestOutputMinus = outputMinus
		return self



input = sys.argv[1]

pipeline = rnaSeq(input)

pipeline\
	.trimNonqualifiedNucleotides(False)\
	.removePolyAtail(False)\
	.align(True)\
	.sam2bam(True)\
	.sortBam(True)\
	.bam2bed(True)\
	.splitBedByStrand(True)\
	.windowCoverage(True)\
	.bedCount2percentageByMax(True)
