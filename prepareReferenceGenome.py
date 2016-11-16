#!/usr/bin/env python

import os
import sys
from glob import glob
from generalUtils import runGet
import fasta

dir = sys.argv[1]


class referenceGenomePreparation():

	def __init__(self, input):
		self.directory = dir
		directoryGiven = os.path.isdir(self.directory)
		if not directoryGiven:
			raise ValueError('given input is not a directory')
		self.basename = os.path.basename(self.directory)
		self.bowtieIndexName = os.path.join(self.directory, self.basename)
		self.fasta = os.path.join(self.directory, self.basename + '.fa')
		self.gff3 = os.path.join(self.directory, self.basename + '.gff3')
		self.chrFasta = os.path.join(self.directory, self.basename + '.chr.fa')
		self.chromSizes = os.path.join(self.directory, self.basename + '.chrom.sizes')
		self.genesBed = os.path.join(self.directory, self.basename + '.chr.genes.bed')
		self.geneList = os.path.join(self.directory, self.basename + '.genes.txt')

	def fa2chrFa(self):
		runGet('sed \'1 s/.*$/>chr/g\' ' + self.fasta + ' >' + self.chrFasta)
		return self

	def fa2genomeSize(self):
		Fasta = fasta.fasta(self.fasta)
		with open(self.chromSizes, 'w') as out:
			out.write(Fasta.singleEntry2chromSize())
		return self

	def gff3toGenesBed(self):
		runGet('grep -P "\\tgene\\t" ' + self.gff3 + ' | cut -f 1,4,5,7 >' + self.genesBed)
		runGet('sed -i s/' + self.basename + '/chr/g ' + self.genesBed)
		runGet('grep -P "\\tgene\\t" ' + self.gff3 + ' | cut -f 9 | cut -d\';\' -f 3 | cut -d\'=\' -f 2 >' + self.geneList)
		return self

	def bowtieBuild(self):
		runGet('bowtie-build '  + self.chrFasta + ' ' + self.bowtieIndexName)
		return self

pipeline = referenceGenomePreparation(dir)
pipeline\
	.fa2chrFa()\
	.fa2genomeSize()\
	.gff3toGenesBed()\
	.bowtieBuild()

