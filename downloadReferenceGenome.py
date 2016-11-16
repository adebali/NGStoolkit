#!/usr/bin/env python

import os
import sys
from glob import glob
from generalUtils import runGet

assemblyId = sys.argv[1]

root = '/nas02/home/a/d/adebali/ncbi'
runGet('mkdir -p ' + os.path.join(root, assemblyId))
workingDirectory = os.path.join(root, assemblyId)
os.chdir(workingDirectory)

#runGet('wget -O ' + os.path.join(root,'bacteria_assembly_summary.txt') + ' ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt')

ftpAddress = runGet('grep -P "^' + assemblyId + '\\t" ' + os.path.join(root, 'bacteria_assembly_summary.txt') + ' | cut -f20').strip()

runGet('wget ' + ftpAddress + '/*_genomic.fna.gz')

runGet('wget ' + ftpAddress + '/*_genomic.gff.gz')
runGet('gunzip *gz')

fnaFiles = glob('*.fna')

i = 0
for fnaFile in fnaFiles:
	newFnaFile = fnaFile.replace('._genomic.fna', '.fa')
	chrFnaFile = newFnaFile.replace('.fa', '.chr.fa')
	runGet('cp ' + fnaFile + ' ' + newFnaFile)
	runGet('sed \'1 s/.*$/>chr\' ' + newFnaFile + ' >' + chrFnaFile)
	runGet('bowtie-build '  + chrFnaFile + ' ' + assemblyId + '.' + str(i))
	i+= 1
