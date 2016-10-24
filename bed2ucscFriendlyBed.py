#!/usr/bin/env python
import os
import sys

def bedLine2UCSCfriendlyLine(line):
	ll = line.strip().split('\t')
	chrom = ll[0]
	chromStart = ll[1]
	chromEnd = ll[2]
	name = ll[3]
	score = ll[4]
	strand = ll[5]
	thickStart = chromStart
	thickEnd = chromEnd
	if strand == '+':
		color = '255,0,0'
	elif strand == '-':
		color = '255,90,0'
	blockCount = '1'
	blockSizes = '1'
	blockStarts = '0'

	# newLine = [chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, color, blockCount, blockSizes, blockStarts]
	newLine = [chrom, chromStart, chromEnd, name, score, strand]
	return 	'\t'.join(newLine)

print('track name="EscapingDamages" description="Consisted Damages" colorByStrand="255,0,0 0,0,255"')

bedFile = sys.argv[1]
filein = open(bedFile, 'r')
for line in filein:
	print(bedLine2UCSCfriendlyLine(line))
