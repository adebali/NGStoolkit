#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='get subsequences having the motif of interest')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-m', required= True, help='motif')
parser.add_argument('-r', required= True, help='right of motif (Integer)')
parser.add_argument('-l', required=  False, help='left of motif(Integer)')

args = parser.parse_args()
motif = args.m
rightNumber = int(args.r)
leftNumber = int(args.l)

lengthOfMotif = len(motif)
fastaDict = fasta.fasta(args.i).read()

for header in fastaDict.keys():
	sequenceCount = 0
	sequence = fastaDict[header]
	for i in range(len(sequence)):
		subseq = sequence[i : i + lengthOfMotif]
		if subseq == motif:
			startPosition = i - leftNumber
			endPosition = i + lengthOfMotif + rightNumber
			currentSequence = sequence[startPosition : endPosition]
			if len(currentSequence) == lengthOfMotif + leftNumber + rightNumber:
				sequenceCount += 1
				print('>' + header + '.' + str(sequenceCount) + '\n' + currentSequence)
