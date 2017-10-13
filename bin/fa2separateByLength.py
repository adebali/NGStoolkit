#!/usr/bin/env python
# Not tested!
import fasta
import argparse

parser = argparse.ArgumentParser(description='get subsequences having the motif of interest')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-l', required=  False, help='comma separated lengths of interest')

args = parser.parse_args()
input = args.i
seqLengths = args.l.split(',')
for seqLength in seqLengths:
	fileDict[seqLength] = open(input + '.' + str(seqLength) + '.fa', 'w')
fastaDict = fasta.fasta(input).read()

for header in fastaDict.keys():
	sequence = fastaDict[header]
	strSequenceLength = str(len(sequence))
	if strSequenceLength in seqLengths:
		fileDict[strSequenceLength].write('>' + header + '\n' + sequence + '\n')
