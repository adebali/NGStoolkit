#!/usr/bin/env python
import os
import bed
import sys
import argparse
import generalUtils
import fasta
from sequence import DNA

parser = argparse.ArgumentParser(description='normalize bed counts by the number of a motif found in the fragment')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-f', required= True, help='<Required> fasta')
parser.add_argument('-m', required= True, help='<Required> motif')
parser.add_argument('-c', required= True, help='<Required> countTab')
parser.add_argument('-s', required= True, help='<Required> strandInfo, +, - or column No')
parser.add_argument('-o', required=  True, help='<Required> output')

args = parser.parse_args()
bedFile =  args.i
outBedFile = args.o
fastaFile = args.f
motif = args.m
printFlag = args.conPrint

fastaInput = fasta.fasta(fastaFile)
headers = fastaInput.read()


def conPrint(text):
    if printFlag:
        print(text)

if args.perNmotif:
    perNmotif = int(args.perNmotif) 
else:
    perNmotif = (1/4)^(len(motif))*1000

def motifCount(sequence, positionStart, positionEnd, strand, motif):
    Sequence = DNA(sequence)
    Subseq = Sequence.subseq(positionStart, positionEnd)
    if strand == '-':
        Subseq = Subseq.reverseComplement()
    motifNum = Subseq.motifCount(motif)
    return motifNum

def line2strand(line):
    if args.s == '+' or args.s == '-':
        return args.s
    else:
        strandColIndex = int(args.s) - 1
        strand = line.split('\t')[strandColIndex]
        if not (strand == '+' or strand == '-'):
            raise ValueError(strand + ' --> Strand information is not correct.')
        return strand

countTabNumber = int(args.c)

filein = open(bedFile, 'r')
out = open(outBedFile, 'w')

for line in filein:
    ll = line.strip().split('\t')
    start = ll[1]
    end = ll[2]
    strand = line2strand(line)

    sequence = headers[ll[0]]
    motifNum = motifCount(sequence, int(start), int(end), strand, motif)
    conPrint(motifNum)

    count = float(ll[countTabNumber - 1])
    if motifNum != 0:
        newCount = (float(count) / motifNum) * float(perNmotif)
    else:
        newCount = 0

    newLine = ''
    for i in range(len(ll)):
        if i == countTabNumber - 1:
            field = str(newCount)
        else:
            field = ll[i]
        newLine += field + '\t'
    newLine = newLine.strip() + '\n'
    out.write(newLine)

out.close()
