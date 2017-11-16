#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random
import sequence
import re

parser = argparse.ArgumentParser(description='produces motif bed file given input bed file')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-fasta', required= True, help='fasta reference')
parser.add_argument('-motif', required= True, help='dna motif')
args = parser.parse_args()
bedFile = args.i
out = args.o
fasta = args.fasta
dna = sequence.DNA(args.motif)
dna_string = dna.getSequence() 
dna_string_rc = dna.reverseComplement(True).getSequence()
motifLength = dna.getLength()

def substring_indexes(substring, string, onebased=True):
    """ 
    Generate indices of where substring begins in string

    >>> list(substring_indexes('me', "The cat says meow, meow"))
    [13, 19]
    """
    last_found = -1  # Begin at -1 so the next position to search from is 0
    while True:
        # Find next index of substring, by starting after its last known position
        last_found = string.find(substring, last_found + 1)
        if last_found == -1:  
            break  # All occurrences have been found
        if onebased:
            yield last_found + 1
        else:
            yield last_found

def reversePositionCounts(positions, sequenceLength, motifLength):
    newPositions = []
    for position in positions:
        newPositions.append(sequenceLength - position + 1 - motifLength)
    return newPositions

bedObject = bed.bed(lines=bedFile.readlines())
for bedline in bedObject.read():
    sequence = bedline.getSequence(args.fasta)
    positions = list(substring_indexes(dna_string, sequence))
    positions_opp = reversePositionCounts(list(substring_indexes(dna_string_rc, sequence)), len(sequence), motifLength)
    strandList = ['+'] * len(positions) + ['-'] * len(positions_opp)
    allPositions = positions + positions_opp
    for i in range(len(allPositions)):
        position = allPositions[i]
        # fields = bedline.fields() + [
        fields = [
            bedline.chromosome(),
            str(bedline.start() + position),
            str(bedline.start() + position + motifLength - 1),
            bedline.name() + '_' + dna_string,
            '.',
            strandList[i]
        ]
        out.write('\t'.join(fields) + '\n')


