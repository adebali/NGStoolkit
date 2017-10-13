#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random
import sequence
import re

parser = argparse.ArgumentParser(description='count strings')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-fasta', required= True, help='fasta reference')
parser.add_argument('-dna', required= True, help='dna sequence to count')

args = parser.parse_args()

def substring_indexes(substring, string, onebased=True):
    """ 
    Generate indices of where substring begins in string

    >>> list(find_substring('me', "The cat says meow, meow"))
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

bedFile = args.i
out = args.o
fasta = args.fasta
dna = sequence.DNA(args.dna)
dna_string = dna.getSequence()
dna_rc = dna.reverseComplement()
dna_rc_string = dna_rc.getSequence()

bedObject = bed.bed(bedFile)

for bedline in bedObject.read():
    intervalSequence = bedline.getSequence(fasta)
    positions = substring_indexes(dna_string, intervalSequence)
    positions_opp = substring_indexes(dna_rc_string, intervalSequence)
    for position in positions:
        fields = bedline.fields()
        fields.append(str(position))
        fields.append('same')
        line = bedline.fields2line(fields)
        out.write(line + '\n')
    for position in positions_opp:
        fields = bedline.fields()
        fields.append(str(position))
        fields.append('opposite')
        line = bedline.fields2line(fields)
        out.write(line + '\n')
    