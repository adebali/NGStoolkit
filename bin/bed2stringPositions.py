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
parser.add_argument('--flanking', action='store_true', default=False, help='''
    if bed file is for the flanking regions of interest such as TSS upstream and TES downstream.
    Input must go line by line: upstream, downstream, up, down ...
    ''')
parser.add_argument('--fixed', action='store_true', default=False, help='should be used when for flankings we want to get absolute range')
parser.add_argument('-w', type=int, default=1, help='bin (window) size')

args = parser.parse_args()

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

bedFile = args.i
out = args.o
fasta = args.fasta
dna = sequence.DNA(args.dna)
dna_string = dna.getSequence()
dna_complementary = dna.reverseComplement(False)
dna_complementary_string = dna_complementary.getSequence()

bedObject = bed.bed(bedFile)

for bedline in bedObject.read():

    name = bedline.name()
    fields = bedline.fields()
    fields += [bedline.chromosome() + bedline.start() + bedline.end() + bedline.name() + bedline.score() + bedline.strand()]
    intersectline = bed.intersectLine(bedline.fields2line(fields))



    if args.flanking:
        if '_upstream' in name:
            region = 'upstream'
        elif '_downstream' in name:
            region = 'downstream'
        else:
            raise ValueError("Expected keyword (up/downstream) not found in name: " + name)
    else:
        region = 'main'


    if region == 'main':
        distancePercentage = intersectLine.getDistancePercentage()
        distance = distancePercentage
    if region != 'main':
        if args.fixed:
            distance = intersectLine.getAbsoluteDistance(region, args.w)
        else:
            distancePercentage = intersectLine.getDistancePercentage()   
            if region == 'upstream':
                distancePercentage -= 100
            elif region == 'downstream':
                distancePercentage += 100
            else:
                raise ValueError('Unexpected region ' + region)
            distance = distancePercentage
    sameStrands = intersectLine.sameStrands()


    intervalSequence = bedline.getSequence(fasta)
    positions = list(substring_indexes(dna_string, intervalSequence))
    positions_opp = list(substring_indexes(dna_complementary_string, intervalSequence))
    strandList = ['same'] * len(positions) + ['opposite'] * len(positions_opp)

    allPositions = positions + positions_opp
    for i in range(len(allPositions)):
        position = allPositions[i]
        newFields = list(fields)
        newFields[7] = newFields[7] + position
        newFields[8] = newFields[7] + position + len(dna_string) - 1
        newIntersectLine = bed.intersectLine('\t'.join(newFields))

        strand = strandList[i]
        length = bedline.length()


        if region == 'main':
            distancePercentage = getDistancePercentage(position, length)
            distance = distancePercentage
        if region != 'main':
            if args.fixed:
                distance = getAbsoluteDistance(position, length)
            else:
                distancePercentage = intersectLine.getDistancePercentage()   
                if region == 'upstream':
                    distancePercentage -= 100
                elif region == 'downstream':
                    distancePercentage += 100
                else:
                    raise ValueError('Unexpected region ' + region)
                distance = distancePercentage
        
        
        fields = bedline.fields()

        fields.append(str(position))
        fields.append(strand)
        line = bedline.fields2line(fields)
        out.write(line + '\n')
    for position in positions_opp:
        fields = bedline.fields()
        fields.append(str(position))
        fields.append('opposite')
        line = bedline.fields2line(fields)
        out.write(line + '\n')
    

    