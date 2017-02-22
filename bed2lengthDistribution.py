#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='prints length distribution of bed windows')
parser.add_argument('-i', required= True, help='bed file')

args = parser.parse_args()

bedFile = args.i

bed.bed(bedFile).lengthDistribution()
# bed.bed(bedFile).printLengths()
