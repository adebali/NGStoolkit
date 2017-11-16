#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='prints the average length of bed segments')
parser.add_argument('-i', required= True, help='bed file')

args = parser.parse_args()

bedFile = args.i

print(bed.bed(bedFile).getAverageLength())
