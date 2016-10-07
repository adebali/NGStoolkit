#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='fix fragment length if possible')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-s', required= True, help='side')
parser.add_argument('-l', required=  True, help='sequence length of interest')

args = parser.parse_args()

bedFile = args.i
side = args.s
length = int(args.l)

bed.bed(bedFile).fixRange(side, length)
