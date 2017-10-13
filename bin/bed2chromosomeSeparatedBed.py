#!/usr/bin/env python
import os
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='separate bed by chromosomes')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-chr', nargs="+", required= True, help='chromosomes - space separated')
parser.add_argument('-o', nargs="+", required= False, help='outputs')

args = parser.parse_args()

bedFile = args.i
chromosomes = args.chr
if args.o:
    outputFiles = args.o
    if len(chromosomes) != len(outputFiles):
        raise ValueError("chromosome and output file list should be at equal lengths")
# def input2outputName(input, chr):
#     if input.endswith(".bed"):
#         output = input[:-4] + "_" + chr + ".bed"
#     else:
#         output = input + "_" + chr + ".bed"
#     return output

# chromosomeSet = set()
# for bedline in bed.bed(bedFile).read():
#     chromosomeSet.add(bedline.chromosome())

# print(chromosomeSet)
i = 0
for chromosome in chromosomes:
    if args.o:
        output = outputFiles[i]
    else:
        output = input2outputName(bedFile, chromosome)
    os.system('grep -P "^' + chromosome + '" ' + bedFile + ' >' + output)
    i += 1