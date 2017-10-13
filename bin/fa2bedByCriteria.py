#!/usr/bin/env python

import fasta
import argparse

parser = argparse.ArgumentParser(description='converts fasta to bed by applying a criterion')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-c', required= True, help='criterion (seq is sequence variable)')

args = parser.parse_args()

fastaObject = fasta.fasta(args.i)
output = args.o
criterion = args.c

def criterionPass(citerion):
    code = 'True if ' + citerion + ' else False' 
    result = eval(code)
    return result

def fastaHeader2bedLine(header):
    chr = header.replace('>','').split(':')[0]
    start = int(header.split(':')[1].split('(')[0].split('-')[0])
    end = int(header.split(':')[1].split('(')[0].split('-')[1])
    strand = header.split(':')[1].split('(')[1][0]
    return chr + '\t' + str(start) + '\t' + str(end) + '\t' + strand

out = open(output, 'w')
seqDicts = fastaObject.stream()
for seqDict in seqDicts:
    header = seqDict['h']
    sequence = seqDict['s']
    if criterionPass(criterion):
        out.write(fastaHeader2bedLine(header) + '\n')
out.close()