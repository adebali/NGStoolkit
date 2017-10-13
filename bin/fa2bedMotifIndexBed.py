#!/usr/bin/env python
import fasta
import argparse
import re
import sys
from sequence import reMotif, consensus
import bed

parser = argparse.ArgumentParser(description='converts fasta to bed formatted intervals of an index of motif of interest')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-r', required= False, help='regex motif that is expected')
parser.add_argument('-c', required= False, help='consensus motif that is expected')
parser.add_argument('--DNAreverseComplement', required= False, action= "store_true",
                    help='DNA reverse complement flag of motif.\
                    Finds the reverse complement motif as well.\
                    Use with caution: Your motif must be symmetric')

args = parser.parse_args()
fastaObject = fasta.fasta(args.i)
out = args.o
if args.r:
    regexMotif = reMotif(args.r)
elif args.c:
    regexMotif = consensus(args.c)
else:
    raise ValueError('-r or -c must be specified')
DNAreverseComplementFlag = args.DNAreverseComplement


def fastaHeader2bedLine(header):
    chr = header.replace('>','').split(':')[0]
    start = int(header.split(':')[1].split('(')[0].split('-')[0])
    end = int(header.split(':')[1].split('(')[0].split('-')[1])
    return chr + '\t' + str(start) + '\t' + str(end)

def seqObject2bedLine(seq, regexMotif):
    lines = []
    indexList = regexMotif.getIndexList(seq.getSequence())
    startList = indexList[0]
    endList = indexList[1]
    for i in range(len(startList)):
        bedLine = bed.bedline(fastaHeader2bedLine(seq.getHeader()))
        lines.append(bedLine.chromosome() + "\t" + str(bedLine.start() + startList[i]) + "\t" + str(bedLine.start() + endList[i]) + "\t" + "+" + "\n")
    if DNAreverseComplementFlag:
        reverseComplementMotif = regexMotif.DNA_complement().reverse()
        indexList = reverseComplementMotif.getIndexList(seq.getSequence())
        startList = indexList[0]
        endList = indexList[1]
        for i in range(len(startList)):
            bedLine = bed.bedline(fastaHeader2bedLine(seq.getHeader()))
            lines.append(bedLine.chromosome() + "\t" + str(bedLine.start() + startList[i]) + "\t" + str(bedLine.start() + endList[i]) + "\t" + "-" + "\n")
    return lines

seqObjects = fastaObject.stream2(100*4096)
for seq in seqObjects:
    for line in seqObject2bedLine(seq, regexMotif):
        out.write(line)
out.close()