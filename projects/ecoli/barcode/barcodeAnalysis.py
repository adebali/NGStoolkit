#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='XR-seq Barcode Analysis')
parser.add_argument('-i', required= True, help='input: fastq (line) file that includes the full length of sequence in the header as the last item')
parser.add_argument('-o', required= True, help='output: results')
args = parser.parse_args()
input = args.i
out = open(args.o, 'w')


def getBarcode(sequence, fullSequence, adapterSequence = 'GGCTCAGTTCGTATGAGTGCCG', barcodeLength = 6):
    sequenceLength = len(sequence)
    ambigiousAdapterLength = len(adapterSequence)
    return fullSequence[sequenceLength + ambigiousAdapterLength : sequenceLength + ambigiousAdapterLength + 6]

# testSequence = 'NNNNNNNNNNNNN' + 'XXXXXXXXXXXXXXXXXXXXXX' + 'AAAAAA' + 'NNNNNNNNNNNNNNNN'
# print(getAmbigiousAdapterSequence('NNNNNNNNNNNNN', testSequence))

def analyzeList(barcodeListPerRead, fastqSeqList):
    #
    totalIdenticalReadCount = len(barcodeListPerRead)
    barcodeSet = set(barcodeListPerRead)
    uniqueBarcodeCount = len(barcodeSet)
    #

    # if uniqueBarcodeCount > 2000:
    #     print(barcodeSet)
    #     print(uniqueBarcodeCount)
    #     sys.exit()

    barcodeCountPerSequence = {}
    for barcode in barcodeSet:
        barcodeCountPerSequence[barcode] = barcodeListPerRead.count(barcode)
    
    topTenAbundantBarcodeCount = [0] * 10
    i = 0
    for key in sorted(barcodeCountPerSequence, key=barcodeCountPerSequence.get, reverse=True):
        topTenAbundantBarcodeCount[i] = barcodeCountPerSequence[key] 
        if i == 9:
            break
        i += 1
    return '\t'.join(fastqSeqList + [str(uniqueBarcodeCount), str(totalIdenticalReadCount), str(topTenAbundantBarcodeCount), str(topTenAbundantBarcodeCount[0]/float(totalIdenticalReadCount)), str(topTenAbundantBarcodeCount[0])]) + '\n'


filein = open(input, 'r')
prevSeq = ''
barcodeListPerRead = []
for line in filein:
    fastqSeqList = line.strip().split('\t')
    sequence = fastqSeqList[1]
    sequenceLength = len(sequence)
    uncutSequence = fastqSeqList[0].split(' ')[-1]
    barcode = getBarcode(sequence, uncutSequence)
    if sequence != prevSeq:
        if prevSeq != '':
            out.write(analyzeList(barcodeListPerRead, fastqSeqList))
        barcodeListPerRead = []
    barcodeListPerRead.append(barcode)
    prevSeq = sequence
else:
    if prevSeq == sequence:
        barcodeListPerRead.append(barcode)
    out.write(analyzeList(barcodeListPerRead, fastqSeqList))
out.close()
