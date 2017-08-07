#!/usr/bin/env python
import os
import sys
import argparse
import math

parser = argparse.ArgumentParser(description='XR-seq Barcode Analysis')
parser.add_argument('-i', required= True, help='input: fastq (line) file that includes the full length of sequence in the header as the last item')
parser.add_argument('-o', required= True, help='output: results')
args = parser.parse_args()
input = args.i
out = open(args.o, 'w')

def totalCount2expectedMaxSingleton(x):
    # print(x)
    # print(math.ceil(50.38*(x**2) - 142.01*x))
    return math.ceil(50.38*(x**2) - 142.01*x)
    # return 1

def artifactFlag(x, totalCount):
    expectedMinCountToGiveTheMaxUniqueCount = totalCount2expectedMaxSingleton(x)
    if expectedMinCountToGiveTheMaxUniqueCount > totalCount:
        return True
    else:
        return False

def getBarcode(sequence, fullSequence, adapterSequence = 'GGCTCAGTTCGTATGAGTGCCG', barcodeLength = 6):
    sequenceLength = len(sequence)
    ambigiousAdapterLength = len(adapterSequence)
    return fullSequence[sequenceLength + ambigiousAdapterLength : sequenceLength + ambigiousAdapterLength + 6]

def reduceArtifacts(barcodeListPerRead, fastqSeqList, out):
    #
    artifacts = []
    totalIdenticalReadCount = len(barcodeListPerRead)
    barcodeSet = set(barcodeListPerRead)
    uniqueBarcodeCount = len(barcodeSet)
    #

    barcodeCountPerSequence = {}
    for barcode in barcodeSet:
        barcodeCountPerSequence[barcode] = barcodeListPerRead.count(barcode)
    
    for barcode in sorted(barcodeCountPerSequence, key=barcodeCountPerSequence.get, reverse=True):
        identicalBarcodeCount = barcodeCountPerSequence[barcode]
        # print(barcodeCountPerSequence)
        if identicalBarcodeCount == 1:
            break
        if artifactFlag(identicalBarcodeCount, totalIdenticalReadCount):
            artifacts.append(barcode)
        else:
            break
    # print("ARTIFACTS:")
    # print(artifacts)

    newFastqSeqList = []
    savedArtifactBarcodes = set()
    for seqLine in fastqSeqList:
        fastqSeqList = seqLine.strip().split('\t')
        sequence = fastqSeqList[1]
        uncutSequence = fastqSeqList[0].split(' ')[-1]
        barcode = getBarcode(sequence, uncutSequence)
        if (barcode not in artifacts) or (barcode not in savedArtifactBarcodes):
            newFastqSeqList.append(seqLine)
            savedArtifactBarcodes.add(barcode)

    for seqLine in newFastqSeqList:
        out.write(seqLine)

filein = open(input, 'r')
prevSeq = ''
barcodeListPerRead = []
fastqLineListPerIdenticalRead = []
for line in filein:
    fastqSeqList = line.strip().split('\t')
    sequence = fastqSeqList[1]
    sequenceLength = len(sequence)
    uncutSequence = fastqSeqList[0].split(' ')[-1]
    barcode = getBarcode(sequence, uncutSequence)
    if sequence != prevSeq:
        if prevSeq != '':
            reduceArtifacts(barcodeListPerRead, fastqLineListPerIdenticalRead, out)
        barcodeListPerRead = []
        fastqLineListPerIdenticalRead = []
    barcodeListPerRead.append(barcode)
    fastqLineListPerIdenticalRead.append(line)
    prevprevSeq = prevSeq
    prevSeq = sequence
else:
        reduceArtifacts(barcodeListPerRead, fastqLineListPerIdenticalRead, out)


out.close()

