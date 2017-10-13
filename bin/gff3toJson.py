#!/usr/bin/env python
import os
import json
import argparse

parser = argparse.ArgumentParser(description='converts gff3 format to json')
parser.add_argument('-i', required=True, help='<Required> input')
parser.add_argument('-o', required=True, help='<Required> output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
features = []

def attributes2object(attributes):
    myDict = {}
    al = attributes.split(';')
    for tagField in al:
        tl = tagField.split('=')
        tag = tl[0]
        tagValue = tl[1]
        if ',' in tagValue:
            myDict[tag] = []
            for tagValueItem in tagValue.split(','):
                if ':' in tagValueItem:
                    l = tagValueItem.split(':')
                    key = l[0]
                    value = ':'.join(l[1:])
                    myDict[tag].append({key: value})
                else:
                    myDict[tag].append(tagValueItem)
        elif ':' in tagValue:
            l = tagValue.split(':')
            key = l[0]
            value = ':'.join(l[1:])
            myDict[tag] = {key: value}
        else:
            myDict[tag] = tagValue
    return myDict

# print(attributes2object('ID=id3;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2'))

filein = open(inputFile, 'r')
i = 0
for line in filein:
    i += 1
    if not line.startswith('#'):
        ll = line.strip().split('\t')
        feature = {
            'seqid'     : ll[0],
            'source'    : ll[1],
            'type'      : ll[2],
            'start'     : int(ll[3]),
            'end'       : int(ll[4]),
            'score'     : str(ll[5]) if ll[5] == '.' else float(ll[5]),
            'strand'    : ll[6],
            'phase'     : str(ll[7]) if ll[7] == '.' else int(ll[7]),
            'attributes': attributes2object(ll[8])
        }
        if 'gbkey' in feature['attributes'].keys() and \
            feature['attributes']['gbkey'] == 'Gene' and \
            feature['attributes']['gene_biotype'] != 'pseudogene' and \
            feature['attributes']['gene_biotype'] != 'misc_RNA' \
            :
            features.append(feature)
    # if i == 1000:
    #     break
# print(len(features))
with open(outputFile, 'w') as fp:
    json.dump(features, fp, indent=4)