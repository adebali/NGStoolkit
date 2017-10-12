#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

# Uses chr as chromosome identifier.
parser = argparse.ArgumentParser(description='converts gff exon format to gene list (TSS to TES) in bed format')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')

# 1	unknown	exon	3631	3913	.	+	.	gene_id "NAC001"; gene_name "NAC001"; p_id "P20202"; transcript_id "NM_099983.2"; tss_id "TSS18959";
# 1	unknown	CDS	3760	3913	.	+	0	gene_id "NAC001"; gene_name "NAC001"; p_id "P20202"; transcript_id "NM_099983.2"; tss_id "TSS18959";

args = parser.parse_args()
gtfFile = args.i

filein = open(gtfFile, 'r')
out = args.o

masterTSdict = {}

for line in filein:
    ll = line.split('\t')
    region = ll[2]
    # if region == "exon":
    if True:
        chromosome = ll[0]
        start = int(ll[3])
        end = int(ll[4])
        strand = ll[6]
        info = ll[8]
        info_units = info.strip().split(';')
        infoDict = {}
        for info_unit in (x for x in info_units if x != ''):
            l = info_unit.strip().split('"')
            infoDict[l[0].strip()] = l[1].strip()
            # break
        # print(infoDict)
        gene_id = infoDict['gene_id'].replace('"','')
        # print(gene_id)
        currentGeneid = masterTSdict.get(gene_id, {'start': start, 'end': end})
        # print(currentGeneid)
        masterTSdict[gene_id] = {
            'chromosome': chromosome,
            'start': min(start, currentGeneid['start']), 
            'end': max(end, currentGeneid['end']),
            'strand': strand
            }
        # print(masterTSdict)
for gene_id in masterTSdict.keys():
    dictionary = masterTSdict[gene_id]
    chromosome = dictionary['chromosome']
    start = dictionary['start']
    end = dictionary['end']
    strand = dictionary['strand']
    out.write(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + gene_id + '\t.\t' + strand + '\n')
out.close()
