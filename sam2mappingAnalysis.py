#! /usr/bin/env python
import os
import sys
import pysam
import inspect
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='analyze sam for mappability')
# parser.add_argument('-i', required=True, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
args = parser.parse_args()

# fileInput = args.i
out = args.o

fileInput = "/proj/sancarlb/users/ogun/scripts/projects/DamageSeq/hiSeq/dataDir/0521/NC-B17N.1.cu.bo.hg19.sam"
# fileInput = "/proj/sancarlb/users/ogun/scripts/projects/DamageSeq/hiSeq/dataDir/0106/NC-B17N.1.cu.bo.hg19.sam"
# f = "dataDir/0106/NC-B17N.1.cu.bo.hg19.sam"
samfile = pysam.AlignmentFile(fileInput, "rb" )

nameDict = defaultdict(int)
total = 0
unmapped = 0
for r in samfile:
    # print(r.query_name)

    # print(r.is_unmapped)
    if (r.is_paired and r.is_read1) or (not r.is_paired):
        nameDict[r.query_name] += 1
        # print(r.query_name)
        total += 1
        if r.is_unmapped:
            unmapped += 1

        if r.is_duplicate:
            raise ValueError("Duplicate: " + r.query_name)
        if r.is_secondary:
            raise ValueError("Secondary: " + r.query_name)
        if total == 2000000:
            break

countDict = defaultdict(int)
for name in nameDict.keys():
    countDict[nameDict[name]] += 1



out.write("Total:" + str(total) + "\n")
out.write("Total unique:" + str(len(nameDict.keys())) + "\n")
out.write("Unmapped:" + str(unmapped) + "\n")
out.write("singletons:" + str(countDict[1]) + "\n")
out.write("mapped 2 times:" + str(countDict[2]) + "\n")

    # print(r.is_read1)
    # print(r.is_read2)
    # print(r.__class__.__dict__)
    # break
samfile.close()