import os
import sys
import pysam
import inspect
from collections import defaultdict

f = "dataDir/0521/NC-B17N.1.cu.bo.hg19.sam"
# f = "dataDir/0106/NC-B17N.1.cu.bo.hg19.sam"
samfile = pysam.AlignmentFile(f, "rb" )

nameDict = defaultdict(int)
total = 0
unmapped = 0
for r in samfile:
    # print(r.query_name)

    # print(r.is_unmapped)
    if (r.is_paired and r.is_read1) or (not r.is_paired):
        nameDict[r.query_name] += 1
        total += 1
        if r.is_unmapped:
            unmapped += 1
        if total == 10000:
            break

countDict = defaultdict(int)
for name in nameDict.keys():
    countDict[nameDict[name]] += 1



print("Total:" + str(total))
print("Total unique:" + str(len(nameDict.keys())))
print("Unmapped:" + str(unmapped))
print("singletons:" + str(countDict[1]))
print("mapped 2 times:" + str(countDict[2]))

    # print(r.is_read1)
    # print(r.is_read2)
    # print(r.__class__.__dict__)
    # break
samfile.close()