import os
import sys
import re

pfms = open("/proj/sancarlb/users/ogun/ENCODE/Txn/jaspar_pfm_vertebrates.txt", "r").read().split("\n\n")
pfmDict = {}
for pfm in pfms:
    if pfm.startswith(">"):
        lines = pfm.split("\n")
        header = lines[0]
        hl = header.strip().split("\t")
        id = hl[0]
        name = hl[1]
        pfmString = "\n".join(lines[1:])
        pfmDict[name] = pfmString
        # print(pfmString)
        

TFs = open("/proj/sancarlb/users/ogun/ENCODE/Txn/GM12878_TFlist.txt", "r").readlines()
for TFline in TFs:
    ll = TFline.strip().split(' ')
    TFname = ll[1]
    count = ll[0]
    pfm = pfmDict.get(TFname, None)
    if pfm != None:
        if int(count) > 10000:
            print('>' + TFname + '\t' + count)
            print(pfm)