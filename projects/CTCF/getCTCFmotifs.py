import os
import sys
import re
from bed import bedline

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
        
TFlist = []
TFs = ['1 CTCF']
for TFline in TFs:
    ll = TFline.strip().split(' ')
    TFname = ll[1]
    count = ll[0]
    pfm = pfmDict.get(TFname, None)
    if pfm != None:
        TFlist.append(TFname)
        # if int(count) > 10000:
        #     print('>' + TFname + '\t' + count)
        #     print(pfm)
fastaInput = "/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.fa"
# out = open("wgEncodeRegTfbsClusteredWithCellsV3.GM12878.strand.selectedExactTFsites_4.bed", "w")
# fastaOut = open("wgEncodeRegTfbsClusteredWithCellsV3.GM12878.strand.selectedExactTFsites_4.fa", "w")
input = sys.argv[1]
output = sys.argv[2]
fastaOutput = output + '.fa'
out = open(output, 'w')
fastaOut = open(fastaOutput, 'w')
TFbedFile = open(input, "r")
# TFbedFile = open("wgEncodeRegTfbsClusteredWithCellsV3.GM12878.strand.bed", "r")
i = 0
for line in TFbedFile:
    if i%1000 == 0:
        print(i)
    i += 1
    bedLine = bedline(line)
    pfm = pfmDict['CTCF']
    newLines = bedLine.getNewLinesWithPfm(fastaInput, pfm, True)
    if len(newLines) < 3:
        for line in newLines:
            singleBedLine = bedline(line)
            fastaOut.write(singleBedLine.getFasta(fastaInput, True, True))
            flankingBedLine = singleBedLine.centerAndExtent(500)
            out.write(flankingBedLine + "\n")
out.close()
