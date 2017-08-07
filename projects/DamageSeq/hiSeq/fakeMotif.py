import os
import sys
import re
from bed import bedline



out = open("/nas/longleaf/home/adebali/ogun/ENCODE/Txn/inputDNA210_noTFregion_TFmotifs_unique_1K.bed", "w")
filein = open("/nas/longleaf/home/adebali/ogun/ENCODE/Txn/inputDNA210_noTFregion_TFmotifs_unique.bed", "r")
for line in filein:
    bedLine = bedline(line)
    out.write(bedLine.centerAndExtent(500) + "\n")
out.close()