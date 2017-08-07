import os


path = "/proj/sancarlb/users/ogun/seq/mm10/geneList.bed"

filein = open(path, 'r')
for line in filein:
    ll = line.replace(",","_").strip().split("\t")
    newList = [ll[0], ll[1], ll[2], ll[4], '.', ll[3]]
    print("\t".join(newList))