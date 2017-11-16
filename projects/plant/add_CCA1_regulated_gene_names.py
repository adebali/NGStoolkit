import os
import sys

d = {}
filein = open('/Users/ogunadebali/SancarLab/plant/CCA1_regulated_geneList.txt')
for line in filein:
    if line.strip():
        ll = line.strip().split('\t')
        d[ll[0]] = ll[1]

out = open('/Users/ogunadebali/SancarLab/plant/meta2d_CAC1high_TS_dcasted_lt0.05_geneNames.csv', 'w')
filein = open('/Users/ogunadebali/SancarLab/plant/meta2d_CAC1high_TS_dcasted_lt0.05.csv')
for line in filein:
    ll = line.strip().split(',')
    geneName = d.get(ll[1].replace('"', ''), "NA")
    ll.append(geneName)
    out.write(','.join(ll) + '\n')

out.close()