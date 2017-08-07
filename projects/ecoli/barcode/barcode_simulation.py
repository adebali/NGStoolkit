#!/usr/bin/env python
import sys
import random
import operator
from collections import Counter 

n = sys.argv[1]
output = 'simulation_' + n + '.txt'
out= open(output, 'w')
def getRandomBarcode():
    barcodeLength = 6
    nucleotideList = ['A', 'T', 'G', 'C']
    return ''.join([random.choice(nucleotideList) for _ in range(barcodeLength)])

def getPoolOfBarcodes(n):
    l = []
    for i in range(n):
        l.append(getRandomBarcode())
    return l

commonDir = {}
for i in range(10000):
    mostCommon = 0
    pool = getPoolOfBarcodes(i)
    poolLength = len(pool)
    # sorted_pool = sorted(pool.items(), key=operator.itemgetter(1), reverse=True)
    # maxCount = sorted_pool[0][1]
    mostCommonList = Counter(pool).most_common(1)
    if len(mostCommonList) != 0:
        mostCommon = mostCommonList[0][1]

    commonDir[mostCommon] = min(commonDir.get(mostCommon, poolLength), poolLength)
    # print(str(poolLength) + '\t' + str(len(set(pool))) + '\t' + str(mostCommon))
for key in sorted(commonDir.keys()):
    out.write(str(key) + '\t' + str(commonDir[key]) + '\n')



