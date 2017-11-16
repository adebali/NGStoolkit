import os
import sys

for line in open(sys.argv[1]):
    ll = line.strip().split('\t')
    ll[0] = ll[0].split('.')[0]
    newLine = '\t'.join(ll)
    print(newLine)