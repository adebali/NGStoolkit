#!/usr/bin/env python
import os
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='slice bed based on the score order')
parser.add_argument('-i', required= True, help='<Required> input')
parser.add_argument('-n', required= False, default=4, type=int, help='total n pieces (default: 4)')
parser.add_argument('-slice', required= True, type=int, help='slice number from big to small')

args = parser.parse_args()
bedFile =  args.i
totalPieces = args.n
sliceNum = args.slice

totalHitNumber = bed.bed(bedFile).getHitNum()

code = 'sort -r -k5,5n ' + bedFile + ' | awk "NR>=' + str((sliceNum-1)*float(totalHitNumber/totalPieces)) +' && NR<' + str(sliceNum*float(totalHitNumber/totalPieces)) + '" | sort -k1,1 -k2,2n -k3,3n'
# print(code)
os.system(code)