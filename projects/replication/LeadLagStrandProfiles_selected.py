#!/usr/bin/env python

import argparse
import os
import sys
import tempfile

PARSER = argparse.ArgumentParser(description='removes neighboring intervals')
PARSER.add_argument('-i', required=True, help='input')
PARSER.add_argument('-a', required=True, help='bed A file')
PARSER.add_argument('-o', required=True, help='output')
PARSER.add_argument('-s', required=True, type=float, help='scaleFactor')
PARSER.add_argument('-g', required=True, type=str, help='genome file')
PARSER.add_argument('-d', required=False, default=1000000, type=int, help='distance (one side)')
PARSER.add_argument('-cutoff', required=False, default=900, type=int, help='cutoff for the called peak/region')

args = PARSER.parse_args()

distance = args.d
cutoff = args.cutoff
genomeFile = args.g

tempfile1_, tempfile1 = tempfile.mkstemp()
tempfile2_, tempfile2 = tempfile.mkstemp()
tempfile3_, tempfile3 = tempfile.mkstemp()
tempfile4_, tempfile4 = tempfile.mkstemp()

def run(code):
    print(code)
    if os.system(code):
        raise ValueError(code + ' cannot be executed')

code = 'bed2removeNeighbors.py -d ' + str(distance) + ' -i ' + args.a + ' >' + tempfile1 
run(code)

code = 'bed2midPointFlankingBed.py -i ' + tempfile1 + ' -w ' + str(distance) + ' -g ' + genomeFile + ' -o ' + tempfile2
run(code)

oneSideWindowNumer = 100

code = 'bed2makeWindows.py -i ' + tempfile2 + ' -w ' + str(distance/oneSideWindowNumer) + ' | cut -f 1-6 '
code += ' | bedtools intersect -a \'stdin\' -b ' + args.i + ' -c > ' + tempfile3
run(code)

windowNumer = oneSideWindowNumer * 2
countCol = 7
code = 'bedCount2normalizedCount.py -i ' + tempfile3 + ' -m ' + str(args.s) + ' -l 1000 -c ' + str(countCol) + ' -o' + tempfile4
run(code)

code = 'bedCountWindows2positions.py -i' + tempfile4 + ' -wnum ' + str(windowNumer) + ' -countCol ' + str(countCol) + ' -o ' + args.o
run(code)

run('rm ' + tempfile1)
run('rm ' + tempfile2)
run('rm ' + tempfile3)
run('rm ' + tempfile4)
