#!/usr/bin/env python
import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='pastes last column of files while keeping the all colums from the first file')
parser.add_argument('-inputs', nargs='+', help='<Required> Input file list (min 2)', required=True)
parser.add_argument('-n', default=1, type=int, help='last N number of columns to merge', required=False)
# parser.add_argument('-o', required= True, help='<Required> Output')

def n2str(n):
    s = '1'
    for i in range(2, n+1):
        s += ',' + str(i)
    return s

args = parser.parse_args()
inputs = args.inputs
# output = args.o


code = ['bash', '-c']
pasteString = 'paste <(cat ' + os.path.realpath(inputs[0]) + ') '
for i in range(1,len(inputs)):
    input = inputs[i]
    pasteString += '<(cat ' + os.path.realpath(input) + ' | rev | cut -f ' + n2str(args.n) + ' | rev) '
code.append(pasteString)
# code.append('>')
# code.append(os.path.realpath(output))

codeString = ' '.join(code)
# print(codeString)
# print(code)
# os.system(codeString)
subprocess.call(code)