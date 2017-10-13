#!/usr/bin/env python
import os
import sys
import bed
import generalUtils
import argparse
import random

parser = argparse.ArgumentParser(description='count strings')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')

args = parser.parse_args()

bedFile = args.i
output = args.o

os.system('uniq -c ' + bedFile + ' | cut -f 1 | awk -F" chr" \'{print $1}\' | sed \'s/ *$//\' | sort -k1,1n | uniq -c | unexpand -a > ' + output)

# cutOff = 1000
# code = 'uniq -c ' + bedFile + ' | sed -r \'s/^( *[^ ]+) +/\\1\\t/\' | awk \' $1 >= ' + str(cutOff) + ' \' > ' + output
# print(code)