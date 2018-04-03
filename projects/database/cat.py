import os
import json
import slurm
import argparse
from conf import samples

parser = argparse.ArgumentParser(description='cat master pipeline')
parser.add_argument('-json', default='sample.json', required= False, help='file for sample configuration in JSON format')
parser.add_argument('-e', required= True, type=int, help='Experiment number')
parser.add_argument('-s', required= False, type=int, nargs='+', help='space separated sample list')
parser.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')

args = parser.parse_args()

sampleClass = samples(args.json)
if args.s:
    sampleNoList = args.s
else:
    sampleNoList = sampleClass.experiment2sampleNoList(args.e)

code = 'rm -f mergedGeneCounts.txt'
os.system(code)

for n in sampleNoList:
    code = 'python pipeline.py cat -e ' + str(args.e) + ' -n ' + str(n)
    print(code)
    os.system(code)



