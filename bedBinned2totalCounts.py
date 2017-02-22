#!/usr/bin/env python
import sys
import bed
import argparse
from operator import add

parser = argparse.ArgumentParser(description='converts binned bed file into total count per window')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('--mergeStrands', action='store_true', help='give out two lines for each strand')
parser.add_argument('-reverseStrand', required=False, help='reverse the specified strand')
parser.add_argument('-n', required=True, help='window number for each interval')
args = parser.parse_args()

filein = args.i
out = args.o
windowNumber = int(args.n)

strands = { '+' : [0] * windowNumber, '-': [0] * windowNumber}
i = 0
for line in filein:
	position = i%windowNumber
	bedLine = bed.bedline(line)
	strand = bedLine.strand()
	count = int(bedLine.fields()[-1])
	strands[strand][position] += count
	i += 1

if args.reverseStrand:
	strand = args.reverseStrand
	strands[strand] = list(reversed(strands[strand]))

if args.mergeStrands:
	lists = [map(add, strands['+'], strands['-'])]
else:
	lists = [strands['+'], strands['-']]

print(strands)
for e in lists:
	for positionCount in e:
		out.write(str(positionCount) + '\t')
	out.write('\n')

