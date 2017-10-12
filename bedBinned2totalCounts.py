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
parser.add_argument('-start', default=0, type=int, required=False, help='start of the relative position')
parser.add_argument('-winLength', default=1, type=int, required=False, help='window length')
parser.add_argument('--writePosition', action='store_true', help='writes out the position as well, should be used with --mergeStrands')
parser.add_argument('--average', action='store_true', help='calculates the average instead of the total count')


args = parser.parse_args()

filein = args.i
out = args.o
windowNumber = int(args.n)
windowLength = args.winLength
startPosition = args.start

strands = { '+' : [0] * windowNumber, '-': [0] * windowNumber}
i = 0
geneNo = 0
for line in filein:
	position = i%windowNumber
	if position == 0:
		geneNo += 1
	bedLine = bed.bedline(line)
	strand = bedLine.strand()
	if strand == '.':
		strand = '+'
	count = float(bedLine.fields()[-1])
	strands[strand][position] += count
	i += 1

if args.reverseStrand:
	strand = args.reverseStrand
	strands[strand] = list(reversed(strands[strand]))

if args.mergeStrands:
	lists = [map(add, strands['+'], strands['-'])]
	separator = "\n"
else:
	lists = [strands['+'], strands['-']]
	separator = "\t"

print(strands)
for e in lists:
	i = 0
	for positionCount in e:
		if args.average:
			positionCount = float(positionCount)/ geneNo
		position = startPosition + (i * windowLength)
		i += 1
		line = ''
		if args.mergeStrands and args.writePosition:
			line = str(position) + '\t'
		line += str(positionCount) + separator
		out.write(line)
	if not args.mergeStrands:
		out.write('\n')

