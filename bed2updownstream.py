#!/usr/bin/env python
import bed
import sys
import argparse

parser = argparse.ArgumentParser(description='get twoo flanking regions based on the average length of the main intervals')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('--average', type=bool, default=True, help='Flag if we want to extend at the average length of the original bed file (default=True)')
parser.add_argument('--fixed', action='store_true', default=False, help='fixed region length flanking')
parser.add_argument('-l', type=int, default=1000, help='length of the flanking regions')



args = parser.parse_args()

bedFile = args.i
lines = []
for line in bedFile:
    lines.append(line)

out = args.o

bedObject = bed.bed(None, lines=lines)

if args.fixed:
    flankingLength = args.l
elif args.average:
    averageLength = round(bedObject.getAverageLength())
    flankingLength = averageLength
else:
    raise ValueError("Wrong input")

for line in lines:
    bedLine = bed.bedline(line)
    start1 = int(bedLine.start() - flankingLength)
    end1 = int(bedLine.start())
    
    start2 = int(bedLine.end())
    end2 = int(bedLine.end() + flankingLength)

    upstreamLine = None
    downstreamLine = None

    upstreamLine = downstreamLine = bedLine
    if bedLine.strand() == '+' or bedLine.strand() == '.':
        upstreamLine = upstreamLine.changeField(2, start1)
        upstreamLine = upstreamLine.changeField(3, end1)
        downstreamLine = downstreamLine.changeField(2, start2)
        downstreamLine = downstreamLine.changeField(3, end2)
    elif bedLine.strand() == '-': 
        upstreamLine = upstreamLine.changeField(2, start2)
        upstreamLine = upstreamLine.changeField(3, end2)
        downstreamLine = downstreamLine.changeField(2, start1)
        downstreamLine = downstreamLine.changeField(3, end1)
        # downstreamLine = bedLine.newline(start1, end1)
        # upstreamLine = bedLine.newline(start2, end2)
    else:
        raise ValueError('strand is not what we expected: ' + bedLine.strand())
    

    upstreamLineName = upstreamLine.name()
    if '_upstream' in upstreamLineName:
        raise ValueError('_upstream word cannot be in the name')
        upstreamLineName = upstreamLine.name()
    
    downstreamLineName = downstreamLine.name()
    if '_downstream' in downstreamLineName:
        raise ValueError('_downstream word cannot be in the name')

    upstreamLine = upstreamLine.changeField(4, upstreamLineName + '_upstream')
    downstreamLine = downstreamLine.changeField(4, downstreamLineName + '_downstream')
    # upBedLine = bed.bedline(upstreamLine)
    # line1 = upBedLine.newline(upBedLine.start(), upBedLine.end(), {4: 'upstream'})


    # downBedLine = bed.bedline(downstreamLine)
    # line2 = downBedLine.newline(downBedLine.start(), downBedLine.end(), {4: 'downstream'})


    # upstreamFields = bed.bedline(upstreamLine).fields()
    # upstreamFields.append('upstream')
    # downstreamFields = bed.bedline(downstreamLine).fields()
    # downstreamFields.append('downstream')
    # line1 = '\t'.join(upstreamFields)
    # line2 = '\t'.join(downstreamFields) 

    out.write(upstreamLine.getLine() + '\n')
    out.write(downstreamLine.getLine() + '\n')