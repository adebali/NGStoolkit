import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='expand windows up to a defined size')
parser.add_argument('-i', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input')
parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output')
parser.add_argument('-l', type=int, required= True, help='minimum length')
args = parser.parse_args()


filein = args.i
out = args.o
length_cutoff = args.l
d = {}

def getLength(row):
    start = int(row[1])
    end = int(row[2])
    length = end-start
    return length

def getCount(row):
    count = float(row[6])
    return count

def getName(row):
    name = row[3]
    return name

def printRow(l):
    length = 0
    count = 0
    for e in l:
        length += getLength(e)
        count += getCount(e)
    name = getName(l[0])
    return '\t'.join([name, str(length), str(count)])


def expand(l):
    # chr1    559200  559600  ERD_DnaseU      1000    .       6
    # chr1    659200  659600  ERD_DnaseU      1000    .       7
    # chr1    759200  759600  ERD_DnaseU      1000    .       5

    newList = []
    currentLength = 0
    for row in l:
        print('$')
        print(row)
        print('%')
        length = getLength(row)
        if length < length_cutoff:
            newList.append(row)
            currentLength += getLength(row)
            if currentLength > length_cutoff:
                    return printRow(newList)
                newList = []
        else:
            return printRow([row])
    return printRow(newList)


for row in csv.reader(filein, delimiter='\t'):
    # chr1    559200  559600  ERD_DnaseU      1000    .       6
    name = row[3]
    d[name] = d.get(name, [])
    d[name].append(row)
    if type(row)!=list:
        raise ValueError('hey' + row)

for name in d.keys():
    print(name)
    l = d[name]
    print('X')
    print(l)
    newList = expand(l)
    for e in newList:
        out.write(expand(e) + '\n')

# python expandWindows.py -i XR/dataDir/1801/XR_CPD_Async_B_Minus.inReCh.txt -l 10000