#!/usr/bin/env python
import os
import sys

def getFirstReadsOnly(file):
	firstReadIdentifier = '/1'
	standAloneList = []
	filein = open(file, 'r')
	for line in filein:
		ll = line.split('\t')
		name = ll[3].strip()
		originalReadName = name[:-2]
		if name.endswith(firstReadIdentifier):
			print(line.strip())
		if originalReadName in standAloneList:
			standAloneList.remove(originalReadName)
		else:
			standAloneList.append(originalReadName)
	if not standAloneList == []:
		sys.exit('This is not a paired end data, exiting...')
	else:
		return True

bedFile = sys.argv[1]
getFirstReadsOnly(bedFile)
