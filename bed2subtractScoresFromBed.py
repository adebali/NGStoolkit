#!/usr/bin/env python
import os
import sys
import argparse
import tempfile
from itertools import izip
from bed import bedline


parser = argparse.ArgumentParser(description='subtraction of the scores between two files: b - a')
parser.add_argument('-a', required= False, help='a of b - a')
parser.add_argument('-b', required= True, help='b of b - a')
parser.add_argument('-o', required= True, help='output')

args = parser.parse_args()

if not args.a:
	os.system("cp " + args.b + " " + args.o)
	sys.exit()

afile = open(args.a, 'r')
bfile = open(args.b, 'r')
out = open(args.o, 'w')

for lineB, lineA in izip(bfile, afile):
	llA = lineA.strip().split('\t')
	stringA = '_'.join([llA[0], llA[1], llA[2], llA[5]])

	llB = lineB.strip().split('\t')
	stringB = '_'.join([llB[0], llB[1], llB[2], llB[5]])

	while stringB != stringA:
		BL_A = bedline(lineA)
		BL_B = bedline(lineB)
		llA = lineA.strip().split('\t')
		print(llA)
		stringA = '_'.join([llA[0], llA[1], llA[2], llA[5]])

		llB = lineB.strip().split('\t')
		stringB = '_'.join([llB[0], llB[1], llB[2], llB[5]])
			
		
		walkA = False
		walkB = False

		if BL_A.chromosome() > BL_B.chromosome():
			walkB = True
			walkA = False
		elif BL_A.chromosome() < BL_B.chromosome():
			walkB = False
			walkA = True
		else: # BL_A.chromosome() == BL_B.chromosome()
			if int(BL_A.start()) > int(BL_B.start()):
				walkB = True
				walkA = False
			elif int(BL_A.start()) < int(BL_B.start()):
				walkB = False
				walkA = True
			else: # int(BL_A.start()) == int(BL_B.start())
				if int(BL_A.end()) > int(BL_B.end()):
					walkB = True
					walkA = False
				elif int(BL_A.end()) < int(BL_B.end()):
					walkB = False
					walkA = True
				else:
					if BL_A.strand() > BL_B.strand():
						walkB = True
						walkA = False
					elif BL_A.strand() < BL_B.strand():
						walkB = False
						walkA = True
					else:
						if stringA != stringB:
							raise ValueError('problem')
		if walkA and walkB:
			raise ValueError('tow files cannot be walked at the same time')
		if walkA:
			try:
				lineA = next(afile)
			except:
				lineA = 'zzz\t9999999999\t9999999999\t.\t99999999999\t+'
		elif walkB:
			out.write(lineB.strip() + '\n')
			try:		
				lineB = next(bfile)
			except:
				out.close()
				sys.exit()
	if stringA != stringB:
		raise ValueError('two strings are not equal: ' + stringA + ' ' + stringB)
	BL_A = bedline(lineA)
	BL_B = bedline(lineB)
	newScore = float(BL_B.score()) - float(BL_A.score())
	fields = BL_B.fields()
	fields[4] = str(newScore)
	out.write('\t'.join(fields) + '\n')

	