import os
import sys
import unittest

class bed:
	def __init__(self, input):
		self.file = input

	def getTotalRegionLength(self):
		totalLength = 0
		filein = open(self.file, 'r')
		for line in filein:
			ll = line.split('\t')
			beg = int(ll[1])
			end = int(ll[2])
			regionLength = end - beg
			totalLength += regionLength
		return totalLength

	def getAverageLength(self):
		totalLength = self.getTotalRegionLength()
		hitNum = self.getHitNum()
		return float(totalLength)/hitNum

	def getHitNum(self):
		hitNum = 0
		filein = open(self.file, 'r')
		for line in filein:
			hitNum += 1
		return hitNum

	def fixRange(self, side, length):
		filein = open(self.file, 'r')
		for line in filein:
			print(bedLine2fixedRangedLine(line, side, length))
		return self

class bedpe(bed):
	def mergeFragments(self):
		filein = open(self.file, 'r')
		for line in filein:
			print(bedpeLine2bedLine(line))
		return self


def bedpeLine2bedLine(bedpeLine):
	# Converts paired bed line to single fragment by getting the maximum range between to genomic locations.
	ll = bedpeLine.strip().split('\t')
	chr1 = ll[0]
	beg1 = int(ll[1])
	end1 = int(ll[2])
	chr2 = ll[3]
	beg2 = int(ll[4])
	end2 = int(ll[5])
	name = ll[6]
	score = ll[7]
	strand1 = ll[8]
	strand2  = ll[9]

	if chr1 != chr2:
		sys.exit('Chromosome is not identical between mates! Exiting...')

	if strand1 == '+':
		beg = beg1
		end = end2
	elif strand1 == '-':
		beg = beg2
		end = end1

	newLine = chr1 + \
			'\t' + str(beg) + \
			'\t' + str(end) + \
			'\t' + name + \
			'\t' + score + \
			'\t' + strand1

	if len(ll) > 10:
		newLine += '\t' + '\t'.join(ll[10:])

	return newLine

def bedLine2fixedRangedLine(bedLine, side, length, strandFlag=True):
	ll = bedLine.strip().split('\t')
	beg = int(ll[1])
	end = int(ll[2])
	strand = ll[5]
	if (strandFlag == True and strand == '+') or (strandFlag == False):
		if side == 'l':
			newEnd = min(beg + length, end)
			newBeg = beg
		elif side == 'r':
			newBeg = max(end - length, beg)
			newEnd = end
	elif strandFlag == True and strand == '-':
		if side == 'l':
			newBeg = max(end - length, beg)
			newEnd = end
		elif side == 'r':
			newEnd = min(beg + length, end)
			newBeg = beg
	if side != 'l' and side != 'r':
		sys.exit('Error: side must be either l or r. Exiting...')
	ll[1] = str(newBeg)
	ll[2] = str(newEnd)
	return '\t'.join(ll)


class bedTests(unittest.TestCase):

	def test_bedpeLine2bedLine(self):
		bedpeLine = 'chr1\t100\t200\tchr1\t5000\t5100\tbedpe_example1\t30\t+\t-\n'
		expectedResult1 = 'chr1\t100\t5100\tbedpe_example1\t30\t+'
		self.assertEqual(bedpeLine2bedLine(bedpeLine), expectedResult1)
		bedpeLineMinus = 'chr1\t1000\t1200\tchr1\t900\t950\tbedpe_example2\t30\t-\t+\n'
		expectedResult2 = 'chr1\t900\t1200\tbedpe_example2\t30\t-'
		self.assertEqual(bedpeLine2bedLine(bedpeLineMinus), expectedResult2)
		bedpeLineCustomFields = bedpeLine.strip() + '\tField1\tField2\n'
		expectedResult3 = expectedResult1 + '\tField1\tField2'
		self.assertEqual(bedpeLine2bedLine(bedpeLineCustomFields), expectedResult3)


	def test_bedLine2fixedRangeLine(self):
		bedLine = 'chr1\t100\t5100\tbedpe_example1\t30\t+\n'
		expectedResult = 'chr1\t100\t120\tbedpe_example1\t30\t+'
		self.assertEqual(bedLine2fixedRangedLine(bedLine, 'l', 20), expectedResult)
		expectedResult = 'chr1\t5075\t5100\tbedpe_example1\t30\t+'
		self.assertEqual(bedLine2fixedRangedLine(bedLine, 'r', 25), expectedResult)

		bedLine2 = 'chr1\t100000212\t100000296\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t+'
		expectedResult2 = 'chr1\t100000212\t100000222\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t+'
		self.assertEqual(bedLine2fixedRangedLine(bedLine2, 'l', 10), expectedResult2)

		bedLine3 = 'chr1\t100\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
		expectedResult3 = 'chr1\t140\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
		self.assertEqual(bedLine2fixedRangedLine(bedLine3, 'l', 10), expectedResult3)

		bedLine4 = 'chr1\t100\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
		expectedResult4 = 'chr1\t100\t110\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
		self.assertEqual(bedLine2fixedRangedLine(bedLine3, 'r', 10), expectedResult4)

	def test_bedClass(self):
		myBed = bed(os.path.join(os.path.dirname(os.path.realpath(__file__)),'testFiles/bedExample.bed'))
		myBed.fixRange('l', 10)

if __name__ == "__main__":
	unittest.main()
