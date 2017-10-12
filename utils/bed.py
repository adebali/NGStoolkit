import os
import sys
import unittest
import random
import pybedtools
import math
from pfm import Pfm

class bedline:
	def __init__(self, line, sep='\t'):
		self.line = line.strip()
		self.separator = sep

	def getLine(self):
		return self.line

	def fields(self):
		return self.line.split(self.separator)

	def fields2line(self, fields, sep='\t'):
		return sep.join(str(e) for e in fields)

	def getField(self, columNo):
		return self.fields()[columNo - 1]

	def chromosome(self, chromosomeCol = 1):
		return self.getField(chromosomeCol)

	def start(self, startCol = 2):
		return int(self.getField(startCol))

	def length(self):
		return self.end() - self.start() 

	def end(self, endCol = 3):
		return int(self.getField(endCol))

	def name(self, nameCol = 4):
		return self.getField(nameCol)

	def score(self, scoreCol = 5):
		return float(self.getField(scoreCol))

	def strand(self, strandCol = 6):
		return self.getField(strandCol)

	def reverseStrand(self, strand):
		if strand == "+":
			return "-"
		elif strand == "-":
			return "+"
		else:
			raise ValueError("Unexpected strand information.")

	def midpoint(self, randomness=False):
		if not randomness:
			addition = 0
		else:
			addition = random.sample([0,1], 1)[0]
		return int((self.start() + self.end() + addition) / 2)

	def centerAndExtent(self, flankingLength):
		'''gets the center point and extent the interval to both directions for each is flankingLength'''
		center = self.midpoint()
		strand = self.strand()
		start = center - flankingLength
		end = center + flankingLength
		if strand == "+":
			end += 1
		elif strand == "-":
			start -= 1
		return self.newline(start, end)

	def newline(self, start, end, otherFields = {}):
		newList = self.fields()
		newList[1] = str(start)
		newList[2] = str(end)
		if otherFields != {}:
			for columnNo in otherFields.keys():
				newList[columnNo - 1] = otherFields[columnNo]
		return self.separator.join(newList)

	def addField(self, fields = []):
		newList = self.fields() + fields
		return self.separator.join(newList)

	def getFasta(self, fastaInput, s=True, name=False):
		b = pybedtools.BedTool(self.getLine(), from_string=True)
		s = b.sequence(fi=fastaInput, s=s, name=name)
		return(open(s.seqfn).read())

	def singleFastaToSequence(self, fastaString):
		lines = fastaString.strip().split("\n")
		return("".join(lines[1:]))

	def getSequence(self, fastaInput):
		return self.singleFastaToSequence(self.getFasta(fastaInput)).upper()

	def countString(self, fastaInput, string):
		sequence = self.getSequence(fastaInput)
		return sequence.count(string)

	def getNewLinesWithPfm(self, fastaInput, matrix, asString = False, strandCol = 6):
		currentStrand = self.strand(strandCol)
		if currentStrand == ".":
			currentStrand = "+"
		newBedLine = bedline(self.newline(self.start(), self.end(), {strandCol: currentStrand}))
		sequence = newBedLine.getSequence(fastaInput)

		if asString:
			pfm = Pfm()
			pfm.takeStringInput(matrix)
		else:
			pfm = Pfm(matrix)
		hits = pfm.getHits(sequence)
		
		newLines = []
		if hits[0] != []:
			strand = currentStrand
			for i in hits[0]:
				newLines.append(newBedLine.newline(self.start() + i, self.start() + i + pfm.length, {strandCol: strand}))
		if hits[1] != []:
			strand = newBedLine.reverseStrand(currentStrand)
			for i in hits[1]:
				newLines.append(newBedLine.newline(self.start() + i, self.start() + i + pfm.length, {strandCol: strand}))
		return newLines

	def changeField(self, fieldCol, value):
		fields = self.fields()
		fields[fieldCol - 1] = value
		return bedline(self.fields2line(fields))

class bed:
	def __init__(self, input=None, **kwargs):
		self.file = input
		if self.file == None:
			if kwargs['lines']:
				self.lines = kwargs['lines']

	def getTotalRegionLength(self):
		totalLength = 0
		if self.file:
			filein = open(self.file, 'r')
			for line in filein:
				ll = line.split('\t')
				beg = int(ll[1])
				end = int(ll[2])
				regionLength = end - beg
				totalLength += regionLength
			return totalLength
		elif self.lines:
			for line in self.lines:
				ll = line.split('\t')
				beg = int(ll[1])
				end = int(ll[2])
				regionLength = end - beg
				totalLength += regionLength
			return totalLength
		return None

	def getColumnNumber(self):
		firstLine = open(self.file, 'r').readline().strip()
		return len(firstLine.split('\t'))

	def getFirstLineLength(self):
		with open(self.file, 'r') as f:
			firstBedLine = bedline(f.readline().strip())
			return firstBedLine.length()

	def getAverageLength(self):
		totalLength = self.getTotalRegionLength()
		hitNum = self.getHitNum()
		return float(totalLength)/hitNum

	def getHitNum(self):
		hitNum = 0
		if self.file:
			filein = open(self.file, 'r')
			for line in filein:
				hitNum += 1
			return hitNum
		if self.lines:
			return(len(self.lines))
		return None
	

	def fixRange(self, side, length):
		filein = open(self.file, 'r')
		for line in filein:
			print(bedLine2fixedRangedLine(line, side, length))
		return self

	def lengthDistribution(self):
		filein = open(self.file, 'r')
		myDict = {}
		for line in filein:
			ll = line.split('\t')
			start = int(ll[1])
			end = int(ll[2])
			length = end - start
			if length in myDict.keys():
				myDict[length] += 1
			else:
				myDict[length] = 1
		sortedKeys = sorted(myDict.keys())
		for i in range(min(sortedKeys), max(sortedKeys) + 1):
			if not i in myDict.keys():
				count = 0
			else:
				count = myDict[i]
			print(str(i) + '\t' + str(count))

	def printLengths(self):
		filein = open(self.file, 'r')
		for line in filein:
			ll = line.split('\t')
			start = int(ll[1])
			end = int(ll[2])
			length = end - start
			print(length)

	def makeWindowsPerLine(self, interval, windowSize):
		start = interval[0]
		end = interval[1]
		length = end - start
		windowNumber = length / windowSize
		newIntervals = []
		for i in range(windowNumber):
			newIntervals.append([start + i * windowSize, start + i * windowSize + windowSize])
		return newIntervals

	def makeWindows(self, windowSize, noShortFlag=False):
		filein = open(self.file, 'r')
		for line in filein:
			ll = line.strip().split('\t')
			interval = [int(ll[1]), int(ll[2])]
			newIntervals = self.makeWindowsPerLine(interval, windowSize)
			for newInterval in newIntervals:
				newLineList = list(ll)
				newLineList[1] = str(newInterval[0])
				newLineList[2] = str(newInterval[1])
				printFlag = True
				if noShortFlag:
					newIntervalLength = newInterval[1] - newInterval[0]
					if newIntervalLength < windowSize:
						printFlag = False
				if printFlag:
					print('\t'.join(newLineList))

	def read(self):
		filein = open(self.file, 'r')
		for line in filein:
			yield(bedline(line))

	def removeNeighbors(self, distance):
		def areNeighbors(previousLine, line, distance):
			bedLine = bedline(line)
			previousBedLine = bedline(previousLine)
			if bedLine.chromosome() != previousBedLine.chromosome():
				return False
			start = bedLine.start()
			previousEnd = previousBedLine.end()
			if abs(start - previousEnd) < distance:
				return True
			else:
				return False
		previousLine = "NA\t" + str(-2*distance) + "\t" + str(-1*distance)
		filein = open(self.file, 'r')
		start = True
		dontPrintNextline = False
		for line in filein:
			linesAreNeighbors = areNeighbors(previousLine, line, distance)
			if linesAreNeighbors == False:
				if not start:
					if not previousLineWasNeigbor:
						print(previousLine.strip())
				dontPrintNextline = False
				previousLineWasNeigbor = False
			else:
				previousLineWasNeigbor = True
			previousLine = line
			start = False
		if linesAreNeighbors == False:
			print(line.strip())

	def sumCount(self, columnNo = 7):
		totalCount = 0
		for bedline in self.read():
			totalCount += int(bedline.getField(columnNo))
		return totalCount

class bedpe(bed):
	def mergeFragments(self):
		filein = open(self.file, 'r')
		for line in filein:
			print(bedpeLine2bedLine(line))
		return self


def bedpeLine2bedLine(bedpeLine):
	'''Converts paired bed line to single fragment by getting the maximum range between to genomic locations.'''
	if bedpeLine.strip() == '':
		return ''
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
	
class bedintersect(bedline):
	'''Assumes bed intersect line with 12 columns:
	chr1    100     200     .       .       +       chr1    150     170     .       .       +
	chr1    330     2400    .       .       -       chr1    350     370     .       .       -
	'''
	def chromosome1(self, c=1):
		return self.chromosome(c)
	def chromosome2(self, c=7):
		return self.chromosome(c)
	def start1(self, c=2):
		return self.start(c)
	def end1(self, c=3):
		return self.start(c)
	def start2(self, c=8):
		return self.start(c)
	def end2(self, c=9):
		return self.start(c)
	def strand1(self, c=6):
		return self.strand(c)
	def strand2(self, c=12):
		return self.strand(c)
	def length1(self):
		return self.end1() - self.start1()
	def length2(self):
		return self.end2() - self.start2()
	
	def sameStrands(self):
		if not self.strand1() in ['+', '-', '.'] or not self.strand2() in ['+', '-', '.']:
			raise ValueError('strands are not defined properly: ' + self.strand1() + ' and ' + self.strand2())
		if self.strand1() == self.strand2():
			return True
		return False

	def position2(self, pointOfB='center'):
		if pointOfB == 'center':
			position = (float(self.start2()) + float(self.end2())) / 2
		elif pointOfB == 'start':
			if self.strand1() == '+':
				position = self.start2()
			elif self.strand1() == '-':
				position = self.end2()
			elif self.strand1() == '.':
				position = self.start2()
			else:
				raise ValueError('Strand is not defined')		
		elif pointOfB == 'end':
			if self.strand1() == '+':
				position = self.end2()
			elif self.strand1() == '-':
				position = self.start2()
			elif self.strand1() == '.':
				position = self.end1()
			else:
				raise ValueError('Strand is not defined')		
		else:
			raise ValueError('pointOfB is not defined as expected (center, start or end): ' + pointOfB)
		return position

	def relativeDistanceOfPosition2(self, pointOfB= 'center', relativeTo='start'):
		position = self.position2(pointOfB)
		if relativeTo == 'start':
			if self.strand1() == '+':
				distance = abs(position - self.start1())
			elif self.strand1() == '-':
				distance = abs(position - self.end1())
			else:
				raise ValueError('relativeTo is not defined as expected (start or end): ' + relativeTo)
		if relativeTo == 'end':
			if self.strand1() == '+':
				distance = abs(position - self.end1())
			elif self.strand1() == '-':
				distance = abs(position - self.start1())
			else:
				raise ValueError('relativeTo is not defined as expected (start or end): ' + relativeTo)
		return distance

	def getDistancePercentage(self, pointOfB = 'center', relativeTo = 'start'):
		return round(100*(self.relativeDistanceOfPosition2(pointOfB, relativeTo) / float(self.length1())))

	def getAbsoluteDistance(self, region = 'main', window= 1, pointOfB = 'center'):
		if region == 'main' or region == 'downstream':
			relativeTo = 'start'
		elif region == 'upstream':
			relativeTo = 'end'
		else:
			raise ValueError('Not a valid region: ' + region)
		
		position = self.position2(pointOfB)
		if relativeTo == 'start':
			if self.strand1() == '+':
				distance = abs(position - self.start1())
			elif self.strand1() == '-':
				distance = abs(position - self.end1())
			else:
				raise ValueError('relativeTo is not defined as expected (start or end): ' + relativeTo)
		if relativeTo == 'end':
			if self.strand1() == '+':
				distance = abs(position - self.end1())
			elif self.strand1() == '-':
				distance = abs(position - self.start1())
			else:
				raise ValueError('relativeTo is not defined as expected (start or end): ' + relativeTo)
		if window != 1:
			distance = math.ceil(distance/window)
		if region == 'upstream':
			distance = 0 - distance
		elif region == 'downstream':
			distance = distance + 100
		return distance



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

	# def test_bedClass(self):
	# 	myBed = bed(os.path.join(os.path.dirname(os.path.realpath(__file__)),'testFiles/bedExample.bed'))
	# 	myBed.fixRange('l', 10)

	# 	myBed = bed(os.path.join(os.path.dirname(os.path.realpath(__file__)),'testFiles/genes.bed'))
	# 	myBed.removeNeighbors(20)

	def test_bedline2(self):
		bedLine = bedline('chr1\t100\t5100\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.fields(), ['chr1','100','5100','bedpe_example1','30','+'])
		self.assertEqual(bedLine.chromosome(), 'chr1')
		self.assertEqual(bedLine.start(), 100)
		self.assertEqual(bedLine.end(), 5100)
		self.assertEqual(bedLine.midpoint(True), 2600)
		self.assertEqual(bedLine.newline(10,20), 'chr1\t10\t20\tbedpe_example1\t30\t+')
		self.assertEqual(bedLine.newline(10,20, {6: "-"}), 'chr1\t10\t20\tbedpe_example1\t30\t-')
		bedLine = bedline('chr1\t100\t5103\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.midpoint(), 2601)
		bedLine = bedline('chr1\t20\t30\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.centerAndExtent(20), 'chr1\t5\t46\tbedpe_example1\t30\t+')
		
		bedLine = bedline('chr1\t20\t30\tbedpe_example1\t30\t-\n')
		self.assertEqual(bedLine.centerAndExtent(20), 'chr1\t4\t45\tbedpe_example1\t30\t-')
		self.assertEqual(bedLine.changeField(4, 'newName').getLine(), 'chr1\t20\t30\tnewName\t30\t-')

		bedLine = bedline('chr1\t2\t5\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'AAG')
		bedLine = bedline('chr1\t2\t5\tbedpe_example1\t30\t-\n')
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'CTT')
		bedLine = bedline('chr1\t0\t5\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'AAAAG')
		bedLine = bedline('chr1\t1\t5\tbedpe_example1\t30\t+\n')
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'AAAG')
		matrixString = '''
			A [1 8 4]
			T [9 8 6]
			G [0 2 0]
			C [0 0 0]'''
		
		bedLine = bedline('chr1\t2\t20\tbedpe_example1\t30\t+\n')
		newLines = bedLine.getNewLinesWithPfm('utils/testFiles/bed.fa', matrixString, True)
		self.assertEqual(newLines, [
			'chr1\t11\t14\tbedpe_example1\t30\t+',
			'chr1\t7\t10\tbedpe_example1\t30\t-',
			'chr1\t8\t11\tbedpe_example1\t30\t-',
			])
		bedLine = bedline(newLines[0])
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'TTT')
		
		matrixString = '''
			A [9 1 4 1 9 9 9]
			T [1 1 6 8 0 0 0]
			G [0 9 0 9 0 0 0]
			C [0 0 9 0 0 0 0]'''
		bedLine = bedline('chr1\t2\t20\tbedpe_example1\t30\t+\n')
		newLines = bedLine.getNewLinesWithPfm('utils/testFiles/bed.fa', matrixString, True)
		self.assertEqual(newLines, [
			'chr1\t3\t10\tbedpe_example1\t30\t+'
			])
		bedLine = bedline(newLines[0])
		self.assertEqual(bedLine.getSequence('utils/testFiles/bed.fa'), 'AGCTAAA')


	def test_bedIntersect(self):
		bedIntersectLine = 'chr1\t100\t500\t.\t.\t+\tchr1\t130\t150\t.\t.\t+\n'
		bedIntersect = bedintersect(bedIntersectLine)
		self.assertEqual(bedIntersect.length1(), 400)
		self.assertEqual(bedIntersect.length2(), 20)
		self.assertEqual(bedIntersect.position2(), 140)
		self.assertEqual(bedIntersect.relativeDistanceOfPosition2(), 40)
		self.assertEqual(bedIntersect.getDistancePercentage(), 10)
		self.assertEqual(bedIntersect.sameStrands(), True)

		bedIntersectLine = 'chr1\t100\t500\t.\t.\t-\tchr1\t130\t150\t.\t.\t+\n'
		bedIntersect = bedintersect(bedIntersectLine)
		self.assertEqual(bedIntersect.length1(), 400)
		self.assertEqual(bedIntersect.length2(), 20)
		self.assertEqual(bedIntersect.position2(), 140)
		self.assertEqual(bedIntersect.relativeDistanceOfPosition2(), 360)
		self.assertEqual(bedIntersect.getDistancePercentage(), 90)
		self.assertEqual(bedIntersect.relativeDistanceOfPosition2('start'), 350)
		self.assertEqual(bedIntersect.relativeDistanceOfPosition2('end'), 370)
		self.assertEqual(bedIntersect.getDistancePercentage(), 90)
		self.assertEqual(bedIntersect.getDistancePercentage('start', 'start'), 88)
		# self.assertEqual(bedIntersect.getDistancePercentage('end', 'start'), 93)
		self.assertEqual(bedIntersect.getDistancePercentage('start', 'end'), 13)
		self.assertEqual(bedIntersect.getDistancePercentage('end', 'end'), 8)
		self.assertEqual(bedIntersect.getDistancePercentage(), 90)

		bedIntersectLine = 'chr1\t100\t500\t.\t.\t-\tchr1\t136\t150\t.\t.\t+\n'
		bedIntersect = bedintersect(bedIntersectLine)
		self.assertEqual(bedIntersect.getDistancePercentage(), 89)
		self.assertEqual(bedIntersect.sameStrands(), False)
		self.assertEqual(bedIntersect.getAbsoluteDistance(), 357)
		self.assertEqual(bedIntersect.getAbsoluteDistance('upstream', 1), -43)
		self.assertEqual(bedIntersect.getAbsoluteDistance('upstream', 10), -5)
		self.assertEqual(bedIntersect.getAbsoluteDistance('upstream', 20), -3)
		self.assertEqual(bedIntersect.getAbsoluteDistance('downstream'), 457)
		self.assertEqual(bedIntersect.getAbsoluteDistance('downstream', 10), 136)

if __name__ == "__main__":
	unittest.main()
