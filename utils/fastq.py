class fastq:
	def __init__(self, input):
		self.file = input

	def getSequenceCount(self):
		count = 0
		filein = open(self.file, 'r')
		for line in filein:
			count += 1
		if count%4 != 0:
			sys.exit('incorrect fastq format, exiting...')
		else:
			return count/4

	def getMaxSeqLength(self):
		maxSeqLength = 0
		with open(self.file) as f:
			for line in islice(f, 1, None, 4):
				seqLength = len(line.strip())
				maxSeqLength = max(maxSeqLength, seqLength)
		return maxSeqLength

	def getLengthDistribution(self, defaultMaxLength=50):
		from itertools import islice
		theDict = {}

		if defaultMaxLength == 0:
			maxSeqLength = getMaxSeqLength(self.file)
		else:
			maxSeqLength = defaultMaxLength

		for i in range(0, maxSeqLength + 1):
			theDict[i] = 0

		with open(self.file) as f:
			for line in islice(f, 1, None, 4):
				seqLength = len(line.strip())
				theDict[seqLength] += 1

		return theDict
