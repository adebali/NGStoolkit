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

	def writeSingleLineFastq(self, outputFile):
		out = open(outputFile, 'w')
		for seqObject in self.stream():
			out.write(seqObject.getSingleLine())
		out.close()
		# count = 0
		# filein = open(self.file, 'r')
		# newLine = ''
		# for line in filein:
		# 	count += 1
		# 	newLine += line.strip() + '\t'				
		# 	if count%4 == 0:
		# 		if newLine != '':
		# 			out.write(newLine.strip() + '\n')
		# 		newLine = ''

	# def stream(self, bufsize=40960):
	# 	def chunk2object(chunk):
	# 		lines = chunk.split('\n')
	# 		seqObject = fastqSeq()
	# 		seqObject.assignFourLines(lines)
	# 		return seqObject

	# 	filein = open(self.file, 'r')
	# 	delimiter = '\n@'
	# 	buf = ''
	# 	justStarted = True
	# 	while True:
	# 		newbuf = filein.read(bufsize)
	# 		if not newbuf:
	# 			yield chunk2object(buf)
	# 			return
	# 		buf += newbuf
	# 		sequenceChunks = buf.split(delimiter)
	# 		for chunk in sequenceChunks[0:-1]:
	# 			if justStarted and chunk.startswith('@'):
	# 				chunk = chunk[1:]
	# 				justStarted = False
	# 			yield chunk2object(chunk)
	# 		buf = sequenceChunks[-1]

	def stream(self, bufsize=40960):
		def fourLines2object(lines):
			seqObject = fastqSeq()
			seqObject.assignFourLines(lines)
			return seqObject

		filein = open(self.file, 'r')
		delimiter = '\n@'
		buf = ''
		justStarted = True
		i = 0
		fourLines = []
		for line in filein:
			i += 1
			fourLines.append(line)
			if i % 4 == 0:
				yield fourLines2object(fourLines)
				fourLines = []


class fastqSeq:
	def __init__(self, sequence = ''):
		self.sequence = sequence

	def assignSequence(self, sequence):
		self.sequence = sequence

	def assignHeader(self, header):
		self.header = header

	def assignQuality(self, quality):
		self.quality = quality

	def assignDescription(self, description):
		self.description = description

	def assignFourLines(self, lines):
		'''lines is a list of fastq lines'''
		header = lines[0].strip()
		sequence = lines[1].strip()
		description = lines[2].strip()
		quality = lines[3].strip()
		self.assignSequence(sequence)
		self.assignHeader(header)
		self.assignDescription(description)
		self.assignQuality(quality)

	def getAsString(self):
		return('\n'.join([self.header, self.sequence, self.description, self.quality]) + '\n')

	def getSingleLine(self):
		return('\t'.join([self.header, self.sequence, self.description, self.quality]) + '\n')
		