import os
import re
import itertools
import unittest

kDNAcomplementaryDictionary = {
	'A':'T', 'T':'A', 'G':'C', 'C':'G',
	'a':'t', 't':'a', 'g':'c', 'c':'g',
	'n':'n', 'N':'N', 'x':'x', 'X':'X',
	'B':'V', 'V':'B', 'D':'H', 'H':'D',
	'b':'v', 'v':'b', 'd':'h', 'h':'d'
}

kDNAmotifDictionary = {
	'A': '(A)',
	'C': '(C)',
	'G': '(G)',
	'T': '(T)',
	'U': '(U)',
	'W': '(A|T)',
	'S': '(C|G)',
	'M': '(A|C)',
	'K': '(G|T)',
	'R': '(A|G)',
	'Y': '(C|T)',
	'B': '(C|G|T)',
	'D': '(A|G|T)',
	'H': '(A|C|T)',
	'V': '(A|C|G)',
	'N': '.'
	}

class Sequence:
	def __init__(self, input, up=True):
		self.sequence = input
		if up:
			self.sequence = self.sequence.upper()
		self.header = None

	def subseq(self, start, end):
		# Zero-based subsequence function
		return DNA(self.sequence[start : end])

	def motifCount(self, motif):
		count = start = 0
		while True:
			start = self.sequence.find(motif, start) + 1
			if start > 0:
				count+=1
			else:
				return count

	def assignHeader(self, header):
		self.header = header

	def getHeader(self):
		return self.header

	def getLength(self):
		return len(self.getSequence())

	def getSequence(self):
		return self.sequence


class DNA(Sequence):
	def reverseComplement(self):
		seq_dict = kDNAcomplementaryDictionary
		return DNA("".join([seq_dict[base] for base in reversed(self.sequence)]))

class reMotif(object):
	def __init__(self, input, ignorecaseFlag=True):
		self.patternAsString = input
		if ignorecaseFlag:
			self.pattern = re.compile(self.patternAsString, re.IGNORECASE)
		else:
			self.pattern = re.compile(self.patternAsString)

	def getPatternAsString(self):
		return self.patternAsString

	def DNA_complement(self, reverseFlag = False):
		def complementIfPossibe(key, dictionary):
			if key in dictionary.keys():
				return dictionary[key]
			return key
		if reverseFlag: # Use with caution. It doesn't work for non-symmetric patterns like: .{3}ATG.{4} will give .{3}GTA.{4}
			string = reverse(self.getPatternAsString())
		else:
			string = self.getPatternAsString()
		return reMotif("".join([complementIfPossibe(base, kDNAcomplementaryDictionary) for base in string]))

	def reverse(self):
		replacements = {"[": "]", "]": "[", "(": ")", ")": "(", "{": "}", "}": "{"}
		reversedString = reversed(self.getPatternAsString())
		replacedString = "".join(replacements.get(c, c) for c in reversedString)
		return reMotif(replacedString)

	def getIndexList(self, subject):
		queryPattern = self.pattern
		matchObjects = re.finditer(queryPattern, subject)
		startList = []
		endList = []
		for matchObject in matchObjects:
			start = matchObject.start()
			end = matchObject.end()
			startList.append(start)
			endList.append(end)
		return startList, endList

class consensus(reMotif):
	def __init__(self, input):
		convertedInput = self.toRegex(input)
		reMotif.__init__(self, convertedInput)

	def toRegex(self, string):
		newString = ''
		for letter in string:
			upperLetter = letter.upper()
			if upperLetter in kDNAmotifDictionary.keys():
				newString += kDNAmotifDictionary[upperLetter]
			else:
				newString += upperLetter
		return newString

class tests(unittest.TestCase):
	def test_reverseComplement(self):
		sequence = DNA('ATCGGCAntAxXGCG')
		self.assertEqual(sequence.reverseComplement().getSequence(), 'CGCXxTanTGCCGAT')

	def test_subseq(self):
		sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
		self.assertEqual(sequence.subseq(2, 5).getSequence(), 'TTT')

	def test_motifCount(self):
		sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
		self.assertEqual(sequence.motifCount('TT'), 6)

	def test_motifClass(self):
		motif = reMotif('(A|T)CGC')
		complementMotif = reMotif('(T|A)GCG')
		reverseMotif = reMotif('CGC(T|A)')
		reverseComplementMotif = reMotif('GCG(A|T)')
		self.assertEqual(motif.DNA_complement().reverse().getPatternAsString(), reverseComplementMotif.getPatternAsString())
		self.assertEqual(motif.reverse().getPatternAsString(), reverseMotif.getPatternAsString())
		self.assertEqual(motif.DNA_complement().getPatternAsString(), complementMotif.getPatternAsString())
		self.assertEqual(motif.getIndexList('acgtACgCaaaACGCaaaa'), ([4, 11], [8,15]))

	def test_consensus(self):
		consensusString = 'NBNBATTTCCSGGAARTSNNN'
		con = consensus(consensusString)
		self.assertEqual(con.getPatternAsString(), '.(C|G|T).(C|G|T)(A)(T)(T)(T)(C)(C)(C|G)(G)(G)(A)(A)(A|G)(T)(C|G)...')

	def test_DNA(self):
		dna = DNA("ATGC")
		self.assertEqual(dna.reverseComplement().getSequence(), "GCAT")

if __name__ == "__main__":
	unittest.main()

