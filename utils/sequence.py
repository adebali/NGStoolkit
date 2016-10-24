import os
import re
import itertools
import unittest

class Sequence:
	def __init__(self, input):
		self.sequence = input

	def subseq(self, start, end):
		return DNA(self.sequence[start - 1 : end - 1])

	def motifCount(self, motif):
		count = start = 0
		while True:
			start = self.sequence.find(motif, start) + 1
			if start > 0:
				count+=1
			else:
				return count

	def getSequence(self):
		return self.sequence


class DNA(Sequence):
	def reverseComplement(self):
		seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
					'a':'t', 't':'a', 'g':'c', 'c':'g',
					'n':'n', 'N':'N', 'x':'x', 'X':'X'}
		return DNA("".join([seq_dict[base] for base in reversed(self.sequence)]))

class tests(unittest.TestCase):
	def test_reverseComplement(self):
		sequence = DNA('ATCGGCAntAxXGCG')
		self.assertEqual(sequence.reverseComplement().getSequence(), 'CGCXxTanTGCCGAT')

	def test_subseq(self):
		sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
		self.assertEqual(sequence.subseq(2, 5).getSequence(), 'ATTT')

	def test_motifCount(self):
		sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
		self. assertEqual(sequence.motifCount('TT'), 6)

if __name__ == "__main__":
	unittest.main()

