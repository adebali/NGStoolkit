import os
import unittest
from fasta import fasta

scriptDir = os.path.dirname(os.path.realpath(__file__))
kTestFilesDir = os.path.join(scriptDir, 'testFiles')
kFastaExample1 = os.path.join(kTestFilesDir, 'fastaExample1.fa')
kFastaExample2 = os.path.join(kTestFilesDir, 'fastaExample2.fa')
kFastaExample4 = os.path.join(kTestFilesDir, 'fastaExample4.fa')

fasta1 = fasta(kFastaExample1)
fasta2 = fasta(kFastaExample2)
fasta4 = fasta(kFastaExample4)

class fastaTests(unittest.TestCase):
	def test_sequenceCount(self):
		self.assertEqual(fasta1.getSequenceCount(), 2)
		self.assertEqual(fasta2.getSequenceCount(), 4)

	def test_firstSeqLength(self):
		self.assertEqual(fasta1.getFirstSeqLength(), 7)
		self.assertEqual(fasta2.getFirstSeqLength(), 11)

	def test_maxSeqLength(self):
		self.assertEqual(fasta1.getMaxSeqLength(), 7)
		self.assertEqual(fasta2.getMaxSeqLength(), 25)

	def test_nucleotideAbundance(self):
		expectedResult = {
			1: { 'A': 1, 'T':0, 'G': 1, 'C': 0},
			2: { 'A': 0, 'T':1, 'G': 0, 'C': 0},
			3: { 'A': 0, 'T':0, 'G': 2, 'C': 0},
			4: { 'A': 0, 'T':1, 'G': 0, 'C': 1},
			5: { 'A': 0, 'T':1, 'G': 1, 'C': 0},
			6: { 'A': 0, 'T':1, 'G': 1, 'C': 0},
			7: { 'A': 2, 'T':0, 'G': 0, 'C': 0}
		}
		self.assertEqual(fasta1.getNucleotideAbundance(), expectedResult)
		self.assertEqual(fasta1.getNucleotideAbundance(7), expectedResult)
		self.assertEqual(fasta1.getKmerAbundance(1), expectedResult)
		del expectedResult[7]
		self.assertEqual(fasta1.getNucleotideAbundance(6), expectedResult)


	def test_getNucleotidePercentages(self):
		expectedResult = {
			1: { 'A': 0.5, 'T':0, 'G': 0.50, 'C': 0},
			2: { 'A': 0, 'T':1, 'G': 0, 'C': 0},
			3: { 'A': 0, 'T':0, 'G': 1, 'C': 0},
			4: { 'A': 0, 'T':0.5, 'G': 0, 'C': 0.5},
			5: { 'A': 0, 'T':0.5, 'G': 0.5, 'C': 0},
			6: { 'A': 0, 'T':0.5, 'G': 0.5, 'C': 0},
			7: { 'A': 1, 'T':0, 'G': 0, 'C': 0}
		}
		self.assertEqual(fasta1.getNucleotidePercentages(fasta1.getNucleotideAbundance()), expectedResult)


	def test_getKmerAbundance(self):
		expectedResult = {
			1: {'AA': 0, 'AT': 1, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 0, 'TG': 0, 'TC': 0,\
				'GA': 0, 'GT': 0, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0, },
			2: {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 0, 'TG': 1, 'TC': 0,\
				'GA': 0, 'GT': 0, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0, },
			3: {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 0, 'TG': 0, 'TC': 0,\
				'GA': 0, 'GT': 1, 'GG': 0, 'GC': 1, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0, },
			4: {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 1, 'TG': 0, 'TC': 0,\
				'GA': 0, 'GT': 0, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 1, 'CC': 0, },
			5: {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 0, 'TG': 1, 'TC': 0,\
				'GA': 0, 'GT': 1, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0, },
			6: {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 1, 'TT': 0, 'TG': 0, 'TC': 0,\
				'GA': 1, 'GT': 0, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0, },
		}

		self.assertEqual(fasta1.getKmerAbundance(2), expectedResult)

	def test_FastaStream(self):
		expectedObjects_1 = [{'h': 'header1', 's': 'ATGCGtA'}, {'h': 'header2', 's':'GXGTTGA'}]
		fastaStream_1 = fasta1.stream(4)
		results_1 = []
		for e in fastaStream_1:
			results_1.append(e)
		self.assertEqual(results_1, expectedObjects_1)
		

		fastaStream_4 = fasta4.stream(6)
		expectedObjects_4 = [{'h': 'header1', 's': 'ATGCGtA'}, {'h': 'header2', 's':'GXGTTGAATGCGTA'}, {'h': 'header3', 's':'headerthree'}]
		results_4 = []
		for e in fastaStream_4:
			results_4.append(e)
		self.assertEqual(results_4, expectedObjects_4)

		fastaStream_4 = fasta4.stream(4096)
		results_4b = []
		for e in fastaStream_4:
			results_4b.append(e)
		self.assertEqual(results_4, results_4b)


if __name__ == "__main__":
	unittest.main()
