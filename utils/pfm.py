import os
import re
import itertools
import unittest
from sequence import Sequence, DNA, reMotif

class Pfm:
    def __init__(self, dictionary = {}):
        self.pfm = dictionary
        self.assignLength()

    def assignLength(self):
        length = 0
        if self.pfm.keys() != []:
            firstKey = self.pfm.keys()[0]
            length = len(self.pfm[firstKey])
        self.length = length

    def values2integers(self, dictionary):
        newDictionary = {}
        for key in dictionary.keys():
            newDictionary[key] = []
            for e in dictionary[key]:
                newDictionary[key].append(int(e.strip()))
        return newDictionary

    def takeStringInput(self, string):
        tempDictionary = {}
        for line in string.strip().split("\n"):
            line = re.sub(' +',' ',line).strip().replace('[', '').replace(']', '').strip()
            if line:
                ll = line.split(' ')
                tempDictionary[ll[0]] = ' '.join(ll[1:]).strip().split(' ')
        self.pfm = self.values2integers(tempDictionary)
        self.assignLength()

    def getMostLikelySequence(self):
        transposedD = self.transpose()
        sequence = ''
        for i in range(self.length):
            sequence += self.getTopKeyFromDictionary(transposedD[i])[0]
        return sequence

    def sumPosition(self, position):
        total = 0
        for key in self.pfm.keys():
            total += self.pfm[key][position]
        return total

    def getRatio(self, letter, position):
        if letter == "N":
            return 1
        if not letter in self.pfm.keys():
            return 0
        return float(self.pfm[letter][position])/self.sumPosition(position)

    def getSequenceScore(self, seq):
        score = 1
        seq = seq.upper()
        for i in range(len(seq)):
            letter = seq[i]
            score *= self.getRatio(letter, i)
        return score

    def invertDictionary(self, dictionary):
        inv_map = {}
        for k, v in dictionary.iteritems():
            inv_map[v] = inv_map.get(v, [])
            inv_map[v].append(k)
        return inv_map

    def getBestHits(self, seq, strand=False):
        scoreDict = self.getScoreDictionary(seq)
        return self.getTopKeyFromDictionaryAndSeq(scoreDict, seq)

    def getScoreDictionary(self, seq):
        seqObject = Sequence(seq)
        scoreDict = {}
        for i in range(0, seqObject.getLength() - self.length):
            subseq = seqObject.subseq(i, i + self.length)
            scoreDict[i] = self.getSequenceScore(subseq.getSequence())
        return scoreDict

    def reverseComplementIndex(self, hits, seq):
        seqLength = len(seq)
        # hits = self.getHitsFromSequence(seq)
        newList = []
        for hit in hits:
            newHit = seqLength - (hit + self.length)
            newList.append(newHit)
        return newList

    def getHitsFromSequence(self, seq, exactMatch = False):
        identicalMatches = self.getIndex(seq)
        if identicalMatches == []:
            if exactMatch == False:
                return self.getBestHits(seq)
            else:
                return identicalMatches
        else:
            return identicalMatches

    def getHits(self, seq):
        d0 = self.getScoreDictionary(seq)
        seq_rc = DNA(seq).reverseComplement().getSequence()
        d1 = self.getScoreDictionary(seq_rc)
        return self.getTopKeyDictionary(d0, d1, seq)

    def getTopKeyFromDictionary(self, d):
        if d != {}:
            invertedD = self.invertDictionary(d)
            topValue = list(reversed(sorted(invertedD.keys())))[0]
            result = invertedD[topValue]
            return result
        return None


    def getTopKeyFromDictionaryAndSeq(self, d, seq, reverseFlag = False):
        if d != {}:
            invertedD = self.invertDictionary(d)
            topValue = list(reversed(sorted(invertedD.keys())))[0]
            result = invertedD[topValue]
            if reverseFlag == True:
                result = self.reverseComplementIndex(result, seq)
            return result
        return None

    def getTopKeyDictionary(self, sortDict0, sortDict1, seq):
        invertedD0 = self.invertDictionary(sortDict0)
        invertedD1 = self.invertDictionary(sortDict1)
        topValue0 = list(reversed(sorted(invertedD0.keys())))[0]
        topValue1 = list(reversed(sorted(invertedD1.keys())))[0]
        hits0 = sorted(invertedD0[topValue0])
        hits1 = sorted(self.reverseComplementIndex(invertedD1[topValue1], seq))
        if topValue0 > topValue1:
            return {0: hits0, 1: []}
        elif topValue1 > topValue0:
            return {0: [], 1: hits1}
        else:
            return {0: invertedD0[topValue0], 1: hits1}
            
            

    def transpose(self):
        newDict = {}
        for i in range(0, self.length):
            newDict[i] = {}
            for key in self.pfm.keys():
                newDict[i][key] = self.pfm[key][i]
        return newDict

    def getIndex(self, sequence):
        def find_offsets(haystack, needle):
            """
            Find the start of all (possibly-overlapping) instances of needle in haystack
            """
            offs = -1
            while True:
                offs = haystack.find(needle, offs+1)
                if offs == -1:
                    break
                else:
                    yield offs
        seq = self.getMostLikelySequence()
        return list(find_offsets(sequence, seq))


class tests(unittest.TestCase):

    def test_pfm(self):

        egPfmString = '''
    A  [3952 7636 7636    3    0 2991   30 7636 7636 4131 ]
    C  [1515  692  331 7636 7636   76 7636  368 1015 1198 ]
    G  [ 160 2431  273    0   22 4645   96  202 1772  861 ]
    T  [2008  438    0    0    0   41   52  184  964 1445 ]'''

        expectedDictionary = {
            'A': [3952, 7636, 7636, 3, 0, 2991, 30, 7636, 7636, 4131],
            'C': [1515, 692, 331, 7636, 7636, 76, 7636, 368, 1015, 1198],
            'G': [160, 2431, 273, 0, 22, 4645, 96, 202, 1772, 861],
            'T': [2008, 438, 0, 0, 0, 41, 52, 184, 964, 1445]
        }
        pfm = Pfm()
        pfm.takeStringInput(egPfmString)
        self.assertEqual(pfm.pfm, expectedDictionary)
        self.assertEqual(pfm.length, 10)
        self.assertEqual(pfm.sumPosition(3), 7639)
        self.assertEqual(pfm.getRatio('A', 3), float(3)/7639)

        shortDict = {
            'A': [1, 2, 3],
            'T': [5, 12, 31],
            'G': [6, 24, 33],
            'C': [8, 32, 35]
        }
        transposedShortDict = {
            0: {'A': 1, 'T': 5, 'G': 6, 'C':8},
            1: {'A': 2, 'T': 12, 'G': 24, 'C':32},
            2: {'A': 3, 'T': 31, 'G': 33, 'C':35},
        }
        pfm = Pfm(shortDict)
        # self.assertEqual(pfm.getTopKeyFromDictionary({}), None)
        self.assertEqual(pfm.getSequenceScore('AAT'), float(1)/20 * 2/70 * 31/102)
        self.assertEqual(pfm.transpose(), transposedShortDict)
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGCAAA' + 'GGGGGGGGGG'
        pfm = Pfm()
        pfm.takeStringInput(egPfmString)
        self.assertEqual(pfm.getSequenceScore('TTTAAATTTG'), 0)
        self.assertEqual(pfm.getBestHits(seq), [10, 30])
        self.assertEqual(pfm.getHitsFromSequence(seq), [10, 30])
        self.assertEqual(pfm.getIndex(seq), [10, 30])
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getBestHits(seq), [10])
        self.assertEqual(pfm.getMostLikelySequence(), 'AAACCGCAAA')
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHitsFromSequence(seq), [10])
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: []})
        seq = 'TTTAAATTTG' + 'AAAAAAAAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHitsFromSequence(seq), [30])
        seq = 'CCCCCCCCCC' + 'TTTCCGGTTT' + 'CAAATTTAAA' + 'TTTTTTTTTT' + 'CAAATTTAAA'
        seq_rc = DNA(seq).reverseComplement().getSequence()
        self.assertEqual(pfm.getHitsFromSequence(seq_rc), [30])
        seq = 'TTTAAATTTG' + 'AAAAAAAAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
        self.assertEqual(pfm.reverseComplementIndex(pfm.getHitsFromSequence(seq_rc), seq_rc), [10])
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: [30]})
        seq = 'TTTAAATTTG' + 'AAACCACAAA' + 'TTTAAATTTG' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [], 1: [30]})
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'TTTGCGGTTA' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: []})
        seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: [20, 30]})
        seq = 'TTTAAATTTG' + 'AAACNgCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: []})
        seq = 'TTTAAATTTG' + 'AAACNgCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [10], 1: []})
        seq = 'TTTAAATTTG' + 'AAACNgXAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
        self.assertEqual(pfm.getHits(seq), {0: [], 1: [20, 30]})
if __name__ == "__main__":
    unittest.main()

