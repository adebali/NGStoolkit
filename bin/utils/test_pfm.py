from pfm import *

def test_pfm():

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
    assert pfm.pfm == expectedDictionary
    assert pfm.length == 10
    assert pfm.sumPosition(3) == 7639
    assert pfm.getRatio('A', 3) == float(3)/7639

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
    # assert pfm.getTopKeyFromDictionary({}), None
    assert pfm.getSequenceScore('AAT') == float(1)/20 * 2/70 * 31/102
    assert pfm.transpose() == transposedShortDict
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGCAAA' + 'GGGGGGGGGG'
    pfm = Pfm()
    pfm.takeStringInput(egPfmString)
    assert pfm.getSequenceScore('TTTAAATTTG') == 0
    assert pfm.getBestHits(seq) == [10, 30]
    assert pfm.getHitsFromSequence(seq) == [10, 30]
    assert pfm.getIndex(seq) == [10, 30]
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
    assert pfm.getBestHits(seq) == [10]
    assert pfm.getMostLikelySequence(), 'AAACCGCAAA'
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
    assert pfm.getHitsFromSequence(seq) == [10]
    assert pfm.getHits(seq) == {0: [10], 1: []}
    seq = 'TTTAAATTTG' + 'AAAAAAAAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
    assert pfm.getHitsFromSequence(seq) == [30]
    seq = 'CCCCCCCCCC' + 'TTTCCGGTTT' + 'CAAATTTAAA' + 'TTTTTTTTTT' + 'CAAATTTAAA'
    seq_rc = DNA(seq).reverseComplement().getSequence()
    assert pfm.getHitsFromSequence(seq_rc) == [30]
    seq = 'TTTAAATTTG' + 'AAAAAAAAAA' + 'TTTAAATTTG' + 'AAACCGGAAA' + 'GGGGGGGGGG'
    assert pfm.reverseComplementIndex(pfm.getHitsFromSequence(seq_rc), seq_rc) == [10]
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [10], 1: [30]}
    seq = 'TTTAAATTTG' + 'AAACCACAAA' + 'TTTAAATTTG' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [], 1: [30]}
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTAAATTTG' + 'TTTGCGGTTA' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [10], 1: []}
    seq = 'TTTAAATTTG' + 'AAACCGCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [10], 1: [20, 30]}
    seq = 'TTTAAATTTG' + 'AAACNgCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [10], 1: []}
    seq = 'TTTAAATTTG' + 'AAACNgCAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [10], 1: []}
    seq = 'TTTAAATTTG' + 'AAACNgXAAA' + 'TTTGCGGTTT' + 'TTTGCGGTTT' + 'GGGGGGGGGG'
    assert pfm.getHits(seq) == {0: [], 1: [20, 30]}