from sequence import *

def test_reverseComplement():
    sequence = DNA('ATCGGCAntAxXGCG', False)
    assert sequence.reverseComplement().getSequence() == 'CGCXxTanTGCCGAT'
    sequence = DNA('ATCGGCAntAxXGCG')
    print (sequence.getSequence())
    assert sequence.reverseComplement().getSequence() == 'CGCXXTANTGCCGAT'

def test_subseq():
    sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
    assert sequence.subseq(2, 5).getSequence() == 'TTT'

def test_motifCount():
    sequence = DNA('AATTTAGCGTTAGCTGCTTTT')	
    assert sequence.motifCount('TT') == 6

def test_motifClass():
    motif = reMotif('(A|T)CGC')
    complementMotif = reMotif('(T|A)GCG')
    reverseMotif = reMotif('CGC(T|A)')
    reverseComplementMotif = reMotif('GCG(A|T)')
    assert motif.DNA_complement().reverse().getPatternAsString() == reverseComplementMotif.getPatternAsString()
    assert motif.reverse().getPatternAsString() == reverseMotif.getPatternAsString()
    assert motif.DNA_complement().getPatternAsString() == complementMotif.getPatternAsString()
    assert motif.getIndexList('acgtACgCaaaACGCaaaa') == ([4, 11], [8,15])

def test_consensus():
    consensusString = 'NBNBATTTCCSGGAARTSNNN'
    con = consensus(consensusString)
    assert con.getPatternAsString() == '.(C|G|T).(C|G|T)(A)(T)(T)(T)(C)(C)(C|G)(G)(G)(A)(A)(A|G)(T)(C|G)...'

def test_DNA():
    dna = DNA("ATGC")
    assert dna.reverseComplement().getSequence() == "GCAT"
