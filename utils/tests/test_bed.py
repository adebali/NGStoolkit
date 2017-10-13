from bed import *


theTestBedFile = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testFiles', 'bed.fa') 

def test_bedpeLine2bedLine():
    bedpeLine = 'chr1\t100\t200\tchr1\t5000\t5100\tbedpe_example1\t30\t+\t-\n'
    expectedResult1 = 'chr1\t100\t5100\tbedpe_example1\t30\t+'
    assert bedpeLine2bedLine(bedpeLine) == expectedResult1
    bedpeLineMinus = 'chr1\t1000\t1200\tchr1\t900\t950\tbedpe_example2\t30\t-\t+\n'
    expectedResult2 = 'chr1\t900\t1200\tbedpe_example2\t30\t-'
    assert bedpeLine2bedLine(bedpeLineMinus) == expectedResult2
    bedpeLineCustomFields = bedpeLine.strip() + '\tField1\tField2\n'
    expectedResult3 = expectedResult1 + '\tField1\tField2'
    assert bedpeLine2bedLine(bedpeLineCustomFields) == expectedResult3


def test_bedLine2fixedRangeLine():
    bedLine = 'chr1\t100\t5100\tbedpe_example1\t30\t+\n'
    expectedResult = 'chr1\t100\t120\tbedpe_example1\t30\t+'
    assert bedLine2fixedRangedLine(bedLine, 'l', 20) == expectedResult
    expectedResult = 'chr1\t5075\t5100\tbedpe_example1\t30\t+'
    assert bedLine2fixedRangedLine(bedLine, 'r', 25) == expectedResult

    bedLine2 = 'chr1\t100000212\t100000296\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t+'
    expectedResult2 = 'chr1\t100000212\t100000222\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t+'
    assert bedLine2fixedRangedLine(bedLine2, 'l', 10) == expectedResult2

    bedLine3 = 'chr1\t100\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
    expectedResult3 = 'chr1\t140\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
    assert bedLine2fixedRangedLine(bedLine3, 'l', 10) == expectedResult3

    bedLine4 = 'chr1\t100\t150\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
    expectedResult4 = 'chr1\t100\t110\tHISEQ:355:C9U5RANXX:4:1310:14307:58329\t255\t-'
    assert bedLine2fixedRangedLine(bedLine3, 'r', 10) == expectedResult4

# def test_bedClass():
# 	myBed = bed(os.path.join(os.path.dirname(os.path.realpath(__file__)),'testFiles/bedExample.bed'))
# 	myBed.fixRange('l', 10)

# 	myBed = bed(os.path.join(os.path.dirname(os.path.realpath(__file__)),'testFiles/genes.bed'))
# 	myBed.removeNeighbors(20)

def test_bedline2():
    bedLine = bedline('chr1\t100\t5100\tbedpe_example1\t30\t+\n')
    assert bedLine.fields() == ['chr1','100','5100','bedpe_example1','30','+']
    assert bedLine.chromosome() == 'chr1'
    assert bedLine.start() == 100
    assert bedLine.end() == 5100
    assert bedLine.midpoint(True) == 2600
    assert bedLine.newline(10,20) == 'chr1\t10\t20\tbedpe_example1\t30\t+'
    assert bedLine.newline(10,20, {6: "-"}) == 'chr1\t10\t20\tbedpe_example1\t30\t-'
    bedLine = bedline('chr1\t100\t5103\tbedpe_example1\t30\t+\n')
    assert bedLine.midpoint() == 2601
    bedLine = bedline('chr1\t20\t30\tbedpe_example1\t30\t+\n')
    assert bedLine.centerAndExtent(20) == 'chr1\t5\t46\tbedpe_example1\t30\t+'
    
    bedLine = bedline('chr1\t20\t30\tbedpe_example1\t30\t-\n')
    assert bedLine.centerAndExtent(20) == 'chr1\t4\t45\tbedpe_example1\t30\t-'
    assert bedLine.changeField(4, 'newName').getLine() == 'chr1\t20\t30\tnewName\t30\t-'

    bedLine = bedline('chr1\t2\t5\tbedpe_example1\t30\t+\n')
    assert bedLine.getSequence(theTestBedFile) == 'AAG'
    bedLine = bedline('chr1\t2\t5\tbedpe_example1\t30\t-\n')
    assert bedLine.getSequence(theTestBedFile) == 'CTT'
    bedLine = bedline('chr1\t0\t5\tbedpe_example1\t30\t+\n')
    assert bedLine.getSequence(theTestBedFile) == 'AAAAG'
    bedLine = bedline('chr1\t1\t5\tbedpe_example1\t30\t+\n')
    assert bedLine.getSequence(theTestBedFile) == 'AAAG'
    matrixString = '''
        A [1 8 4]
        T [9 8 6]
        G [0 2 0]
        C [0 0 0]'''
    
    bedLine = bedline('chr1\t2\t20\tbedpe_example1\t30\t+\n')
    newLines = bedLine.getNewLinesWithPfm(theTestBedFile, matrixString, True)
    assert newLines == [ 
        'chr1\t11\t14\tbedpe_example1\t30\t+',
        'chr1\t7\t10\tbedpe_example1\t30\t-',
        'chr1\t8\t11\tbedpe_example1\t30\t-',
        ]
    bedLine = bedline(newLines[0])
    assert bedLine.getSequence(theTestBedFile) == 'TTT'
    
    matrixString = '''
        A [9 1 4 1 9 9 9]
        T [1 1 6 8 0 0 0]
        G [0 9 0 9 0 0 0]
        C [0 0 9 0 0 0 0]'''
    bedLine = bedline('chr1\t2\t20\tbedpe_example1\t30\t+\n')
    newLines = bedLine.getNewLinesWithPfm(theTestBedFile, matrixString, True)
    assert newLines == [ 
        'chr1\t3\t10\tbedpe_example1\t30\t+'
        ]
    bedLine = bedline(newLines[0])
    assert bedLine.getSequence(theTestBedFile) == 'AGCTAAA'


def test_bedIntersect():
    bedIntersectLine = 'chr1\t100\t500\t.\t.\t+\tchr1\t130\t150\t.\t.\t+\n'
    bedIntersect = bedintersect(bedIntersectLine)
    assert bedIntersect.length1() == 400
    assert bedIntersect.length2() == 20
    assert bedIntersect.position2() == 140
    assert bedIntersect.relativeDistanceOfPosition2() == 40
    assert bedIntersect.getDistancePercentage() == 10
    assert bedIntersect.sameStrands() == True

    bedIntersectLine = 'chr1\t100\t500\t.\t.\t-\tchr1\t130\t150\t.\t.\t+\n'
    bedIntersect = bedintersect(bedIntersectLine)
    assert bedIntersect.length1() == 400
    assert bedIntersect.length2() == 20
    assert bedIntersect.position2() == 140
    assert bedIntersect.relativeDistanceOfPosition2() == 360
    assert bedIntersect.getDistancePercentage() == 90
    assert bedIntersect.relativeDistanceOfPosition2('start') == 350
    assert bedIntersect.relativeDistanceOfPosition2('end') == 370
    assert bedIntersect.getDistancePercentage() == 90
    assert bedIntersect.getDistancePercentage('start', 'start') == 88
    # assert bedIntersect.getDistancePercentage('end', 'start') == 93
    assert bedIntersect.getDistancePercentage('start', 'end') == 13
    assert bedIntersect.getDistancePercentage('end', 'end') == 8
    assert bedIntersect.getDistancePercentage() == 90

    bedIntersectLine = 'chr1\t100\t500\t.\t.\t-\tchr1\t136\t150\t.\t.\t+\n'
    bedIntersect = bedintersect(bedIntersectLine)
    assert bedIntersect.getDistancePercentage() == 89
    assert bedIntersect.sameStrands() == False
    assert bedIntersect.getAbsoluteDistance() == 357
    assert bedIntersect.getAbsoluteDistance('upstream', 1) == -43
    assert bedIntersect.getAbsoluteDistance('upstream', 10) == -5
    assert bedIntersect.getAbsoluteDistance('upstream', 20) == -3
    assert bedIntersect.getAbsoluteDistance('downstream') == 457
    assert bedIntersect.getAbsoluteDistance('downstream', 10) == 136
