from pipeTools import *

def test_codeList2multiCodeList():
    codeList = ['code.sh', '-m', 'm1', '-i', ['in1', 'in2', 'in3'], '-o', ['out1', 'out2', 'out3']]
    parallelJobLists = codeList2multiCodeList(codeList)
    results = [
        ['code.sh', '-m', 'm1', '-i', 'in1', '-o', 'out1'],
        ['code.sh', '-m', 'm1', '-i', 'in2', '-o', 'out2'],
        ['code.sh', '-m', 'm1', '-i', 'in3', '-o', 'out3']
    ]
    assert results == parallelJobLists

def test_list2gappedString():
    theList = ['a', 'b', 1, 89]
    expectedResult= 'a b 1 89'
    assert list2gappedString(theList) == expectedResult

def test_in2out():
    input = '~/dir1/dir2/somefile.otherNames.extension1'
    expectedResult= '~/dir1/dir2/somefile.otherNames.extension2'
    assert in2out(input, 'extension1', 'extension2') == expectedResult

def test_replaceLast():
    source_string = 'someNameReplacehereOtherNameReplaceherename2Replaceherename3'
    replace_what = 'Replacehere'
    replace_with = 'ThisIsReplaced'
    expectedResult = 'someNameReplacehereOtherNameReplaceherename2ThisIsReplacedname3'
    assert replaceLast(source_string, replace_what, replace_with) == expectedResult

def test_funIn2out():		
    assert funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa') == 'input.method1.getNucAbu.csv'
    assert funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 1) == 'input.method1.gNA.csv'
    assert funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 2) == 'input.method1.geNuAb.csv'
    assert funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 4) == 'input.method1.getNuclAbun.csv'
    assert funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', 'Extra') == 'input.method1.getNucAbuExtra.csv'

def test_listOperation():
    def basicSumFunction(e, *args):
        result = e
        for arg in args:
            print(arg)
            result += arg
        return result

    theList = [3,4,5,6]
    assert listOperation(basicSumFunction, theList, 1,2,3) == [9,10,11,12]