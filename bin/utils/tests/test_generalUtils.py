import os
import generalUtils


def test_list2gappedString():
	theList = ['a', 'b', 1, 89]
	expectedResult= 'a b 1 89'
	assert generalUtils.list2gappedString(theList) == expectedResult

def test_in2out():
	input = '~/dir1/dir2/somefile.otherNames.extension1'
	expectedResult= '~/dir1/dir2/somefile.otherNames.extension2'
	assert generalUtils.in2out(input, 'extension1', 'extension2') == expectedResult

def test_replaceLast():
	source_string = 'someNameReplacehereOtherNameReplaceherename2Replaceherename3'
	replace_what = 'Replacehere'
	replace_with = 'ThisIsReplaced'
	expectedResult = 'someNameReplacehereOtherNameReplaceherename2ThisIsReplacedname3'
	assert generalUtils.replaceLast(source_string, replace_what, replace_with) == expectedResult

def test_linesAreDuplicate():
	line1 = 'chr1\t15\t30\tchr1\t45\t65\tname1\t+\n'
	line2 = 'chr1\t15\t30\tchr1\t45\t65\tname2\t+\n'
	assert generalUtils.linesAreDuplicate(line1, line2, '1,2,3,4,5,6,8') == True
	assert generalUtils.linesAreDuplicate(line1, line2, '1,2,3,4,5,6,7') == False

def test_line2switchColumnsBasedOnValue():
	line = 'chr\t100\t135\t+\t30\t40'
	expectedResult = 'chr\t100\t135\t+\t40\t30'
	valueColNo = 3
	value = '+'
	col1no = 4
	col2no = 5
	assert generalUtils.line2switchColumnsBasedOnValue(line, valueColNo, value, col1no, col2no, separator = '\t') == expectedResult

def test_reverseComplement():
	assert generalUtils.dna2reverseComplement('ATCGGCAntAxXGCG') == 'CGCXxTanTGCCGAT'

def test_subseq():
	assert generalUtils.subseq('ATGCGCA', 2, 5) == 'TGCG'

def test_motifCount():
	assert generalUtils.motifCount('AATTTAGCGTTAGCTGCTTTT', 'TT') == 6


def test_table2dictionary():
	scriptDir = os.path.dirname(os.path.realpath(__file__))
	kTestFilesDir = os.path.join(scriptDir, 'testFiles')
	table1 = os.path.join(kTestFilesDir, 'table.csv')
	dict1 = {'key1': [{'key':'key1','value': '1'}], 'key2': [{'key':'key2','value': '2'}]}
	table2 = os.path.join(kTestFilesDir, 'table.2.csv')
	dict2 = {'A': [{'key':'key1','value': '1', 'group': 'A'}, {'key':'key3','value': '3', 'group': 'A'}], 'B': [{'key':'key2','value': '2', 'group': 'B'},{'key':'key4','value': '4', 'group': 'B'}]}
	assert generalUtils.table2dictionary(table1, 'key', ',') == dict1
	assert generalUtils.table2dictionary(table2, 'group', ',') == dict2
	assert generalUtils.dictionary2table(dict1) == open(table1).read()

def test_funIn2out():		
	assert generalUtils.funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa') == 'input.method1.getNucAbu.csv'
	assert generalUtils.funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 1) == 'input.method1.gNA.csv'
	assert generalUtils.funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 2) == 'input.method1.geNuAb.csv'
	assert generalUtils.funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 4) =='input.method1.getNuclAbun.csv'
	assert generalUtils.funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', 'Extra') == 'input.method1.getNucAbuExtra.csv'

