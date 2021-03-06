import os

# Requirements
# bedtools in PATH

scriptDir = os.path.dirname(os.path.realpath(__file__))
kTestFilesDir = os.path.join(scriptDir, 'testFiles')
kBedExample = os.path.join(kTestFilesDir, 'bedExample.bed')
kBedFasta = os.path.join(kTestFilesDir, 'bed.fa')
tempOutput = os.path.join(kTestFilesDir, 'temp.out')


def test_getfasta():
	os.system('rm ' + kBedFasta + '.fai')
	os.system('bedtools getfasta -fi ' + kBedFasta + ' -bed ' + kBedExample + ' -fo ' + tempOutput + ' -s')
	results = open(tempOutput, 'r').read()
	expectedResults = '>chr1:10-20(+)\nAtttAAGCCC\n>chr1:10-20(-)\nGGGCTTaaaT\n>chr2:0-20(+)\nGAAAACAAAAATGCATGXAA\n>chr2:0-20(-)\nTTXCATGCATTTTTGTTTTC\n>chr3:1-5(+)\nAAAA\n>chr3:1-5(-)\nTTTT\n'
	assert results == expectedResults
	os.system('rm ' + tempOutput)

