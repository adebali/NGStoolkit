import os
import sys
from glob import glob
import generalUtils

dataDir = '/proj/sancarlb/users/ogun/scripts/projects/DamageSeq/hiSeq/dataDir'
chromHMMdir = os.path.join(dataDir, 'chromHMM')
filesToBeMerged = sorted(glob(os.path.join(dataDir, "*_mergedChromHMM.txt")))
output = os.path.join(dataDir,'chromHMM_allMerged.csv')
headerLine = os.path.join(dataDir,'chromHMM_allMerged_headers.csv')
headerOut = open(headerLine, 'w')
tmpFile= os.path.join(dataDir,'temp.txt')
kHomeDir = '/proj/sancarlb/users/ogun'

chromHMMfiles = generalUtils.sorted_nicely(glob(os.path.join(kHomeDir, 'chromHMM', '*_chromHMM')))

def run(code):
	print(code)
	if '--run' in sys.argv:
		os.system(code)

for i in range(len(filesToBeMerged)):
	file = filesToBeMerged[i]
	fileBaseName = os.path.basename(file)
	for chromFile in chromHMMfiles:
		chromBaseName = os.path.basename(chromFile).replace('_NHFL_chromHMM', '')
		headerOut.write(fileBaseName.replace('mergedChromHMM.txt', '') + chromBaseName +  '\t')
	if i == 0:
		code = 'cp ' + file + ' ' + tmpFile
		run(code)
	elif i != 0:
		code = 'paste ' + tmpFile + ' ' + file + ' >' + output
		run(code) 
		code = 'cp ' + output + ' ' + tmpFile
		run(code)
headerOut.write('\n')
headerOut.close()

code = 'cat ' + headerLine + ' ' + tmpFile + ' >' + output
os.system(code)
