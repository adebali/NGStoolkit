import os
import sys
import generalUtils

fileList = sys.argv[1:len(sys.argv)]

latestOutput = fileList[0]

for i in range(1,len(fileList)):
	bed1 = latestOutput
	bed2 = fileList[i]
	if i == len(fileList):
		output = 'iterativeCov.' + 'final' + '.bed'
	else:
		output = 'iterativeCov.' + str(i) + '.bed'

	tempOutput = 'tempOutput.bed'

	codeList = [
		'bedtools',
		'coverage',
		'-counts',
		'-a', bed1, # Primary genomic locations
		'-b', bed2, # Read file, to get the read counts overlapping with the primary genomic locations
		'>', output
	]

	generalUtils.run(codeList)

	# codeList2 = [
	# 	'bedCount2rmNoHits.py',
	# 	tempOutput,
	# 	'>', output
	# ]
	#
	# generalUtils.run(codeList2)

	latestOutput = output
