


import argparse

parser = argparse.ArgumentParser(description='returns abundance of nucleotide(s) near by a motif of interest')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
parser.add_argument('-m', required= True, help='motif')
parser.add_argument('-p', required= True, help='position of interest (eg -1, 2 etc)')
parser.add_argument('-l', required=  False, help='sequence length of interest')
parser.add_argument('--percentage', action='store_true', help = 'Write percentages instead if actual counts')

args = parser.parse_args()

sequence = 'ATATTGCCGGCGAATTGCGTGGCGTTAAGGGATTTGAGATGTGATTTGAGAGATTGC'

lengthOfMotif = len(motif)
positionOfInterest = -1
lengthOfInterest = 2
myDict = {}

for i in range(len(sequence)):
	if positionOfInterest > 0:
		startPosition = i + lengthOfMotif
	elif positionOfInterest < 0:
		startPosition = i + positionOfInterest

	subseq = sequence[i : i + lengthOfMotif]
	if subseq == motif:
		nearSubstring = sequence[startPosition : startPosition + lengthOfInterest]
		if nearSubstring in myDict.keys():
			myDict[nearSubstring] += 1
		else:
			myDict[nearSubstring] = 1

print(myDict)
