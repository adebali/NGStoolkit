import os
import re
import itertools
import unittest

def metaLine2metaDict(metaLine):
	metaDict = {}
	fields = metaLine.split(';')
	for field in fields:
		fl = field.split('=')
		subfieldName = fl[0]
		fieldInfo = fl[1]
		metaDict[subfieldName] = fieldInfo
	return metaDict

def getGeneInformationFromGFFline(line, field):
	result = False
	if not line.startswith('#'):
		ll = line.strip().split('\t')
		if len(ll) > 2 and ll[2] == 'gene':
			start = ll[3]
			end = ll[4]
			strand = ll[6]
			metaLine = ll[8]
			metaDict= metaLine2metaDict(metaLine)
			result = metaDict[field]
	return result

class tests(unittest.TestCase):
	def test_metaLine2metaDict(self):
		metaLine = 'ID=gene0;Dbxref=EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001,JW4367;locus_tag=b0001'
		metaDict = { 
			'ID':'gene0',
			'Dbxref':'EcoGene:EG11277,GeneID:944742',
			'Name':'thrL',
			'gbkey':'Gene',
			'gene':'thrL',
			'gene_biotype':'protein_coding',
			'gene_synonym':'ECK0001,JW4367',
			'locus_tag':'b0001'
		}
		self.assertEqual(metaLine2metaDict(metaLine), metaDict)
		
	def test_getGeneInformationFromGFFline(self):
		fullLine = "NC_000913.2	RefSeq	gene	190	255	.	+	.	ID=gene0;Dbxref=EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001,JW4367;locus_tag=b0001"
		self.assertEqual(getGeneInformationFromGFFline(fullLine, 'Name'), 'thrL')


if __name__ == "__main__":
	unittest.main()

