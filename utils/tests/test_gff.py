from gff import *

def test_metaLine2metaDict():
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
    assert metaLine2metaDict(metaLine) == metaDict
    
def test_getGeneInformationFromGFFline():
    fullLine = "NC_000913.2	RefSeq	gene	190	255	.	+	.	ID=gene0;Dbxref=EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001,JW4367;locus_tag=b0001"
    assert getGeneInformationFromGFFline(fullLine, 'Name') == 'thrL'

