import os


def downloadChromosome(accession):
    print("Downloading " + accession)
    os.system("efetch -db nuccore -format fasta -id " + accession + ' >' + '~/ogun/ENCODE/hg19/nucleosome/' + accession + '.fa')

def catFasta(accession):
    os.system("cat ~/ogun/ENCODE/hg19/nucleosome/" + accession + ".fa >> ~/ogun/ENCODE/hg19/nucleosome/hg19.fa" )

filein = open('nucleosome_data_referenceGenome.txt', 'r')
os.system("touch ~/ogun/ENCODE/hg19/nucleosome/hg19.fa" )



for line in filein:
    accession = line.strip()
    #downloadChromosome(accession)
    catFasta(accession)

