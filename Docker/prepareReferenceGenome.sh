#!/usr/bin/bash

cd /data

wget -r --no-parent -A 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz' ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/

cd ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna
mkdir -p /data/genomes/GRCh38/Bowtie2

bowtie2-build --threads 4 -f \
Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz,\
Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz \
/data/genomes/GRCh38/Bowtie2/genome

zcat *gz >/data/genomes/GRCh38/genome.fa && cd /data/genomes/GRCh38 && samtools faidx genome.fa
wget -O /data/genomes/GRCh38/genes.gff3.gz ftp://ftp.ensembl.org/pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz && gunzip /data/genomes/GRCh38/genes.gff3.gz && gff2bed.py -i /data/genomes/GRCh38/genes.gff3 -r gene -o /data/genomes/GRCh38/genes.bed && rm /data/genomes/GRCh38/genes.gff3