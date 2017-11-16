#!/usr/bin/env bash

#SBATCH --mem=128000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com
#SBATCH --output=TT.out
#SBATCH --error=TT.err

# bedtools slop -i ~/ogun/seq/TAIR10/Hetzel_sorted_noEdge10K.bed -b 10000 -g /proj/seq/data/TAIR10_NCBI/Sequence/WholeGenomeFasta/genome.fa.fai |  bed2motifBed.py -fasta ~/ogun/seq/TAIR10/genome.fa -motif TT | sort -k1,1 -k2,2n -k3,3n >~/ogun/seq/TAIR10/Hetzel_sorted_noEdge10K_s10K_TT.bed

bed2motifBed.py -i ~/ogun/seq/TAIR10/genome.bed -fasta ~/ogun/seq/TAIR10/genome.fa -motif TT | sort -k1,1 -k2,2n -k3,3n >~/ogun/seq/TAIR10/genome_TT.bed