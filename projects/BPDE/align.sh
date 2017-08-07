#!/bin/bash

#SBATCH --mem=64000
#SBATCH -n 4
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com


# bowtie -x -q --nomaqround --phred33-quals -m 4 -n 2 -e 70 -l 20 --best -p 4 --seed=123 /proj/sancarlb/users/wentao/HG19_UCSC/Bt1Male/male_bt1.hg19 -f /proj/sancarlb/users/wentao/20150428jinchuancpd4h10j/NHFCPD4H11_out.fasta -S dataDir/0308/NHFCPD4H11.sam
# bowtie -x -q --nomaqround --phred33-quals -m 4 -n 2 -e 70 -l 20 --best -p 4 --seed=123 /proj/sancarlb/users/wentao/HG19_UCSC/Bt1Male/male_bt1.hg19 -f /proj/sancarlb/users/wentao/20150428jinchuancpd4h10j/NHFCPD4H6_out.fasta -S dataDir/0308/NHFCPD4H6.sam
# bowtie -x -q --nomaqround --phred33-quals -m 4 -n 2 -e 70 -l 20 --best -p 4 --seed=123 /proj/sancarlb/users/wentao/HG19_UCSC/Bt1Male/male_bt1.hg19 -f /proj/sancarlb/users/wentao/20150428jinchuancpd4h10j/1H6-4-10_out.fasta -S dataDir/0308/1H6-4-10.sam
# bowtie -x -q --nomaqround --phred33-quals -m 4 -n 2 -e 70 -l 20 --best -p 4 --seed=123 /proj/sancarlb/users/wentao/HG19_UCSC/Bt1Male/male_bt1.hg19 -f /proj/sancarlb/users/wentao/20150428jinchuancpd4h10j/1H6-4-6_out.fasta -S dataDir/0308/1H6-4-6.sam


# samtools view -b -S -o dataDir/0308/NHFCPD4H11.sam.bam dataDir/0308/NHFCPD4H11.sam
# samtools view -b -S -o dataDir/0308/NHFCPD4H6.sam.bam dataDir/0308/NHFCPD4H6.sam
# samtools view -b -S -o dataDir/0308/1H6-4-10.sam.bam dataDir/0308/1H6-4-10.sam
# samtools view -b -S -o dataDir/0308/1H6-4-6.sam.bam dataDir/0308/1H6-4-6.sam

# bedtools bamtobed -i dataDir/0308/NHFCPD4H11.sam.bam >NHFCPD4H11.bed
# bedtools bamtobed -i dataDir/0308/NHFCPD4H6.sam.bam >NHFCPD4H6.bed
# bedtools bamtobed -i dataDir/0308/1H6-4-10.sam.bam >1H6-4-10.bed
# bedtools bamtobed -i dataDir/0308/1H6-4-6.sam.bam >1H6-4-6.bed