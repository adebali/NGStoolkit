#!/usr/bin/env bash

# Define the sample basename (eg. you have fastq file named as runSample.fastq)
SAMPLE=runSample

# Cut adapter
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG 
-o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq

# Align with the reference genome
bowtie2 -x /data/genomes/GCRh38/Bowtie2/genome -U ${SAMPLE}_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam

# Convert to bam
samtools view -b -q 20 -o ${SAMPLE}_cutadapt.sam ${SAMPLE}_cutadapt.bam

# Convert to bed
bedtools bamtobed -i ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.bed

# Sort bed file. Use -u to remove duplicates
sort -u -k1,1 -k2,2n –k3,3n ${SAMPLE}_cutadapt.bed > ${SAMPLE}_cutadapt_sorted.bed

# Count number of mapped reads
grep -c “^” ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_readCount.txt 

# Get the read length distribution of the aligned and deduplicated reads
awk '{print $3-$2}' ${SAMPLE}_cutadapt_sorted.bed | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' | sort -k1,1n

# Get the certain-sized reads (eg 27)
awk '{ if ($3-$2 == 27) { print } }' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_27.bed

# Retrieve sequences in fasta format
bedtools getfasta -fi GCRh38.fa -bed ${SAMPLE}_cutadapt_sorted_27.bed -fo ${SAMPLE}_cutadapt_sorted_27.fa

# Get the dinucleotide content of the reads
fa2kmerAbundanceTable.py -i ${SAMPLE}_cutadapt_sorted_27.fa -k 2 -o ${SAMPLE}_cutadapt_sorted_27_dinucleotideTable.txt

# Convert bed to bedgraph
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted.bed -g genome.fa.fai -bg -scale $(${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print 1000000/$1}') >${SAMPLE}_cutadapt_sorted.bdg

# Convert bedgraph to bigwig
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted.bdg genome.fa.fai ${SAMPLE}_cutadapt_sorted.bw

# Count read values for transcribed (TS) and nontranscribed (NTS) strands
bedtools intersect -a geneList.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -S -F 0.5 > ${SAMPLE}_cutadapt_sorted_ TScount.txt
bedtools intersect -a geneList.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -s -F 0.5 > ${SAMPLE}_cutadapt_sorted_NTScount.txt