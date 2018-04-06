#!/usr/bin/env bash

echo "Define the sample basename (eg. you have fastq file named as runSample.fastq)"
SAMPLE=runSample

SAMPLE=SRR1976056
fastq-dump $SAMPLE

echo "Cut adapter"
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq

echo "Align with the reference genome"
BOWTIE2_IND=/data/genomes/GCRh38/Bowtie2/genome
bowtie2 -p 4 -x $BOWTIE2_IND -U ${SAMPLE}_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam

echo "Convert to bam"
samtools view -q 20 -b -o ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.sam

echo "Convert to bed"
bedtools bamtobed -i ${SAMPLE}_cutadapt.bam >${SAMPLE}_cutadapt.bed

echo "Sort bed file. Use -u to remove duplicates"
sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE}_cutadapt.bed >${SAMPLE}_cutadapt_sorted.bed

echo "Count number of mapped reads"
grep -c "^" ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_readCount.txt 

echo "Get the read length distribution of the aligned and deduplicated reads"
awk '{print $3-$2}' ${SAMPLE}_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}'

echo "Get the certain-sized reads (eg 26)"
awk '{ if ($3-$2 == 26) { print } }' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_26.bed

echo "Retrieve sequences in fasta format"
bedtools getfasta -fi /data/genomes/GRCh38/genome.fa -bed ${SAMPLE}_cutadapt_sorted_26.bed -fo ${SAMPLE}_cutadapt_sorted_26.fa

echo "Get the dinucleotide content of the reads"
fa2kmerAbundanceTable.py -i ${SAMPLE}_cutadapt_sorted_26.fa -k 2 -o ${SAMPLE}_cutadapt_sorted_26_dinucleotideTable.txt

echo "Convert bed to bedgraph"
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted.bed -g /data/genomes/GRCh38/genome.fa.fai -bg -scale $(cat ${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print 1000000/$1}') >${SAMPLE}_cutadapt_sorted.bdg

echo "Convert bedgraph to bigwig"
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted.bdg /data/genomes/GRCh38/genome.fa.fai ${SAMPLE}_cutadapt_sorted.bw

echo "Count read values for transcribed (TS) and nontranscribed (NTS) strands"
bedtools intersect -a /data/genomes/GRCh38/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -S -F 0.5 > ${SAMPLE}_cutadapt_sorted_TScount.txt
bedtools intersect -a /data/genomes/GRCh38/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -s -F 0.5 > ${SAMPLE}_cutadapt_sorted_NTScount.txt