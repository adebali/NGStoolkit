#!/bin/bash

#SBATCH --mem=8000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

HELA_REPLICATION_SEGMENTS_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53984/suppl/GSE53984_GSM923449_Helas3_Rep1_segments.bed.gz"
HELA_REPLICATION_SEGMENTS_GZ="repSegments.bed.gz"
HELA_REPLICATION_SEGMENTS="repSegments.bed"
HELA_REPLICATION_SEGMENTS_PRE="repSegments_prepared.bed"
HELA_CHMM_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz"
HELA_CHMM_GZ="chmm.bed.gz"
HELA_CHMM="chmm.bed"
HELA_CHMM_PRE="chmm_prepared.bed"
OUTPUT="repChmm.txt"

wget -O $HELA_REPLICATION_SEGMENTS_GZ $HELA_REPLICATION_SEGMENTS_URL
wget -O $HELA_CHMM_GZ $HELA_CHMM_URL

gunzip $HELA_REPLICATION_SEGMENTS_GZ
gunzip $HELA_CHMM_GZ

cut -f 1-6 $HELA_CHMM >$HELA_CHMM_PRE
addColumns.py -i $HELA_REPLICATION_SEGMENTS -c . .>$HELA_REPLICATION_SEGMENTS_PRE


bedtools intersect -a $HELA_CHMM_PRE -b $HELA_REPLICATION_SEGMENTS_PRE -wa -wb -f 1 | awk -v OFS='\t' '{print $1,$2,$3,$10"_"$4,$5,$6}' >$OUTPUT