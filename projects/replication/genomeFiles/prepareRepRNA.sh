#!/bin/bash

#SBATCH --mem=8000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

HELA_REPLICATION_SEGMENTS_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53984/suppl/GSE53984_GSM923449_Helas3_Rep1_segments.bed.gz"
HELA_REPLICATION_SEGMENTS_GZ="repSegments.bed.gz"
HELA_REPLICATION_SEGMENTS="repSegments.bed"
HELA_REPLICATION_SEGMENTS_PRE="repSegments_prepared.bed"

HELA_RNAcount_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM958nnn/GSM958734/suppl/GSM958734_hg19_wgEncodeCaltechRnaSeqHuvecR2x75Il200TranscriptGencV3cRep1V3.gtf.gz"

HELA_RNAcount_GZ="RNAcount.gtf.gz"
HELA_RNAcount="RNAcount.gtf"
HELA_RNAcount_BED="RNAcount_prepared.bed"
OUTPUT="/proj/sancarlb/users/ogun/seq/hg19/repRNAcount.bed"

# wget -O $HELA_REPLICATION_SEGMENTS_GZ $HELA_REPLICATION_SEGMENTS_URL
# wget -O $HELA_RNAcount_GZ $HELA_RNAcount_URL

# gunzip $HELA_REPLICATION_SEGMENTS_GZ
# gunzip $HELA_RNAcount_GZ

# cat $HELA_RNAcount | grep -P "\\ttranscript\\t" | gtf2bed | cut -f 1-6 | grep -vP "\\t\.$" >$HELA_RNAcount_BED

bed2getLongestById.py -i $HELA_RNAcount_BED | bedtools intersect -a 'stdin' -b $HELA_REPLICATION_SEGMENTS_PRE -wa -wb -f 1 | awk -v OFS='\t' '{print $1,$2,$3,$10,$5,$6}' | sort -u -k1,1 -k2,2n -k3,3n >$OUTPUT