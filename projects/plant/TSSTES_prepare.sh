#!/bin/bash

#SBATCH --mem=16000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com


# ROOT=/nas/longleaf/home/adebali/ogun/seq/TAIR9
# FILE_BASE=transcripts
ROOT=/nas/longleaf/home/adebali/ogun/seq/TAIR10
FILE_BASE=Hetzel2016_noHeader_sorted_NascentRNA
DISTANCE=2000
MINLENGTH=2000
DISTANCE_WORD=noOv2K
MINLENGTH_WORD=gt2K
WINDOW=100
WINDOW_WORD=w100
GENEBODY=$MINLENGTH
EXTENSION=8000

# bed2removeNeighbors.py -i ${ROOT}/${FILE_BASE}.bed -d $DISTANCE | awk "\$3-\$2>${MINLENGTH}" >${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}.bed

# bedExpression2TxnSites.py -i ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}.bed -u $EXTENSION -d $GENEBODY -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s start >${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_TSS.bed
# bedExpression2TxnSites.py -i ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}.bed -u $GENEBODY -d $EXTENSION -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s end >${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_TES.bed

# bed2makeWindows.py -i ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_TSS.bed -w $WINDOW >${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_${WINDOW_WORD}_TSS.bed
# bed2makeWindows.py -i ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_TES.bed -w $WINDOW >${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_${WINDOW_WORD}_TES.bed



HIT_NUMBER="$(grep -c '^' ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}.bed)"
declare -i QUART_I
declare -i QUART_II
declare -i QUART_III
((QUART_I=$HIT_NUMBER/4))
((QUART_II=($HIT_NUMBER/4)*2))
((QUART_III=($HIT_NUMBER/4)*3))
echo $QUART_I

FILE_SORT_SCORE=${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}_sortByScore.bed
sort -k5,5n ${ROOT}/${FILE_BASE}_${DISTANCE_WORD}_${MINLENGTH_WORD}.bed >${FILE_SORT_SCORE}

awk "NR>=1 && NR<$QUART_I" ${FILE_SORT_SCORE} | sort -k1,1 -k2,2n -k3,3n > ${FILE_SORT_SCORE}_Q4.bed
awk "NR>=${QUART_I} && NR<${QUART_II}" ${FILE_SORT_SCORE} | sort -k1,1 -k2,2n -k3,3n > ${FILE_SORT_SCORE}_Q3.bed
awk "NR>=${QUART_II} && NR<${QUART_III}" ${FILE_SORT_SCORE} | sort -k1,1 -k2,2n -k3,3n > ${FILE_SORT_SCORE}_Q2.bed
awk "NR>=${QUART_III} && NR<=${HIT_NUMBER}" ${FILE_SORT_SCORE} | sort -k1,1 -k2,2n -k3,3n > ${FILE_SORT_SCORE}_Q1.bed


bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q1.bed -u $EXTENSION -d $GENEBODY -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s start >${FILE_SORT_SCORE}_Q1.bed_TSS.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q1.bed -u $GENEBODY -d $EXTENSION -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s end >${FILE_SORT_SCORE}_Q1.bed_TES.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q2.bed -u $EXTENSION -d $GENEBODY -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s start >${FILE_SORT_SCORE}_Q2.bed_TSS.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q2.bed -u $GENEBODY -d $EXTENSION -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s end >${FILE_SORT_SCORE}_Q2.bed_TES.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q3.bed -u $EXTENSION -d $GENEBODY -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s start >${FILE_SORT_SCORE}_Q3.bed_TSS.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q3.bed -u $GENEBODY -d $EXTENSION -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s end >${FILE_SORT_SCORE}_Q3.bed_TES.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q4.bed -u $EXTENSION -d $GENEBODY -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s start >${FILE_SORT_SCORE}_Q4.bed_TSS.bed
bedExpression2TxnSites.py -i ${FILE_SORT_SCORE}_Q4.bed -u $GENEBODY -d $EXTENSION -g /proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9/Sequence/WholeGenomeFasta/genome.fa.fai -s end >${FILE_SORT_SCORE}_Q4.bed_TES.bed

bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q1.bed_TSS.bed -w $WINDOW >${FILE_SORT_SCORE}_Q1.bed_${WINDOW_WORD}_TSS.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q1.bed_TES.bed -w $WINDOW >${FILE_SORT_SCORE}_Q1.bed_${WINDOW_WORD}_TES.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q2.bed_TSS.bed -w $WINDOW >${FILE_SORT_SCORE}_Q2.bed_${WINDOW_WORD}_TSS.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q2.bed_TES.bed -w $WINDOW >${FILE_SORT_SCORE}_Q2.bed_${WINDOW_WORD}_TES.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q3.bed_TSS.bed -w $WINDOW >${FILE_SORT_SCORE}_Q3.bed_${WINDOW_WORD}_TSS.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q3.bed_TES.bed -w $WINDOW >${FILE_SORT_SCORE}_Q3.bed_${WINDOW_WORD}_TES.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q4.bed_TSS.bed -w $WINDOW >${FILE_SORT_SCORE}_Q4.bed_${WINDOW_WORD}_TSS.bed
bed2makeWindows.py -i ${FILE_SORT_SCORE}_Q4.bed_TES.bed -w $WINDOW >${FILE_SORT_SCORE}_Q4.bed_${WINDOW_WORD}_TES.bed