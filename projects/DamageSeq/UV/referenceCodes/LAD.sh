#!/usr/bin bash

DIR=/nas/longleaf/home/adebali/ogun/seq/hg19

WINDOW=100
FLANKING_SIZE=5000
OUTPUT_WORD="w100bp"

cd $DIR
echo getting center
bed2midPointFlankingBed.py -i LAD.bed -o LADcenter.bed -w $FLANKING_SIZE
echo making windows
bed2makeWindows.py -i LADcenter.bed -w $WINDOW | grep -v "^chrY" > LADcenter_w100bp_noY.bed
echo making bed6
cat LADcenter_w1kb_noY.bed | awk '{print $1"\t"$2"\t"$3"\t.\t"$4"\t."}' >LADcenter_w100bp_noY_bed6.bed