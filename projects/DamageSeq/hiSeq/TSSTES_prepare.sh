#!/bin/bash

#SBATCH --mem=64000
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com


SCORE=300
DISTANCE=6000
#MINLENGTH=10000
MINLENGTH=15000
GENEBODY=15000
FLANKING=5000
WINDOW=100

#cat hg19_expression_noHeader.bed | awk "\$5>${SCORE}" >hg19_expression_noHeader_sgt${SCORE}.bed
bed2removeNeighbors.py -i hg19_expression_noHeader_sgt${SCORE}.bed -d $DISTANCE  | awk '$8=="protein_coding"' | awk "\$3-\$2>${MINLENGTH}" > hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}.bed
cat hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}.bed | bedExpression2TxnSites.py -s start -g ~/ogun/seq/hg19/HG19_UCSC/genome.chrSizes.bed -u $FLANKING -d $GENEBODY >hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TSS.bed
cat hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}.bed | bedExpression2TxnSites.py -s end -g ~/ogun/seq/hg19/HG19_UCSC/genome.chrSizes.bed -u $GENEBODY -d $FLANKING >hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TES.bed

bed2makeWindows.py -i hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TSS.bed -w $WINDOW >hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TSS_win${WINDOW}.bed
bed2makeWindows.py -i hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TES.bed -w $WINDOW >hg19_expression_noHeader_sgt${SCORE}_noNeighIn${DISTANCE}_min${MINLENGTH}_TES_win${WINDOW}.bed