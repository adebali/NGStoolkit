#!/usr/bin/env bash

#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com
#SBATCH --output=chromatin.out
#SBATCH --error=chromatin.err

OUT_DIR=/nas/longleaf/home/adebali/ogun/seq/TAIR10/Liu2017

# rm -f ${OUT_DIR}/chrStates.bed
# for i in `seq 1 36`;
# do
#    wget -qO- http://systemsbiology.cau.edu.cn/chromstates/download/At_segments_S${i} | addColumns.py -c . . | sed 's/chr//g' >>${OUT_DIR}/chrStates_unsorted.bed
# done  

#  sort -k1,1 -k2,2n -k3,3n ${OUT_DIR}/chrStates_unsorted.bed | grep -Pv "^C" | grep -Pv "^M" >${OUT_DIR}/chrStates.bed
#  rm ${OUT_DIR}/chrStates_unsorted.bed
#  bed2removeChromosomeEdges.py -i ${OUT_DIR}/chrStates.bed -o ${OUT_DIR}/chrStates_noEdge.bed -g ${OUT_DIR}/../genome.fa.fai --fixed -l 1000

 bedtools shuffle -i ${OUT_DIR}/chrStates_noEdge.bed -g ${OUT_DIR}/../genome.fa.fai | sort -k1,1 -k2,2n -k3,3n | bed2removeChromosomeEdges.py -g ${OUT_DIR}/../genome.fa.fai --fixed -l 1000 >${OUT_DIR}/chrStates_noEdge_shuffled.bed