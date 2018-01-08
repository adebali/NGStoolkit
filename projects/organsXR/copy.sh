#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --error=err.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

TARGET=/proj/yuchaojlab/sancarlab
DIR=/proj/sancarlb/users/ogun/sancarlabutils/projects

cp $DIR/organsXR/dataDir/raw/*fastq $TARGET/XRseq/
cp $DIR/organsXR/samples.csv $TARGET/XRseq/

cp $DIR/DamageSeq/organs/dataDir/raw/*fastq $TARGET/damageseq/
cp $DIR/DamageSeq/organs/samples.csv $TARGET/damageseq/

cp /nas/longleaf/home/adebali/ogun/sancarlabutils/projects/RNA-seq/organs/dataDir/raw/*fastq $TARGET/RNAseq/
cp $DIR/RNA-seq/organs/samples.csv $TARGET/RNAseq/

