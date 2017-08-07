#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com


HOME=/nas/longleaf/home/adebali/ogun
PROJECT_NAME=plant

TARGET=${HOME}/${PROJECT_NAME}
RAW=${TARGET}/raw

SCRIPTS=/nas/longleaf/home/adebali/ogun/scripts
PROJECT_SCRIPT=${SCRIPTS}/projects/${PROJECT_NAME}

# mkdir -p ${TARGET}
# mkdir -p ${RAW}
# mkdir -p ${TARGET}/0707
# mkdir -p ${PROJECT_SCRIPT}/dataDir
#ln -s ${TARGET} ${PROJECT_SCRIPT}/log

# DIR=/proj/sancarlb/HTSF/170630_UNC17-D00216_0455_ACB689ANXX/NERSEQ
# gunzip -c ${DIR}/ONUR1_ATCACG_S13_L003_R1_001.fastq.gz > ${RAW}/ONUR1.fastq

DIR=/proj/sancarlb/HTSF/170725_UNC17-D00216_0460_BCBF6MANXX
gunzip -c ${DIR}/ATCPDrep2_ATCACG_S31_L006_R1_001.fastq.gz > ${RAW}/ATCPDrep2.fastq




