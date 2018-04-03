#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

DATA_HOME=/nas/longleaf/home/adebali/ogun
PROJECT_NAME=lowXR
SUB_DIR=1802

TARGET=${DATA_HOME}/${PROJECT_NAME}
RAW=${TARGET}/raw
REPOSITORY_ROOT=/nas/longleaf/home/adebali/ogun/sancarlabutils
PROJECT_SCRIPT=${REPOSITORY_ROOT}/projects/${PROJECT_NAME}

echo "Creating directory at ${TARGET}"
mkdir -p ${TARGET}
echo "Creating directory at ${RAW}"
mkdir -p ${RAW}
echo "Creating directory at ${TARGET}/${SUB_DIR}"
mkdir -p ${TARGET}/${SUB_DIR}
echo "Creating directory at ${PROJECT_SCRIPT}/data"
mkdir -p ${PROJECT_SCRIPT}/log
echo "Creating shortcut to ${TARGET} as ${PROJECT_SCRIPT}/data"
ln -sf ${TARGET} ${PROJECT_SCRIPT}/data

DIR=/proj/sancarlb/HTSF/180207_UNC32-K00270_0079_BHN5CVBBXX
gunzip -c ${DIR}/JXR_10KA_TTAGGC_S5_L001_R1_001.fastq.gz > ${RAW}/JXR_10KA.fastq
gunzip -c ${DIR}/JXR_10KB_TGACCA_S4_L001_R1_001.fastq.gz > ${RAW}/JXR_10KB.fastq