#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

DATA_HOME=/nas/longleaf/home/adebali/ogun
PROJECT_NAME=G4/DS
SUB_DIR=1712

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

DIR=/proj/sancarlb/HTSF/170925_UNC31-K00269_0083_AHKMYFBBXX/NERSEQ

gunzip -c ${DIR}/HaCat_A-_AGTTCC_S12_L002_R1_001.fastq.gz  > ${RAW}/HaCat_A-_AGTTCC.fastq           
gunzip -c ${DIR}/HaCat_B-_ATGTCA_S13_L002_R1_001.fastq.gz  > ${RAW}/HaCat_B-_ATGTCA.fastq           
gunzip -c ${DIR}/HaCat_BG4_cis_B-_GAGTGG_S37_L006_R1_001.fastq.gz  > ${RAW}/HaCat_BG4_cis_B-_GAGTGG.fastq   
gunzip -c ${DIR}/HaCat_BG4_cont_B-_GTGAAA_S42_L006_R1_001.fastq.gz  > ${RAW}/HHaCat_BG4_cont_B-_GTGAAA.fastq  
gunzip -c ${DIR}/HaCat_cis_B-_AGTTCC_S41_L006_R1_001.fastq.gz  > ${RAW}/HaCat_cis_B-_AGTTCC.fastq       
gunzip -c ${DIR}/HaCat_cont_B-_AGTCAA_S40_L006_R1_001.fastq.gz  > ${RAW}/HaCat_cont_B-_AGTCAA.fastq

DIR=/proj/sancarlb/HTSF/170630_UNC17-D00216_0455_ACB689ANXX/NERSEQ

gunzip -c ${DIR}/HaCat_BG4_cis_CCGTCC_S61_L008_R1_001.fastq.gz  > ${RAW}/HaCat_BG4_cis_CCGTCC.fastq
gunzip -c ${DIR}/HaCat_BG4_cont_ATGTCA_S59_L008_R1_001.fastq.gz  > ${RAW}/HaCat_BG4_cont_ATGTCA.fastq
gunzip -c ${DIR}/HaCat_cis_AGTTCC_S60_L008_R1_001.fastq.gz  > ${RAW}/HaCat_cis_AGTTCC.fastq
gunzip -c ${DIR}/HaCat_cont_AGTCAA_S62_L008_R1_001.fastq.gz  > ${RAW}/HaCat_cont_AGTCAA.fastq
