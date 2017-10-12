#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com
#SBATCH --error=err.txt


HOME=/nas/longleaf/home/adebali/ogun/DamageSeq
PROJECT_NAME=ratio
TARGET=${HOME}/${PROJECT_NAME}
RAW=${TARGET}/raw
SCRIPTS=/nas/longleaf/home/adebali/ogun/scripts
PROJECT_SCRIPT=${SCRIPTS}/projects/${PROJECT_NAME}
DIR=/proj/sancarlb/HTSF/170925_UNC31-K00269_0083_AHKMYFBBXX/NERSEQ

gunzip -c ${DIR}/X20J7-_CAGATC_S39_L006_R1_001.fastq.gz > ${RAW}/X20J7-_CAGATC.fastq
ln -s ${RAW}/X20J7-_CAGATC.fastq ${RAW}/X20J7-_CAGATC_ec.fastq
gunzip -c ${DIR}/XNOJ5-_ACAGTG_S38_L006_R1_001.fastq.gz> ${RAW}/XNOJ5-_ACAGTG.fastq
ln -s ${RAW}/XNOJ5-_ACAGTG.fastq ${RAW}/XNOJ5-_ACAGTG_ec.fastq




