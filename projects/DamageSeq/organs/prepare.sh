#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/scripts/projects/DamageSeq/organs/dataDir/raw

# DIR=/proj/sancarlb/HTSF/170623_UNC13-SN749_0634_BHJYFWBCXY
# gunzip -c ${DIR}/DLung_B_CGTACG_S9_L002_R1_001.fastq.gz > ${TARGET}/DLung_B_CGTACG.fastq
# gunzip -c ${DIR}/DLung_C_GAGTGG_S11_L002_R1_001.fastq.gz > ${TARGET}/DLung_C_GAGTGG.fastq
# gunzip -c ${DIR}/DSpleen_B_ACTGAT_S10_L002_R1_001.fastq.gz > ${TARGET}/DSpleen_B_ACTGAT.fastq
# gunzip -c ${DIR}/DSpleen_C_ATTCCT_S12_L002_R1_001.fastq.gz > ${TARGET}/DSpleen_C_ATTCCT.fastq

# DIR=/proj/sancarlb/HTSF/170623_UNC13-SN749_0634_BHJYFWBCXY
# gunzip -c ${DIR}/DLung_B_CGTACG_S9_L002_R1_001.fastq.gz > ${TARGET}/DLung_B_CGTACG.fastq
# gunzip -c ${DIR}/DLung_C_GAGTGG_S11_L002_R1_001.fastq.gz > ${TARGET}/DLung_C_GAGTGG.fastq
# gunzip -c ${DIR}/DSpleen_B_ACTGAT_S10_L002_R1_001.fastq.gz > ${TARGET}/DSpleen_B_ACTGAT.fastq
# gunzip -c ${DIR}/DSpleen_C_ATTCCT_S12_L002_R1_001.fastq.gz > ${TARGET}/DSpleen_C_ATTCCT.fastq

DIR=/proj/sancarlb/HTSF/170925_UNC31-K00269_0083_AHKMYFBBXX/NERSEQ
gunzip -c ${DIR}/HaCat_BG4_cis_B-_GAGTGG_S37_L006_R1_001.fastq.gz > ${TARGET}/HaCat_BG4_cis_B-_GAGTGG.fastq
gunzip -c ${DIR}/HaCat_BG4_cont_B-_GTGAAA_S42_L006_R1_001.fastq.gz > ${TARGET}/HaCat_BG4_cont_B-_GTGAAA.fastq
gunzip -c ${DIR}/HaCat_cis_B-_AGTTCC_S41_L006_R1_001.fastq.gz > ${TARGET}/HaCat_cis_B-_AGTTCC.fastq
gunzip -c ${DIR}/HaCat_cont_B-_AGTCAA_S40_L006_R1_001.fastq.gz > ${TARGET}/HaCat_cont_B-_AGTCAA.fastq
