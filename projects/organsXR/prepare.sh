#!/usr/bin/env bash

DIR=/proj/sancarlb/HTSF/170215_UNC17-D00216_0422_AC8BE5ANXX
TARGET=/nas/longleaf/home/adebali/ogun/organsXR/raw


gunzip -c ${DIR}/KIDNEY-CISPL_ACTTGA_S19_L005_R1_001.fastq.gz > ${TARGET}/KIDNEY-CISPL_ACTTGA.fastq
gunzip -c ${DIR}/LIVERA_AGTTCC_S21_L005_R1_001.fastq.gz > ${TARGET}/LIVERA_AGTTCC.fastq
gunzip -c ${DIR}/LIVER-CISPL_GCCAAT_S20_L005_R1_001.fastq.gz > ${TARGET}/LIVER-CISPL_GCCAAT.fastq
gunzip -c ${DIR}/LUNG-CISPL_CAGATC_S22_L005_R1_001.fastq.gz > ${TARGET}/LUNG-CISPL_CAGATC.fastq
gunzip -c ${DIR}/SPLEEN-CISPL_GATCAG_S23_L005_R1_001.fastq.gz > ${TARGET}/SPLEEN-CISPL_GATCAG.fastq

