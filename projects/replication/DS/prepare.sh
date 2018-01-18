#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/replication/DS/raw/
DIR=/proj/sancarlb/HTSF/170526_UNC18-D00493_0419_AC97LFANXX

gunzip -c ${DIR}/HDA64A1_ATCACG_S18_L007_R1_001.fastq.gz > ${TARGET}/HDA64A1_ATCACG.fastq
gunzip -c ${DIR}/HDA64B19_GTGAAA_S19_L007_R1_001.fastq.gz > ${TARGET}/HDA64B19_GTGAAA.fastq
gunzip -c ${DIR}/HDACA6_GCCAAT_S20_L007_R1_001.fastq.gz > ${TARGET}/HDACA6_GCCAAT.fastq
gunzip -c ${DIR}/HDACB23_GAGTGG_S21_L007_R1_001.fastq.gz > ${TARGET}/HDACB23_GAGTGG.fastq
gunzip -c ${DIR}/HDE64A4_TGACCA_S22_L007_R1_001.fastq.gz > ${TARGET}/HDE64A4_TGACCA.fastq
gunzip -c ${DIR}/HDE64B20_GTGGCC_S29_L007_R1_001.fastq.gz > ${TARGET}/HDE64B20_GTGGCC.fastq
gunzip -c ${DIR}/HDECA10_TAGCTT_S24_L007_R1_001.fastq.gz > ${TARGET}/HDECA10_TAGCTT.fastq
gunzip -c ${DIR}/HDECB25_ACTGAT_S25_L007_R1_001.fastq.gz > ${TARGET}/HDECB25_ACTGAT.fastq
gunzip -c ${DIR}/HDL64A5_ACAGTG_S28_L007_R1_001.fastq.gz > ${TARGET}/HDL64A5_ACAGTG.fastq
gunzip -c ${DIR}/HDL64B22_CGTACG_S27_L007_R1_001.fastq.gz > ${TARGET}/HDL64B22_CGTACG.fastq
gunzip -c ${DIR}/HDLCA12_CTTGTA_S26_L007_R1_001.fastq.gz > ${TARGET}/HDLCA12_CTTGTA.fastq
gunzip -c ${DIR}/HDLCB27_ATTCCT_S23_L007_R1_001.fastq.gz > ${TARGET}/HDLCB27_ATTCCT.fastq