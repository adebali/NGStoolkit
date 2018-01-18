#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/replication/XR/raw/

DIR=/proj/sancarlb/HTSF/170606_UNC17-D00216_0448_ACABYGANXX/NERSEQ

gunzip -c ${DIR}/HXA64A1_ATCACG_S49_L008_R1_001.fastq.gz > ${TARGET}/HXA64A1_ATCACG.fastq
gunzip -c ${DIR}/HXACA4_TGACCA_S44_L008_R1_001.fastq.gz > ${TARGET}/HXACA4_TGACCA.fastq
gunzip -c ${DIR}/HXE64A2_CGATGT_S45_L008_R1_001.fastq.gz > ${TARGET}/HXE64A2_CGATGT.fastq
gunzip -c ${DIR}/HXECA5_ACAGTG_S46_L008_R1_001.fastq.gz > ${TARGET}/HXECA5_ACAGTG.fastq
gunzip -c ${DIR}/HXL64A3_TTAGGC_S47_L008_R1_001.fastq.gz > ${TARGET}/HXL64A3_TTAGGC.fastq
gunzip -c ${DIR}/HXLCA6_GCCAAT_S50_L008_R1_001.fastq.gz > ${TARGET}/HXLCA6_GCCAAT.fastq

DIR=/proj/sancarlb/HTSF/170526_UNC18-D00493_0419_AC97LFANXX

gunzip -c ${DIR}/HXA64B7_CAGATC_S32_L008_R1_001.fastq.gz > ${TARGET}/HXA64B7_CAGATC.fastq
gunzip -c ${DIR}/HXACB10_TAGCTT_S33_L008_R1_001.fastq.gz > ${TARGET}/HXACB10_TAGCTT.fastq
gunzip -c ${DIR}/HXE64B8_ACTTGA_S34_L008_R1_001.fastq.gz > ${TARGET}/HXE64B8_ACTTGA.fastq
gunzip -c ${DIR}/HXECB11_GGCTAC_S35_L008_R1_001.fastq.gz > ${TARGET}/HXECB11_GGCTAC.fastq
gunzip -c ${DIR}/HXL64B9_GATCAG_S30_L008_R1_001.fastq.gz > ${TARGET}/HXL64B9_GATCAG.fastq
gunzip -c ${DIR}/HXLCB12_CTTGTA_S36_L008_R1_001.fastq.gz > ${TARGET}/HXLCB12_CTTGTA.fastq