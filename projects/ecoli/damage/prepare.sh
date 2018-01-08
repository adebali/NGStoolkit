#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/ecoli/damage/raw/

DIR=/proj/sancarlb/HTSF/171205_UNC32-K00270_0068_AHN7T3BBXX

gunzip -c ${DIR}/DA0NAK_GGCTAC_S1_L006_R1_001.fastq.gz > ${TARGET}/DA0NAK_GGCTAC_1.fastq
gunzip -c ${DIR}/DA0NAK_GGCTAC_S1_L006_R2_001.fastq.gz > ${TARGET}/DA0NAK_GGCTAC_2.fastq

gunzip -c ${DIR}/DA20CEL_ACTTGA_S2_L006_R1_001.fastq.gz > ${TARGET}/DA20CEL_ACTTGA_1.fastq
gunzip -c ${DIR}/DA20CEL_ACTTGA_S2_L006_R2_001.fastq.gz > ${TARGET}/DA20CEL_ACTTGA_2.fastq

gunzip -c ${DIR}/DA20NAK_CTTGTA_S3_L006_R1_001.fastq.gz > ${TARGET}/DA20NAK_CTTGTA_1.fastq
gunzip -c ${DIR}/DA20NAK_CTTGTA_S3_L006_R2_001.fastq.gz > ${TARGET}/DA20NAK_CTTGTA_2.fastq

gunzip -c ${DIR}/DA600CEL_TAGCTT_S4_L006_R1_001.fastq.gz > ${TARGET}/DA600CEL_TAGCTT_1.fastq
gunzip -c ${DIR}/DA600CEL_TAGCTT_S4_L006_R2_001.fastq.gz > ${TARGET}/DA600CEL_TAGCTT_2.fastq

gunzip -c ${DIR}/DA80CEL_GATCAG_S5_L006_R1_001.fastq.gz > ${TARGET}/DA80CEL_GATCAG_1.fastq
gunzip -c ${DIR}/DA80CEL_GATCAG_S5_L006_R2_001.fastq.gz > ${TARGET}/DA80CEL_GATCAG_2.fastq