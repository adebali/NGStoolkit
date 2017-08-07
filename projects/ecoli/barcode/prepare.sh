#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/ecoli/barcode/raw/

# DIR=/proj/sancarlb/HTSF/170630_UNC17-D00216_0455_ACB689ANXX/CLOCK

# gunzip -c ${DIR}/CPSMfdGMIP_CGTACG_S6_L002_R1_001.fastq.gz > ${TARGET}/CPSMfdGMIP.fastq
# gunzip -c ${DIR}/CPSMfdGPIM_GTTTCG_S7_L002_R1_001.fastq.gz > ${TARGET}/CPSMfdGPIM.fastq
# gunzip -c ${DIR}/CPSWTGMIP_GTGGCC_S8_L002_R1_001.fastq.gz > ${TARGET}/CPSWTGMIP.fastq
# gunzip -c ${DIR}/CPSWTGPIM_GTCCGC_S5_L002_R1_001.fastq.gz > ${TARGET}/CPSWTGPIM.fastq
# gunzip -c ${DIR}/CPSWTGPIP_GTGAAA_S9_L002_R1_001.fastq.gz > ${TARGET}/CPSWTGPIP.fastq
# gunzip -c ${DIR}/ppGppGMIP_GGTAGC_S4_L001_R1_001.fastq.gz > ${TARGET}/ppGppGMIP.fastq
# gunzip -c ${DIR}/ppGppGPIM_GAGTGG_S1_L001_R1_001.fastq.gz > ${TARGET}/ppGppGPIM.fastq
# gunzip -c ${DIR}/ppGppMfdGMIP_ATGAGC_S2_L001_R1_001.fastq.gz > ${TARGET}/ppGppMfdGMIP.fastq
# gunzip -c ${DIR}/ppGppMfdGPIM_ACTGAT_S3_L001_R1_001.fastq.gz > ${TARGET}/ppGppMfdGPIM.fastq

DIR=/proj/sancarlb/HTSF/170725_UNC17-D00216_0460_BCBF6MANXX/

# gunzip -c ${DIR}/CPSMG_TGACCA_S13_L003_R1_001.fastq.gz > ${TARGET}/CPSMG.fastq
# gunzip -c ${DIR}/CPSMI_ACAGTG_S16_L003_R1_001.fastq.gz > ${TARGET}/CPSMI.fastq
# gunzip -c ${DIR}/CPSppG_GCCAAT_S27_L005_R1_001.fastq.gz> ${TARGET}/CPSppG.fastq
# gunzip -c ${DIR}/CPSppI_CAGATC_S28_L005_R1_001.fastq.gz> ${TARGET}/CPSppI.fastq
gunzip -c ${DIR}/CPSppMG_ACTTGA_S29_L005_R1_001.fastq.gz> ${TARGET}/CPSppMG.fastq
# gunzip -c ${DIR}/CPSppMI_GATCAG_S30_L005_R1_001.fastq.gz> ${TARGET}/CPSppMI.fastq
# gunzip -c ${DIR}/CPSRminus_TAGCTT_S25_L005_R1_001.fastq.gz> ${TARGET}/CPSRminus.fastq
# gunzip -c ${DIR}/CPSRplus_GGCTAC_S26_L005_R1_001.fastq.gz> ${TARGET}/CPSRplus.fastq
# gunzip -c ${DIR}/CPSWTG_ATCACG_S15_L003_R1_001.fastq.gz> ${TARGET}/CPSWTG.fastq
# gunzip -c ${DIR}/CPSWTGI_TTAGGC_S14_L003_R1_001.fastq.gz> ${TARGET}/CPSWTGI.fastq
# gunzip -c ${DIR}/CPSWTI_CGATGT_S17_L003_R1_001.fastq.gz> ${TARGET}/CPSWTI.fastq

