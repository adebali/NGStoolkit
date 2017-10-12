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

# DIR=/proj/sancarlb/HTSF/170725_UNC17-D00216_0460_BCBF6MANXX/

# gunzip -c ${DIR}/CPSMG_TGACCA_S13_L003_R1_001.fastq.gz > ${TARGET}/CPSMG.fastq
# gunzip -c ${DIR}/CPSMI_ACAGTG_S16_L003_R1_001.fastq.gz > ${TARGET}/CPSMI.fastq
# gunzip -c ${DIR}/CPSppG_GCCAAT_S27_L005_R1_001.fastq.gz> ${TARGET}/CPSppG.fastq
# gunzip -c ${DIR}/CPSppI_CAGATC_S28_L005_R1_001.fastq.gz> ${TARGET}/CPSppI.fastq
# gunzip -c ${DIR}/CPSppMG_ACTTGA_S29_L005_R1_001.fastq.gz> ${TARGET}/CPSppMG.fastq
# gunzip -c ${DIR}/CPSppMI_GATCAG_S30_L005_R1_001.fastq.gz> ${TARGET}/CPSppMI.fastq
# gunzip -c ${DIR}/CPSRminus_TAGCTT_S25_L005_R1_001.fastq.gz> ${TARGET}/CPSRminus.fastq
# gunzip -c ${DIR}/CPSRplus_GGCTAC_S26_L005_R1_001.fastq.gz> ${TARGET}/CPSRplus.fastq
# gunzip -c ${DIR}/CPSWTG_ATCACG_S15_L003_R1_001.fastq.gz> ${TARGET}/CPSWTG.fastq
# gunzip -c ${DIR}/CPSWTGI_TTAGGC_S14_L003_R1_001.fastq.gz> ${TARGET}/CPSWTGI.fastq
# gunzip -c ${DIR}/CPSWTI_CGATGT_S17_L003_R1_001.fastq.gz> ${TARGET}/CPSWTI.fastq

# DIR=/proj/sancarlb/HTSF/170828_UNC13-SN749_0643_AHT27MBCXY

# gunzip -c ${DIR}/CPSA0_CTTGTA_S8_L002_R1_001.fastq.gz >${TARGET}/CPSA0.fastq
# gunzip -c ${DIR}/CPSA10_AGTTCC_S9_L002_R1_001.fastq.gz >${TARGET}/CPSA10.fastq
# gunzip -c ${DIR}/CPSA30_ATGTCA_S11_L002_R1_001.fastq.gz >${TARGET}/CPSA30.fastq
# gunzip -c ${DIR}/CPSA3_AGTCAA_S10_L002_R1_001.fastq.gz >${TARGET}/CPSA3.fastq
# gunzip -c ${DIR}/CPSA90_CCGTCC_S12_L002_R1_001.fastq.gz >${TARGET}/CPSA90.fastq
# gunzip -c ${DIR}/CPSB0_GTAGAG_S1_L001_R1_001.fastq.gz >${TARGET}/CPSB0.fastq
# gunzip -c ${DIR}/CPSB10_GTGAAA_S2_L001_R1_001.fastq.gz >${TARGET}/CPSB10.fastq
# gunzip -c ${DIR}/CPSB30_GTGGCC_S4_L001_R1_001.fastq.gz >${TARGET}/CPSB30.fastq
# gunzip -c ${DIR}/CPSB3_GTCCGC_S3_L001_R1_001.fastq.gz >${TARGET}/CPSB3.fastq
# gunzip -c ${DIR}/CPSB90_GTTTCG_S5_L001_R1_001.fastq.gz >${TARGET}/CPSB90.fastq
# gunzip -c ${DIR}/CPSC0_CGTACG_S13_L002_R1_001.fastq.gz >${TARGET}/CPSC0.fastq
# gunzip -c ${DIR}/CPSC90_GAGTGG_S14_L002_R1_001.fastq.gz >${TARGET}/CPSC90.fastq
# gunzip -c ${DIR}/CPSD0_GGTAGC_S6_L001_R1_001.fastq.gz >${TARGET}/CPSD0.fastq
# gunzip -c ${DIR}/CPSD90_ACTGAT_S7_L001_R1_001.fastq.gz >${TARGET}/CPSD90.fastq

DIR=/proj/sancarlb/HTSF/170901_UNC18-D00493_0441_ACBFRBANXX/CLOCK

gunzip -c ${DIR}/CPSA0B_ATCACG_S34_L005_R1_001.fastq.gz >${TARGET}/CPSA0B.fastq
gunzip -c ${DIR}/CPSA10B_TTAGGC_S35_L005_R1_001.fastq.gz >${TARGET}/CPSA10B.fastq
gunzip -c ${DIR}/CPSA30B_TGACCA_S36_L005_R1_001.fastq.gz >${TARGET}/CPSA30B.fastq
gunzip -c ${DIR}/CPSA3B_GGCTAC_S37_L005_R1_001.fastq.gz >${TARGET}/CPSA3B.fastq
gunzip -c ${DIR}/CPSA90B_ACAGTG_S38_L005_R1_001.fastq.gz >${TARGET}/CPSA90B.fastq
gunzip -c ${DIR}/CPSB0B_GCCAAT_S41_L006_R1_001.fastq.gz >${TARGET}/CPSB0B.fastq
gunzip -c ${DIR}/CPSB10B_ACTTGA_S47_L006_R1_001.fastq.gz >${TARGET}/CPSB10B.fastq
gunzip -c ${DIR}/CPSB30B_GATCAG_S46_L006_R1_001.fastq.gz >${TARGET}/CPSB30B.fastq
gunzip -c ${DIR}/CPSB3B_CAGATC_S45_L006_R1_001.fastq.gz >${TARGET}/CPSB3B.fastq
gunzip -c ${DIR}/CPSB90B_TAGCTT_S44_L006_R1_001.fastq.gz >${TARGET}/CPSB90B.fastq
gunzip -c ${DIR}/CPSC0B_ATGAGC_S39_L005_R1_001.fastq.gz >${TARGET}/CPSC0B.fastq
gunzip -c ${DIR}/CPSC90B_ATTCCT_S40_L005_R1_001.fastq.gz >${TARGET}/CPSC90B.fastq
gunzip -c ${DIR}/CPSD0B_CAAAAG_S43_L006_R1_001.fastq.gz >${TARGET}/CPSD0B.fastq
gunzip -c ${DIR}/CPSD90B_CAACTA_S42_L006_R1_001.fastq.gz >${TARGET}/CPSD90B.fastq