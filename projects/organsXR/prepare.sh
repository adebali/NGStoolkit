#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

TARGET=/nas/longleaf/home/adebali/ogun/organsXR/raw

# DIR=/proj/sancarlb/HTSF/170215_UNC17-D00216_0422_AC8BE5ANXX
# gunzip -c ${DIR}/KIDNEY-CISPL_ACTTGA_S19_L005_R1_001.fastq.gz > ${TARGET}/KIDNEY-CISPL_ACTTGA.fastq
# gunzip -c ${DIR}/LIVERA_AGTTCC_S21_L005_R1_001.fastq.gz > ${TARGET}/LIVERA_AGTTCC.fastq
# gunzip -c ${DIR}/LIVER-CISPL_GCCAAT_S20_L005_R1_001.fastq.gz > ${TARGET}/LIVER-CISPL_GCCAAT.fastq
# gunzip -c ${DIR}/LUNG-CISPL_CAGATC_S22_L005_R1_001.fastq.gz > ${TARGET}/LUNG-CISPL_CAGATC.fastq
# gunzip -c ${DIR}/SPLEEN-CISPL_GATCAG_S23_L005_R1_001.fastq.gz > ${TARGET}/SPLEEN-CISPL_GATCAG.fastq

# DIR=/proj/sancarlb/HTSF/170426_UNC13-SN749_0624_AH77GVBCXY
# gunzip -c ${DIR}/AYKIDNEY-CISPL_CAGATC_S37_L002_R1_001.fastq.gz > ${TARGET}/AYKIDNEY-CISPL_CAGATC.fastq
# gunzip -c ${DIR}/AYLIVER-CISPL_GCCAAT_S36_L002_R1_001.fastq.gz > ${TARGET}/AYLIVER-CISPL_GCCAAT.fastq
# gunzip -c ${DIR}/AYLUNG-CISPL_ACTTGA_S35_L002_R1_001.fastq.gz > ${TARGET}/AYLUNG-CISPL_ACTTGA.fastq
# gunzip -c ${DIR}/AYSPLEEN-CISPL_GATCAG_S34_L002_R1_001.fastq.gz > ${TARGET}/AYSPLEEN-CISPL_GATCAG.fastq

# DIR=/proj/sancarlb/HTSF/170620_UNC31-K00269_0067_AHKFV7BBXX
# gunzip -c ${DIR}/XRKidney_A-_GGCTAC_S26_L003_R1_001.fastq.gz > ${TARGET}/XRKidney_A_GGCTAC.fastq
# gunzip -c ${DIR}/XRKidney_B-_AGTTCC_S25_L003_R1_001.fastq.gz > ${TARGET}/XRKidney_B_AGTTCC.fastq
# gunzip -c ${DIR}/XRLung_A-_CTTGTA_S24_L003_R1_001.fastq.gz > ${TARGET}/XRLung_A_CTTGTA.fastq
# gunzip -c ${DIR}/XRLung_B-_ATGTCA_S29_L003_R1_001.fastq.gz > ${TARGET}/XRLung_B_ATGTCA.fastq
# gunzip -c ${DIR}/XRSpleen_A-_AGTCAA_S31_L003_R1_001.fastq.gz > ${TARGET}/XRSpleen_A_AGTCAA.fastq
# gunzip -c ${DIR}/XRSpleen_B-_CCGTCC_S30_L003_R1_001.fastq.gz > ${TARGET}/XRSpleen_B_CCGTCC.fastq
 

DIR=/proj/sancarlb/HTSF/170925_UNC31-K00269_0083_AHKMYFBBXX/NERSEQ
# gunzip -c ${DIR}/XRKidney_C-_GGCTAC_S10_L002_R1_001.fastq.gz > ${TARGET}/XRKidney_C-_GGCTAC.fastq
# gunzip -c ${DIR}/XRKidney_D-_CCGTCC_S9_L002_R1_001.fastq.gz > ${TARGET}/XRKidney_D-_CCGTCC.fastq
# gunzip -c ${DIR}/XRLung_C-_CTTGTA_S8_L002_R1_001.fastq.gz > ${TARGET}/XRLung_C-_CTTGTA.fastq
# gunzip -c ${DIR}/XRLung_D-_GTAGAG_S6_L002_R1_001.fastq.gz > ${TARGET}/XRLung_D-_GTAGAG.fastq
# gunzip -c ${DIR}/XRSpleen_C-_AGTCAA_S7_L002_R1_001.fastq.gz > ${TARGET}/XRSpleen_C-_AGTCAA.fastq
# gunzip -c ${DIR}/XRSpleen_D-_GTCCGC_S11_L002_R1_001.fastq.gz > ${TARGET}/XRSpleen_D-_GTCCGC.fastq
# gunzip -c ${DIR}/HaCat_A-_AGTTCC_S12_L002_R1_001.fastq.gz > ${TARGET}/HaCat_A-_AGTTCC.fastq
# gunzip -c ${DIR}/HaCat_B-_ATGTCA_S13_L002_R1_001.fastq.gz > ${TARGET}/HaCat_B-_ATGTCA.fastq
