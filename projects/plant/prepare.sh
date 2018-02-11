#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com


HOME=/nas/longleaf/home/adebali/ogun
PROJECT_NAME=plant

TARGET=${HOME}/${PROJECT_NAME}
RAW=${TARGET}/raw

SCRIPTS=/nas/longleaf/home/adebali/ogun/scripts
PROJECT_SCRIPT=${SCRIPTS}/projects/${PROJECT_NAME}

# mkdir -p ${TARGET}
# mkdir -p ${RAW}
# mkdir -p ${TARGET}/0707
# mkdir -p ${PROJECT_SCRIPT}/dataDir
#ln -s ${TARGET} ${PROJECT_SCRIPT}/log

# DIR=/proj/sancarlb/HTSF/170630_UNC17-D00216_0455_ACB689ANXX/NERSEQ
# gunzip -c ${DIR}/ONUR1_ATCACG_S13_L003_R1_001.fastq.gz > ${RAW}/ONUR1.fastq

# DIR=/proj/sancarlb/HTSF/170725_UNC17-D00216_0460_BCBF6MANXX
# gunzip -c ${DIR}/ATCPDrep2_ATCACG_S31_L006_R1_001.fastq.gz > ${RAW}/ATCPDrep2.fastq


# DIR=/proj/sancarlb/HTSF/170925_UNC31-K00269_0083_AHKMYFBBXX/NERSEQ

# gunzip -c ${DIR}/OO_XR_cir10-_GATCAG_S22_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir10-_GATCAG.fastq
# gunzip -c ${DIR}/OO_XR_cir13-_TAGCTT_S21_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir13-_TAGCTT.fastq
# gunzip -c ${DIR}/OO_XR_cir16-_GGCTAC_S20_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir16-_GGCTAC.fastq
# gunzip -c ${DIR}/OO_XR_cir19-_CTTGTA_S19_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir19-_CTTGTA.fastq
# gunzip -c ${DIR}/OO_XR_cir1-_GCCAAT_S15_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir1-_GCCAAT.fastq
# gunzip -c ${DIR}/OO_XR_cir22-_AGTCAA_S18_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir22-_AGTCAA.fastq
# gunzip -c ${DIR}/OO_XR_cir4-_CAGATC_S17_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir4-_CAGATC.fastq
# gunzip -c ${DIR}/OO_XR_cir7-_ACTTGA_S16_L003_R1_001.fastq.gz >${RAW}/OO_XR_cir7-_ACTTGA.fastq

# DIR=/proj/sancarlb/HTSF/171208_UNC13-SN749_0659_AH3NTGBCX2

# gunzip -c ${DIR}/OO_XR_cir2-10_GCCAAT_S2_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-10_GCCAAT.fastq
# gunzip -c ${DIR}/OO_XR_cir2-13_CAGATC_S3_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-13_CAGATC.fastq
# gunzip -c ${DIR}/OO_XR_cir2-16_ACTTGA_S4_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-16_ACTTGA.fastq
# gunzip -c ${DIR}/OO_XR_cir2-19_GATCAG_S5_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-19_GATCAG.fastq
# gunzip -c ${DIR}/OO_XR_cir2-1_TTAGGC_S1_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-1_TTAGGC.fastq
# gunzip -c ${DIR}/OO_XR_cir2-22_TAGCTT_S6_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-22_TAGCTT.fastq
# gunzip -c ${DIR}/OO_XR_cir2-4_TGACCA_S7_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-4_TGACCA.fastq
# gunzip -c ${DIR}/OO_XR_cir2-7_ACAGTG_S8_L001_R1_001.fastq.gz >${RAW}/OO_XR_cir2-7_ACAGTG.fastq

DIR=/proj/sancarlb/HTSF/180207_UNC32-K00270_0079_BHN5CVBBXX

gunzip -c ${DIR}/Onur_6-4_A_GTGGCC_S8_L002_R1_001.fastq.gz >${RAW}/Onur_6-4_A_GTGGCC.fastq
