import os

dir = "/netscr/adebali/GEO/oadebali@gmail.com/"

fileList = [
    'WT04_TGACC_L001_R1_001.fastq',
    'STLWT5_ACAGTG_L008_R1_001.fastq',
    'PHR11_GGCTAC_L007_R1_001.fastq',
    'UVRD10_TAGCTT_L007_R1_001.fastq',
    'MFD09_GATCAG_L007_R1_001.fastq',
    'STLWT5_ACAGTG_L008_R1_001.fastq',
    'STLPHR8_ACTTGA_L008_R1_001.fastq',
    'STLMFD6_GCCAAT_L008_R1_001.fastq',
    'STLUVRD7_CAGATC_L008_R1_001.fastq'
]

for file in fileList:
    print(file)
    # os.system("cp "+ "dataDir/" + file + " " + dir)
    # os.system("gzip " + dir + "*")
    os.system("cp "+ "dataDir/" + file.replace('.fastq','') + ".cutadapt13.2900K.bowtieNotall.sorted.bed.TT.*.gene.Cov.bed" + " " + dir)

# os.system("md5sum " + dir + "*.gz >GEO_md5sum.txt")
os.system("md5sum " + dir + "*.bed >GEO_md5sum_bed.txt")

