#!/bin/bash

#SBATCH --mem=8000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

# echo -e 'position\tcount\tTF\tproduct\ttreatment\treplicate\tmappedReadNumber\tstrand' >dataDir/mergedTFmotif.txt
#cat dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.inWiTFm*.inReNeToPo.adTr.ad*.txt | grep -vP "_6-4\tnakedDNA\tB\t" | grep -vP "CPD\tnakedDNA\tB\t" | grep -vP "\tA0\t" | grep -vP "\tB0\t" | grep -vP "\tB2\t" | sed 's/\tC\t/\tB\t/' >>dataDir/mergedTFmotif.txt

echo -e 'position\tcount\tTF\tproduct\ttreatment\treplicate\tmappedReadNumber\tstrand' >dataDir/mergedFakeMotif.txt
cat dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.inWiFaMo*.inReNeToPo.adTr.ad*.txt | grep -vP "_6-4\tnakedDNA\tB\t" | grep -vP "CPD\tnakedDNA\tB\t" | grep -vP "\tA0\t" | grep -vP "\tB0\t" | grep -vP "\tB2\t" | sed 's/\tC\t/\tB\t/' >>dataDir/mergedFakeMotif.txt




