#!/bin/bash

#SBATCH --mem=8000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

# echo -e 'position\tcount\ttreatment_title\tcell\tproduct\ttreatment\treplicate\ttotalMappedReads' >dataDir/mergedDNase.txt
# cat dataDir/0106/*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.inWiDNa.paPo.coPo.fiGaWiPo.adAlMe.txt | grep -v "GM12878" | grep -vP "CPD\t48h\tB\t" | grep -vP "\tB2\t" | sed 's/\tB3\t/\tB\t/' >>dataDir/mergedDNase.txt

echo -e 'position\tcount\ttreatment_title\tcell\tproduct\ttreatment\treplicate\ttotalMappedReads' >dataDir/mergedDNase_cellNaked.txt
cat dataDir/0106/G*1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.inWiDNa.paPo.coPo.fiGaWiPo.adAlMe.txt | grep -vP "\tA0\t" | grep -vP "\tB0\t" | grep -vP "\tB2\t" | grep -vP "CPD\tB\t" | sed 's/\tC\t/\tB\t/'>>dataDir/mergedDNase_cellNaked.txt