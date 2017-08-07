#!/bin/bash

#SBATCH --mem=8000
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oadebali@gmail.com

echo -e 'position\tcount\ttreatment_title\tcell\tproduct\ttreatment\treplicate\ttotalMappedReads' >dataDir/mergedDNase.txt
cat dataDir/0131/*.cu.bo.hg19.coToBa.coToBe.unSo.inWiDNa.paPo.coPo.fiGaWiPo.adAlMe.txt >>dataDir/mergedDNase.txt

