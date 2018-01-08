!/usr/bin/env bash

SBATCH --mem=4000
SBATCH --time=24:00:00
SBATCH --mail-type=END,FAIL
SBATCH --mail-user=oadebali@gmail.com

SOURCE=/proj/sancarlb/users/ogun/sancarlabutils/projects/plant/dataDir/raw
SOURCE2=/proj/sancarlb/users/ogun/sancarlabutils/projects/plant/dataDir/0707
TARGET=/pine/scr/a/d/adebali/GEO/oadebali@gmail.com

cp $SOURCE/OO*fastq $TARGET/
cp $SOURCE2/*.ge.al.toBa.soBa.toBe.geFa.geTT.so.su.suBa.geCo_*.bdg $TARGET/

cp $SOURCE2/*ge.al.toBa.soBa.toBe.geFa.geTT.so.bed $TARGET/

cd $TARGET


# for f in *.ge.al.toBa.soBa.toBe.geFa.geTT.so.su.suBa.geCo_*.bdg; do 
# mv -- "$f" "${f%.bdg}.txt"
# done

md5sum * >md5sum.txt
