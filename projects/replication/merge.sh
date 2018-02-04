#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

# filename='merged_1kCounts.txt'
# filename='merged_leadLag.txt'
ARRAY=(merged_leadLag_RIZ_1.0MB_c900.txt merged_leadLag_RIZ_0.5MB_c900.txt merged_leadLag_RIZ_0.1MB_c900.txt merged_leadLag_RTZ_1.0MB_c900.txt merged_leadLag_RTZ_0.5MB_c900.txt merged_leadLag_RTZ_0.1MB_c900.txt merged_leadLag_selectedRIZ_1.0MB.txt)

cd XR
python pipeline.py cat
cd ../DS
python pipeline.py cat
cd ..
for filename in ${ARRAY[*]}; do
    # echo $filename
    cat XR/dataDir/$filename  >dataDir/$filename
    tail -n +2 DS/dataDir/$filename >>dataDir/$filename
done