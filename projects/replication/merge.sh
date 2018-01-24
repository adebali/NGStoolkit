#!/usr/bin/env bash
#
# SBATCH --mem=4000
# SBATCH --time=24:00:00
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=oadebali@gmail.com

cat XR/dataDir/merged_1kCounts.txt DS/dataDir/merged_1kCounts.txt >dataDir/merged_1kCounts.txt