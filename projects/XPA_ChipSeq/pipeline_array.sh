#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=XPApipe
#SBATCH -n 1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=dataDir/log/pipeline.%A_%a.out
#SBATCH --error=dataDir/log/pipeline.%A_%a.err

python pipeline.py -n $SLURM_ARRAY_TASK_ID