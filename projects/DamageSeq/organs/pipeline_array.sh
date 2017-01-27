#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=organ
#SBATCH -n 8
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=./log/%A_%a.out
#SBATCH --error=./log/%A_%a.err

N=$SLURM_ARRAY_TASK_ID
python pipeline.py -n $N

