#!/bin/bash
##SBATCH --array=1-28
##SBATCH --array=36-39
##SBATCH --array=33
##SBATCH --job-name=NHF1
##SBATCH --job-name=Cisplatin
##SBATCH -n 8
#SBATCH -n 1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=./log/pipeline.%A_%a.out
#SBATCH --error=./log/pipeline.%A_%a.err

#SBATCH --array=29-36,39-42
N=$SLURM_ARRAY_TASK_ID
# python pype_pipeline_nucleosome.py -n $N

python pipeline.py -n $N

