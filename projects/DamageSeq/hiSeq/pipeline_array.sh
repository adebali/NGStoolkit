#!/bin/bash

##SBATCH --job-name=Cisplatin
#SBATCH -n 1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=./log/%A_%a.out
#SBATCH --error=./log/%A_%a.err

# #SBATCH --array=29-36,39-42
# # python pype_pipeline_nucleosome.py -n $N
# python pipeline.py -n $N

# #SBATCH --array=1-28
# N=$SLURM_ARRAY_TASK_ID
# python pipeline.py -n $N

#SBATCH --array=31,32,35,36,39,40,43,44,47-50
##SBATCH --array=1-28
##SBATCH --array=53-56
##SBATCH --array=57-60
N=$SLURM_ARRAY_TASK_ID
python pipeline.py -n $N

