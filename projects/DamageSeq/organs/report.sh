#!/bin/bash

##SBATCH --job-name=report
#SBATCH -n 1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=./log/%A_%a.out
#SBATCH --error=./log/%A_%a.err

python pipeline_array.py -g 0 --report