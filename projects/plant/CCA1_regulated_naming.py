import os
import sys
sys.path.append('..')
from referenceGenomePath import referenceGenomePath
import slurm

paths = referenceGenomePath()
reference = paths.get('TAIR10')

code = 'sort -k1,1n -k2,2n -k3,3n genomeFiles/CAC1_regulated.bed | bedtools closest -s -iu -d -k 1 -D b -b ' +reference['transcripts'] + ' -a stdin | grep -P "\t0" | awk \'{ print $1"_"$2"_"$3"_"$6"\t"$10 }\''  
os.system(code)

# job = slurm.Slurm()