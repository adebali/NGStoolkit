import generalUtils
from slurm import Slurm
from glob import glob
import os

# files = glob("/proj/sancarlb/HTSF/170515_UNC13-SN749_0627_AHJM3CBCXY/RNA*.fastq.gz")
files = glob("/proj/sancarlb/HTSF/170908_UNC18-D00493_0443_BCBD4MANXX/RNA*.fastq.gz")

for file in files:
    code = "gunzip -c " + file + " >dataDir/raw/" + os.path.splitext(os.path.basename(file))[0]
    job = Slurm(code)
    job.assignParams({'--mem ': '8000'})
    job.printScript()
    job.run()
