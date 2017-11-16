import os
import sys
import slurm

# filenames = ['WRKY40-SRX2039034']
filenames = []
finalIndividualJobs = []
for line in open('plantdb_AtList.txt'):
    if not line.strip().startswith("#"):
        filenames.append(line.rstrip())

os.chdir('/nas/longleaf/home/adebali/ogun/seq/TAIR10/Liu2017/epigeneticMarks')
os.system('mkdir -p log')
# os.system('rm *.bed')

for filename in filenames:
    wmParams = {
    '--mem=': 8000,
    '-n ': 1,
    '-t ': '24:00:00',
    '--job-name=': 'plantdb',
    '-e ': 'log/err_plantdb' + filename + '.txt',
    '-o ': 'log/out_plantdb' + filename + '.txt'
    }

    code = 'wget http://systemsbiology.cau.edu.cn/chromstates/At_bwfile/' + filename + '.bw'
    job1 = slurm.Slurm(code)
    job1.assignParams(wmParams)
    # job1_id = job1.run()
    
    code = 'bigWigToBedGraph ' + filename + '.bw ' + filename + '.bdg'
    job2 = slurm.Slurm(code)
    job2.assignParams(wmParams)
    # job2.setDependencies([job1_id])
    # job2_id = job2.run()

    code = 'macs2 bdgpeakcall -i ' + filename + '.bdg -o ' + filename + '.txt'
    job3 = slurm.Slurm(code)
    job3.assignParams(wmParams)
    # job3.setDependencies([job2_id])
    # job3_id = job3.run()

    code = 'tail -n +2 ' + filename + '.txt | sed \'s/chr//g\' | grep -vP "^C" | grep -vP "^M" | awk \'{print $1"\\t"$2"\\t"$3"\\t' + filename + '\\t"$5"\\t"$6}\' | bed2removeChromosomeEdges.py -g ../../genome.fa.fai --fixed -l 1000 | sort -k1,1 -k2,2n -k3,3n > ' + filename + '.bed'
    job4 = slurm.Slurm(code)
    job4.assignParams(wmParams)
    # job4.setDependencies([job3_id])
    # job4_id = job4.run()

    # finalIndividualJobs.append(job4_id)

code = 'cat *.bed | sort -k1,1 -k2,2n -k3,3n >epigeneticMarkers.bed'
catJob = slurm.Slurm(code)
catJob.assignParams(wmParams)
# catJob.setDependencies(finalIndividualJobs)
# catJobId = catJob.run()

code = 'bedtools shuffle -i epigeneticMarkers.bed -g ../../genome.fa.fai | bed2removeChromosomeEdges.py -g ../../genome.fa.fai --fixed -l 1000 | sort -k1,1 -k2,2n -k3,3n >epigeneticMarkers_shuffled.bed'
shuffleJob = slurm.Slurm(code)
shuffleJob.assignParams(wmParams)
# shuffleJob.setDependencies([catJobId])
shuffleJobId = shuffleJob.run()