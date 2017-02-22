import generalUtils
from slurm import Slurm

samples = generalUtils.table2dictionary('dataDir/XRseq_SRRlist.csv', 'sample')
for sample in samples.keys():
    code = ''
    for subsample in samples[sample]:
        # print(subsample)
        code += '/nas/longleaf/apps/sratoolkit/2.8.0/sra-tools/bin/fastq-dump ' + subsample['accession_no'] + ' -O ' + 'dataDir/raw/' + subsample['accession_no'] + '.fastq && '

    code += ' cat '
    for subsample in samples[sample]:
        code += 'dataDir/raw/' + subsample['accession_no'] + '.fastq/' + subsample['accession_no'] + '.fastq '
    code += '>dataDir/raw/' + sample + '.fastq'
    print(code)
    job = Slurm(code,8)
    # job.printScript()
    job.run()