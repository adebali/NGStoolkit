import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob
import argparse

# Required tools in path
# bowtie
# bedtools
# samtools
# subsample

SAMPLE_STAT_FILE = 'dataDir/samples.csv'

parser = argparse.ArgumentParser(description='XR-seq Pipeline')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-n', required= False, help='output')
parser.add_argument('--mock', required= False, action='store_true', help='mock flag')
parser.add_argument('--noPrint', required= False, action='store_true', help='prints no code when stated')
args = parser.parse_args()
inputFile = args.i
inputIndex = args.n


if inputFile == None and inputIndex == None:
    raise ValueError("No input or index was given! Exiting ...")
    sys.exit()
elif inputFile != None:
    fileBaseName = os.path.basename(inputFile)
    samples = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'sample')
    sampleDictionary = samples[fileBaseName][0]
    input = generalUtils.file(inputFile)
else:
    indexes = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'no')
    sampleDictionary = indexes[inputIndex][0]
    input = generalUtils.file(os.path.realpath(os.path.join('dataDir', 'raw', sampleDictionary['sample'])))


class myPipe(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        self.execM = self.execMwm
        outputDirName = '0220'
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', outputDirName))
        os.system('mkdir -p ' + self.outputDir)
        self.treatment = sampleDictionary['treatment_title']
        self.group = sampleDictionary['group']
        self.saveInput([self.input])
        self.bowtie_reference = '/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome'
        self.fasta_reference = generalUtils.file('/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa')
        self.chromosome_limits = generalUtils.file('/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai')
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.treatment + '.txt',
            '-e ': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams

    def cutadapt_fastq2fastq(self):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)
        adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'
        codeList = [
            'cutadapt',
            '-a', adapter,
            '-o', output,
            input
        ]
        self.execM(codeList)
        return self


    def bowtie_fastq2sam(self):
        noCpus = 1
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0]]
        codeList = [
            'bowtie',
            '-t', self.bowtie_reference,
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            '-n', 2, # No more than 2 mismatches
            '-e', 70, # The sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed E 
            '-m 4', # Do not report the reads that are mapped on to more than 4 genomic locations
            '-p', noCpus,
            '--seed 123', # Randomization parameter in bowtie,
            self.input,
            self.output
        ]
        self.execM(codeList)
        return self
   
    def convertToBam_sam2bam(self):
        codeList = [
            'samtools',
            'view',
            '-Sb'
            '-o',
            self.output,
            self.input
        ]
        self.execM(codeList)
        return self
   
    def convertToBed_bam2bed(self):
        codeList = [
            'bedtools',
            'bamtobed',
            '-i', self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def uniqueSort_bed2bed(self):
        codeList = [
            'sort',
            '-u',
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertBedToFasta_bed2fa(self):
        codeList = [
            'bedtools',
            'getfasta',
            '-fi', self.fasta_reference,
            '-bed', self.input,
            '-fo', self.output,
            '-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
        ]
        self.execM(codeList)
        return self

    def getNucleotideAbundanceTable_fa2csv(self):
        nucleotideOrder = 'GCTA'
        codeList = [
            'fa2nucleotideAbundanceTable.py',
            '-i', self.input,
            '-o', self.output,
            '-n', nucleotideOrder,
            '--percentage'
        ]
        self.execM(codeList)
        return self

    def plotNucleotideAbundance_csv2pdf(self):
        codeList = [
            'plotNucleotideAbundance.r',
            self.input,
            self.treatment
        ]
        self.execM(codeList)
        return self

    def lengthDistribution_bed2csv(self):
        codeList = [
            'bed2lengthDistribution.py',
            '-i', self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def get26mer_bed2bed(self):
        codeList = [
            'bed2getCertainIntervalLengths.py',
            '-i', self.input,
            '-o', self.output,
            '-l', 27
        ]
        self.execM(codeList)
        return self

    def plotLengthDistribution_csv2pdf(self):
        codeList = [
            'plotLengthDistribution.R',
            self.input,
            self.treatment
        ]
        self.execM(codeList)
        return self

###########################################################
#  Pipeline
###########################################################
p = myPipe(input)
(p
    .run(p.cutadapt_fastq2fastq, False)
    .run(p.bowtie_fastq2sam, False)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)
    
    .branch()
        .run(p.lengthDistribution_bed2csv, False)
        .run(p.plotLengthDistribution_csv2pdf, False)
    .stop()

    .branch()
        .run(p.get26mer_bed2bed, False)
        .run(p.convertBedToFasta_bed2fa, False)
        .run(p.getNucleotideAbundanceTable_fa2csv, False)
        .run(p.plotNucleotideAbundance_csv2pdf, True)
    .stop()
)
