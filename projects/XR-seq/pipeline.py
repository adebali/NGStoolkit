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
        outputDirName = '0131'
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', outputDirName))
        os.system('mkdir -p ' + self.outputDir)
        self.treatment = sampleDictionary['treatment_title']
        self.group = sampleDictionary['group']
        self.cell = sampleDictionary['cell']
        self.saveInput([self.input])
        self.bowtie_reference = '/proj/seq/data/HG19_UCSC/Sequence/BowtieIndex/genome'
        self.referenceNickname = 'hg19'
        self.fasta_reference = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.fa')
        self.chromosome_limits = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.chrSizes.bed')
        self.dnaseBedList = ['/proj/sancarlb/users/ogun/ENCODE/dnase/BJ.bed']
        self.defaultWmParams = {
            '--mem': 32000,
            '-n': 1,
            '-t': '24:00:00',
            '--job-name': 'XR-seq',
            '-o': 'log_' + self.treatment + '.txt',
            '-e': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams
        self.fifteenChromatinStates = '/proj/sancarlb/users/ogun/ENCODE/chromatinStates/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed'

        self.TSSbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min10000_TSS_win100.bed'
        self.TESbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min10000_TES_win100.bed'

    def ls_fastq2txt(self):
        codeList = ['ls && sleep 1m']
        self.execMwm(codeList)
        return self

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
        output = [self.output[0]]
        self.saveOutput([self.addExtraWord(output[0], '.' + self.referenceNickname)])
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
        nucleotideOrder = 'TCGA'
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

    def intersectWithDNase_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.dnaseBedList[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-1000,1000]
        self.execM(codeList)
        return self

    def intersectToPositions_txt2txt(self):
        codeList = [
            'cat', self.input,
            '|',
            'bedIntersect2parsedPosition.py',
            '|',
            'list2countTable.py',
            '|',
            'countTable2filledGaps.py',
            '-s', self.frame[0],
            '-e', self.frame[1],
            '|',
            'cut -f2'
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def separateStrands_bed2bed(self):
        strand = ['+', '-']
        output = [
            self.addExtraWord(self.output[0], '_Plus'), 
            self.addExtraWord(self.output[0], '_Minus')
        ]
        self.saveOutput(output)
        codeList = [
            'grep',
            strand,
            self.input[0],
            '>',
            self.output
        ]
        self.execM(codeList)
        return self

    def countChmm_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.fifteenChromatinStates,
            '-b', self.input,
            '-c',
            '-F', 0.49,
            '|',
            'cut',
            '-f', '1-4,9-10',
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countTSS_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.TSSbedFile,
            '-b', self.input,
            '-c',
            '-F', 0.5,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countTES_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.TESbedFile,
            '-b', self.input,
            '-c',
            '-F', 0.5,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def binnedCountsToPositions_bed2txt(self):
        codeList = [
            'bedBinned2totalCounts.py',
            '-i', self.input,
            '-o', self.output,
            '-n', 200,
            '-reverseStrand', '"-"'
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
    .stop()

    .branch()
        .run(p.convertBedToFasta_bed2fa, False)
        .run(p.getNucleotideAbundanceTable_fa2csv, False)
        .run(p.plotNucleotideAbundance_csv2pdf, False)
    .stop()

    .branch(False) # CHMM
        .run(p.countChmm_bed2bed, False)
    .stop()

    .branch(True)
        .run(p.separateStrands_bed2bed, False)

        .branch(False) # DNase
            .run(p.intersectWithDNase_bed2txt, False)
            .run(p.intersectToPositions_txt2txt, False)
        .stop()

        .branch(True) # TSS
            .run(p.countTSS_bed2bed, False)
            .run(p.binnedCountsToPositions_bed2txt, True)
        .stop()

        .branch(True) # TES
            .run(p.countTES_bed2bed, False)
            .run(p.binnedCountsToPositions_bed2txt, True)
        .stop()

    .stop()
)
