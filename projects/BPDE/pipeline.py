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

SAMPLE_STAT_FILE = 'samples.csv'

parser = argparse.ArgumentParser(description='XR-seq Pipeline for BPDE')
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
    input = inputFile
else:
    indexes = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'no')
    sampleDictionary = indexes[inputIndex][0]
    input = 'dataDir/0308/' + sampleDictionary['sample'] + '.bed'


class myPipe(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        # self.execM = self.execMwm
        outputDirName = '0308'
        self.outputDir = os.path.dirname(self.input)
        os.system('mkdir -p ' + self.outputDir)
        self.treatment = sampleDictionary['treatment_title']
        self.group = sampleDictionary['group']
        self.version = sampleDictionary['version']
        self.cell = sampleDictionary['cell']
        self.saveInput([self.input])
        self.bowtie_reference = '/proj/sancarlb/users/wentao/HG19_UCSC/Bt1Male/male_bt1.hg19'
        self.referenceNickname = 'hg19'
        self.fasta_reference = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.fa')
        self.chromosome_limits = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.chrSizes.bed')
        self.dnaseBedList = ['/proj/sancarlb/users/ogun/ENCODE/dnase/BJ.bed']
        self.defaultWmParams = {
            '--mem': 32000,
            '-n': 1,
            '-t': '5-00:00:00',
            '--job-name': 'BPDE_XR-seq',
            '-o': 'log_' + self.treatment + '.txt',
            '-e': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams
        if self.cell == "GM12878":
            self.fifteenChromatinStates = '/proj/sancarlb/users/ogun/ENCODE/chromatinStates/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed'
        elif self.cell == "NHF1":
            self.fifteenChromatinStates = '/proj/sancarlb/users/ogun/ENCODE/chromatinStates/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed'

        self.TSSbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min15000_TSS_win100.bed'
        self.TESbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min15000_TES_win100.bed'


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
    .run(p.uniqueSort_bed2bed, False)

    .branch(True) # CHMM
        .run(p.countChmm_bed2bed, True)
    .stop()

    .branch(False)
        .run(p.separateStrands_bed2bed, True)

        .branch(True) # TSS
            .run(p.countTSS_bed2bed, True)
            .run(p.binnedCountsToPositions_bed2txt, True)
        .stop()

        .branch(True) # TES
            .run(p.countTES_bed2bed, True)
            .run(p.binnedCountsToPositions_bed2txt, True)
        .stop()

    .stop()
)
