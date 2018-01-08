#!/usr/bin/env python
import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob
import argparse
import argument
sys.path.append('../..')
from referenceGenomePath import referenceGenomePath

# Required tools in path
# bowtie
# bedtools
# samtools

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        SAMPLE_STAT_FILE = 'samples.csv'    
        OUTPUT_DIR = "1801"
        pipe.__init__(self, input, args)
        self.saveInput(os.path.realpath(os.path.join(os.path.curdir, 'dataDir', OUTPUT_DIR, input)))
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.paths = referenceGenomePath()
        self.reference = self.paths.get('S288C_R64_2_1')


    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        os.system(code)

    def fullPath2wildcard(self, fullPath):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join(["*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def prettyOutput(self):
        newOutputs = []
        for o in self.output:
            if 'Plus' in o:
                extraWord = '_Plus'
            elif 'Minus' in o:
                extraWord = '_Minus'
            else:
                extraWord = ''
            extension = pipeTools.getExtension(o)
            newOutputs.append(os.path.join(os.path.dirname(o),self.title + extraWord + '.' + extension))
        return newOutputs

    def cutadapt_fastq2fastq(self):
        adapter = 'ATCTCGTATGCCGTCTTCTGCTTG'
        codeList = [
            'cutadapt',
            '-a', adapter,
            '-o', self.output,
            self.input
        ]
        self.execM(codeList)
        return self

    def bowtie_fastq2sam(self):
        noCpus = 8
        self.mutateWmParams({'-n ': str(noCpus)})
        codeList = [
            'bowtie',
            '-t', self.reference['bowtie'],
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            # '-n', 2, # No more than 2 mismatches
            '-e', 70, # The sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed E 
            # '-m 4', # Do not report the reads that are mapped on to more than 4 genomic locations
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

    def fastqc_fastq2html(self):
        
        codeList= [
            'fastqc',
            self.input
        ]
        self.execM(codeList)
        return self

    def tophat_fastq2bam(self):
        newOutput = [self.input + '_1']

        codeList = [
            'ln -fs',
            self.input,
            newOutput
        ]

        self.execM(codeList)
        self.saveInput(newOutput)

        self.tophatDirectory = os.path.join(self.outputDir, 'tophat', self.title)

        codeList = [
            'mkdir -p',
            self.tophatDirectory,
            '&',
            'tophat',
            '--bowtie1',
            '-o', self.tophatDirectory,
            self.reference['bowtie'],
            self.input,
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
        self.finalBed = self.output
        self.execM(codeList)
        return self

    def geneMap_bed2txt(self):
        # newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        # self.saveOutput(newOutput)
        # strandParameters = ['-S', '-s'] #[different strand, same strand]
        # self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["genes"],
            '-b', self.input,
            '-wa',
            '-c',
            # strandParameters,
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def normalizeCounts_txt2txt(self):
        if self.runFlag and self.runMode:
            self.scaleFactor = float(1000000)/self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number')
        else:
            self.scaleFactor = 1
        codeList = [
            'bedCount2normalizedCount.py',
            '-i', self.input,
            '-c', 7,
            '-m', self.scaleFactor,
            '-l', 1000,
            # '--bypassLength',
            '-o', self.output
        ]
        self.execM(codeList)
        return self
    
    def addTreatment_txt2txt(self):
        columns = self.list2attributes(self.attributes)
        # columnStringList = [' '.join(columns + ['TS']), ' '.join(columns + ['NTS'])]
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', ' '.join(columns)
        ]
        self.execM(codeList)
        return self

    def writeTotalMappedReads_bed2txt(self):
        codeList = [
            'grep -c "^"',
            self.input,
            '>',
            self.output
        ]
        self.totalMappedReadsFile = self.output[0]
        self.execM(codeList)
        return self

    def mergeGeneCounts(self, extraWord = ""):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        os.system("Rscript dcast.r " + output + ' ' + output.replace('.txt', '_dcasted.txt'))
        return self
   
    def toBg_bed2bg(self):
        try:
            totalMappedReads = int(open(self.totalMappedReadsFile, 'r').readline())
        except:
            if self.runFlag == True and self.runMode == True:
                raise ValueError('no read count file')
            totalMappedReads = 1000000
        scaleFactor = float(1000000) / totalMappedReads
        # if self.runFlag and self.runMode:
        #     scaleFactor = float(1000000)/self.internalRun(bed.bed(self.input[0]).getHitNum, [], self.runFlag, 'get hit number')
        # else:
        #     scaleFactor = 1
        codeList = [
            'bedtools',
            'genomecov',
            '-i', self.input,
            '-g', self.reference["limits"],
            '-bg',
            '-scale', scaleFactor,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def toBw_bg2bw(self):
        codeList = [
            'bedGraphToBigWig',
            self.input,
            self.reference["limits"],
            self.output
        ]
        self.execM(codeList)
        return self

    def makeScoreBed6_txt2bed(self):
        codeList = [
            'awk \'{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$7"\\t"$6}\'',
            self.input,
            '>',
            self.output
        ]
        self.execM(codeList)
        return self

def getArgs():
    parser = argparse.ArgumentParser(description='XR-seq ZT Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-n', required= True, help='output')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-n', required= False, default="1", help='output')
    parser_cat.add_argument('--mock', required= False, default=True, action='store_true', help='mock flag')
    parser_cat.add_argument('--noPrint', required= False, default=True, action='store_true', help='prints no code when stated')

    args = parser.parse_args()
    return argument.args(args)

def sampleIO(fileName, in_, by_, out_):
    d1 = generalUtils.table2dictionary(generalUtils.file(fileName), by_)
    d2 = d1[in_][0]
    return d2[out_]

def getInputFromIndex(n):
    SAMPLE_STAT_FILE = 'samples.csv'
    return sampleIO(SAMPLE_STAT_FILE, n, 'no', 'sample')


args = getArgs()
inputIndex = args.get("n")
input = getInputFromIndex(inputIndex)
###########################################################
#  Pipeline
###########################################################
p = myPipe(input, args)
(p
    .branch(False)
        .run(p.fastqc_fastq2html, False)
    .stop()

    .run(p.cutadapt_fastq2fastq, False)
    .run(p.bowtie_fastq2sam, False)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)

    .branch(True)
        .run(p.writeTotalMappedReads_bed2txt, False)
    .stop()

    .branch(True)
        .run(p.geneMap_bed2txt, False)
        .run(p.normalizeCounts_txt2txt, False)
        .run(p.addTreatment_txt2txt, False)
        .branch(True)
            .run(p.makeScoreBed6_txt2bed, True)
        .stop()
        .cat(p.mergeGeneCounts, False)
    .stop()

    .run(p.toBg_bed2bg, False)
    .run(p.toBw_bg2bw, False)

    .stop()
)
