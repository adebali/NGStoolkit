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

# Required tools in path
# bowtie
# bedtools
# samtools

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        SAMPLE_STAT_FILE = 'samples.csv'    
        OUTPUT_DIR = "0519"
        pipe.__init__(self, input, args)
        self.input = os.path.realpath(os.path.join(os.path.curdir, 'dataDir', 'raw', input))
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        # input = generalUtils.file(os.path.realpath(os.path.join('dataDir', 'raw', input)))
        input = os.path.realpath(os.path.join('dataDir', 'raw', input))
        self.saveInput([self.input, self.input.replace("_R1_001.fastq", "_R2_001.fastq")])

        self.mm10 = {
            "name": "",
            "bowtie": "/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome",
            "fasta": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa",
            "limits": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai",
            # "genes": "/proj/sancarlb/users/ogun/seq/mm10/geneList_6.bed",
            "genes": "/proj/sancarlb/users/ogun/seq/mm10/geneList_yy_6.bed",
            "chmm": {
                "liver": "NA",
                "spleen": "NA",
                "kidney": "NA",
                "testes": "NA"
            }
        }
        self.mm9 = {
            "name": "_mm9",
            "bowtie": "/proj/seq/data/MM9_UCSC/Sequence/BowtieIndex/genome",
            "fasta": "/proj/seq/data/MM9_UCSC/Sequence/WholeGenomeFasta/genome.fa",
            "limits": "/proj/seq/data/MM9_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai",
            "genes": "NA",
            "chmm": {
                "liver": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/liver_cStates_HMM",
                "spleen": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/spleen_cStates_HMM",
                "kidney": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/kidney_cStates_HMM",
                "testes": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/testes_cStates_HMM"
            }
        }

        self.readNumber = 1000000

    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "tail --lines=+2 " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
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
            newOutputs.append(os.path.join(os.path.dirname(o),self.treatment_title + extraWord + '.' + extension))
        return newOutputs


    def bowtie_fastq2sam(self, referenceDictionary):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput([self.addExtraWord(output[0], referenceDictionary["name"])])
        noCpus = 1
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0]]
        codeList = [
            'bowtie',
            '-t', referenceDictionary["bowtie"],
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

    def fastqc_fastq2html(self):
        codeList= [
            'fastqc',
            self.input
        ]
        self.execM(codeList)
        return self

    def tophat_fastq2bam(self):
        newOutput = [self.input[0] + '_1', self.input[0] + '_2']

        codeList = [
            'ln -fs',
            self.input[0],
            newOutput[0],
            '&',
            'ln -fs',
            self.input[1],
            newOutput[1]
        ]

        self.execM(codeList)
        self.saveInput(newOutput)

        self.tophatDirectory = os.path.join(self.outputDir, 'tophat', self.treatment_title)

        codeList = [
            'mkdir -p',
            self.tophatDirectory,
            '&',
            'tophat',
            '--bowtie1',
            '-o', self.tophatDirectory,
            self.mm10['bowtie'],
            self.input[0],
            self.input[1]
        ]
        self.execM(codeList)
        self.referenceDictionary = self.mm10
        return self
   
    def convertToBed_bam2bed(self):
        self.saveInput([os.path.join(self.tophatDirectory, 'accepted_hits.bam')])
        self.output = [os.path.join(self.outputDir, self.treatment_title + '.bed')]
        self.saveOutput(self.output)
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
        self.finalBed = self.output[0]
        self.execM(codeList)
        return self

    def geneMap_bed2txt(self):
        print(self.output)
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.referenceDictionary["genes"],
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
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

    def mergeGeneCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts.txt')
        self.catFiles(wildcard, headers, output)
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
    .branch()
        .run(p.fastqc_fastq2html, False)
    .stop()

    .run(p.tophat_fastq2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)

    .branch(True)
        .run(p.geneMap_bed2txt, True)
        .run(p.normalizeCounts_txt2txt, True)
        .run(p.addTreatment_txt2txt, True)
        .cat(p.mergeGeneCounts, True)
    .stop()
)
