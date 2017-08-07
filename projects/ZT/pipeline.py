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
        OUTPUT_DIR = "0512"
        pipe.__init__(self, input, args)
        self.input = os.path.realpath(os.path.join(os.path.curdir, 'dataDir', 'raw', input))
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        input = generalUtils.file(os.path.realpath(os.path.join('dataDir', 'raw', input)))
        self.saveInput([self.input])
        
        self.mm10 = {
            "name": "",
            "bowtie": "/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome",
            "fasta": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa",
            "limits": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai",
            # "genes": "/proj/sancarlb/users/ogun/seq/mm10/geneList_6.bed",
            "genes": "dataDir/mm10_geneList_yy_6.bed",
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

    def addTreatmentAndStrand_txt2txt(self):
        columns = self.list2attributes(self.attributes)
        columnStringList = [' '.join(columns + ['TS']), ' '.join(columns + ['NTS'])]
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', columnStringList
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
        self.finalBed = self.output[0] 
        self.execM(codeList)
        return self

    def countGene_bed2txt(self, referenceDictionary):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', referenceDictionary["genes"],
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

    def mergeGeneCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S")
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes + ['TSNTS']
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def countChmm_bed2bed(self, reference):
        chmmFile = reference["chmm"][self.organ]
        codeList = [
            'bedtools',
            'intersect',
            '-a', chmmFile,
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
    .run(p.bowtie_fastq2sam, False, p.mm10)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)

    .branch(True)
        .run(p.countGene_bed2txt, False, p.mm10)
        .run(p.normalizeCounts_txt2txt, True)
        
        .branch()
            .run(p.addTreatmentAndStrand_txt2txt, True)
            .cat(p.mergeGeneCounts, True)
        .stop()
    .stop()
)

# Chromatin State Analysis on MM9
p = myPipe(input, args)
(p
    .run(p.bowtie_fastq2sam, False, p.mm9)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)
    .run(p.countChmm_bed2bed, False, p.mm9)
)

# def f(x):
#     return (2e-17) * x**6 - (1e-13) * x**5 + (5e-10)*x**4 - (8e-07)*x**3 + (0.0008)*x**2 + (0.8324)*x + 20.877
# # y = 2E-17x6 - 1E-13x5 + 5E-10x4 - 8E-07x3 + 0.0008x2 + 0.8324x + 20.877
# y = 1.8501x - 542.29
# 240.49e0.0013x
