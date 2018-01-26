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
import tempfile

# Required tools in path
# bowtie
# bedtools
# samtools

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        SAMPLE_STAT_FILE = 'samples.csv'    
        OUTPUT_DIR = "1801"
        pipe.__init__(self, input, args)
        self.saveInput(os.path.realpath(os.path.join(os.path.curdir, 'dataDir', 'raw', input)))
        self.input = [self.input]
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input[0]), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'SRA_id')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.paths = referenceGenomePath()
        self.reference = self.paths.get('S288C_R64_2_1')


    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        os.system(code)

    def fullPath2wildcard(self, fullPath, prefix = ''):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join([prefix + "*"] + basename.split(".")[1:])
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

    def dump_2fastq(self):
        output = [self.SRA_id + '.fastq']
        self.saveOutput([os.path.join(os.path.realpath(self.outputDir), os.path.basename(output[0]))])
        print(self.output)
        codeList = [
            'fastq-dump',
            self.SRA_id,
            '-O', os.path.realpath(self.outputDir)
        ]
        self.execM(codeList)
        return self

    def cutadapt_fastq2fastq(self):
        # self.output = [self.output]
        # self.saveOutput([os.path.join(os.path.realpath(self.outputDir), os.path.basename(self.output[0]))])
        adapters = [
            'ATCACCGACTGCCCATAGAGAGGC',
            'ATCCTCTTCTGAGTCGGAGACACGCAGGGATGAGATGGC',
            'CCATCTCATCCCTGCGTGTCTCCGACTCAGAAGAGGATNNNNNN',
            'ATCACGAACTGAGTCGGAGACACGCAGGGATGAGATGGC',
            'CCATCTCATCCCTGCGTGTCTCCGACTCAGTTCGTGATNNNNNN',
            'ATCTCAGGCTGAGTCGGAGACACGCAGGGATGAGATGGC',
            'CCATCTCATCCCTGCGTGTCTCCGACTCAGCCTGAGATNNNNNN',
            'CCATCTCATCCCTGCGTGTCTCCGAC',
            'CCTCTCTATGGGCAGTCGGTGATT'
        ]
        codeList = [
            'cutadapt',
            '-a', ' -a '.join(adapters),
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
        if len(self.output) > 1:
            raise ValueError('Final bed should be a single file')
        self.finalBed = self.output[0]
        self.execM(codeList)
        return self

    def fastqc_fastq2html(self):
        
        codeList= [
            'fastqc',
            self.input
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
        if len(self.output) > 1:
            raise ValueError('Final bed should be a single file')
        self.finalBed = self.output[0]
        self.execM(codeList)
        return self

    def sort_bed2bed(self):
        codeList = [
            'sort',
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        if len(self.output) > 1:
            raise ValueError('Final bed should be a single file')
        self.finalBed = self.output[0]
        self.execM(codeList)
        return self

    def geneMap_bed2txt(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-s', '-S'] #[same strand, different strand] --> Opposite of the damage-seq
        self.input = self.input * 2
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["genes"],
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

    def addTreatment_txt2txt(self, *nargs):
        columns = self.list2attributes(self.attributes) + list(nargs)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', ' '.join(columns)
        ]
        self.execM(codeList)
        columnHeaders = self.attributes
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

    def convertBedToFasta_bed2fa(self):
        codeList = [
            'bedtools',
            'getfasta',
            '-fi', self.reference['fasta'],
            '-bed', self.input,
            '-fo', self.output,
            '-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
        ]
        self.execM(codeList)
        return self
    
    def getNucleotideAbundanceTable_fa2csv(self, kmer=1):
        # nucleotideOrder = 'GCTA'
        headerList = []
        columns = self.list2attributes(self.attributes)
        for i in range(len(columns)):
            headerList.append(sorted(self.attributes)[i] + ',' + columns[i])
        codeList = [
            'fa2kmerAbundanceMeltedData.py',
            '-i', self.input,
            '-o', self.output,
            '-k', kmer,
            '--percentage'
        ]
        self.execM(codeList)
        return self

    def addTreatment_csv2txt(self):
        columns = self.list2attributes(self.attributes)
        columnStringList = ' '.join(columns)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', columnStringList
        ]
        self.execM(codeList)
        return self

    def mergeNucleotideAbundance(self, extraWord = ""):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['position', 'sequence', 'value'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_nucAbu' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self   

    def getDamageSites_fa2bed(self):
        self.motifRegex = '\'.{4}(A|a|G|g)(A|a|G|g).{4}\'' # Get GG only at the positions 5 and 6.

        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', self.motifRegex 
        ]
        self.execM(codeList)
        return self

    def getDimerAbundanceTable_fa2csv(self):
        codeList = [
            'fa2kmerAbundanceMeltedData.py',
            '-i', self.input,
            '-o', self.output,
            '-k', 2,
            '--percentage'
        ]
        self.execM(codeList)
        return self

    def slopBed_bed2bed(self):
        slopB = 6
        self.addWordToOutput('b6')
        codeList = [
            'bedtools',
            'slop',
            '-i', self.input,
            '-g', self.reference['limits'], # Chomosomal lengths, needed in case region shift goes beyond the chromosomal limits.
            '-b', slopB,
            '-s', # Take strand into account
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertToFixedRange_bed2bed(self):
        side = 'left'
        sideInput = side[0]
        fixedRange = 10
        self.addWordToOutput(str(fixedRange))

        codeList = [
            'bed2fixedRangeBed.py',
            '-i', self.input,
            '-s', sideInput,
            '-l', fixedRange,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def addTSNTS_txt2txt(self, *sites):
        columns = self.list2attributes(self.attributes)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', [ ' '.join(columns) +  ' TS ' + ' '.join(sites), 
                    ' '.join(columns) +  ' NTS ' + ' '.join(sites)]
        ]
        self.execM(codeList)
        columnHeaders = self.attributes
        return self

    def transcriptIntersect_bed2txt(self, args={}):
        # transcripts = self.reference[args.get('genes', 'genes')]
        transcripts = self.reference[args.get('genes', 'scoredGenes')]
        keyword = args.get('keyword', '')
        output = [self.addExtraWord(self.output[0], keyword)]
        self.saveOutput(output)
        
        sliceFlag = args.get('slice', False)
        if sliceFlag:
            numberOfSlices = args.get('n')
            currentSlice = args.get('sliceNum')
            keyword = args.get('keyword')
            sliceTemp_, sliceTemp = tempfile.mkstemp()
            codeList = [
                'bed2sliceByScore.py',
                '-i', transcripts,
                '-n', numberOfSlices,
                '-slice', currentSlice,
                '>', sliceTemp
            ]
            self.execM(codeList)
            transcripts = sliceTemp

        random = args.get('random', False)
        if random:
            output = [
                self.addExtraWord(self.output[0], '_random'),
            ]
            self.saveOutput(output)

        averageLength = totalRecord = totalMappedReads = perNmappedReads= 1
        if self.runMode == True and self.runFlag == True:
            inputBed = bed.bed(self.input[0])
            totalMappedReads = inputBed.getHitNum()
            bedA = bed.bed(transcripts)
            averageLength = bedA.getAverageLength()
            totalRecord = bedA.getHitNum()
            perNmappedReads = 1000000
        
        shuffleCode = [
            'bedtools',
            'shuffle',
            '-i', transcripts,
            '-g', self.reference['limits'],
            '|',
            'sort',
            '-k1,1',
            '-k2,2n',
            '-k3,3n'
        ]
        if random:
            codeList = list(shuffleCode)
        else:
            codeList = [
                'cat',
                transcripts
            ]

        codeList += [
            '|',
            'bedtools',
            'intersect',
            '-a', 'stdin',
            '-b', self.input,
            '-wa',
            '-wb',
            '-F', 0.50,
            '|',
            'bedIntersect2positionTxt.py',
            '|',
            'bedIntersectPositionCount.py',
            '-count', 7,
            '-cat', 8,
            '-ends', 'double',
            '-scale',
            1/float(totalRecord),
            100/float(averageLength),
            perNmappedReads/float(totalMappedReads),
            '-o', self.output
        ]
        self.execM(codeList)
        
        _temp, temp = tempfile.mkstemp()

        oneSideFlankingLength = 2000
        windowLength = oneSideFlankingLength/100

        if random:
            prepareAbedCodeList = list(shuffleCode)
        else:
            prepareAbedCodeList = [
                'cat',
                transcripts
            ]
        prepareAbedCodeList += [
            '|',
            'bed2updownstream.py',
            '--fixed',
            '-l', oneSideFlankingLength,
            '|',
            'bed2removeChromosomeEdges.py',
            '--fixed',
            '-l', 10000,
            '-g', self.reference['limits'],
            '>', temp
        ]
        self.execM(prepareAbedCodeList)

        totalRecord = 1
        if self.runMode == True and self.runFlag == True:
            bedA = bed.bed(temp)
            totalRecord = bedA.getHitNum()/2

        flankingCodeList = [
            'bedtools',
            'intersect',
            # '-a', 'stdin',
            '-a', temp,
            '-b', self.input,
            '-wa',
            '-wb',
            '-F', 0.50,
            '|',
            'bedIntersect2positionTxt.py',
            '--flanking',
            '--fixed',
            '-w', windowLength,
            '|',
            'bedIntersectPositionCount.py',
            '-count', 7,
            '-cat', 8,
            '-ends', 'remove',
            '-scale', 
            1/float(totalRecord), 
            1/float(windowLength),
            perNmappedReads/float(totalMappedReads),
            '>>', self.output,
            '&&',
            'rm', temp
        ]
        self.execM(flankingCodeList)
        return self


    def mergeTCR(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_Q4", "*")
        headers = ['pos', 'count', 'cat'] + self.attributes +  ['label']
        output = os.path.join(self.outputDir, '..', 'merged_TCR' + extraWord + '.txt')
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
    return sampleIO(SAMPLE_STAT_FILE, n, 'no', 'SRA_id')


args = getArgs()
inputIndex = args.get("n")
input = getInputFromIndex(inputIndex)
###########################################################
#  Pipeline
###########################################################
p = myPipe(input, args)
(p
    # .branch(True)
    #     .run(p.fastqc_fastq2html, True)
    # .stop()
    .run(p.dump_2fastq, False)
    .run(p.cutadapt_fastq2fastq, False)
    .run(p.bowtie_fastq2sam, False)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bed, False)
    .run(p.sort_bed2bed, False)
    .run(p.slopBed_bed2bed, False)
    .run(p.convertToFixedRange_bed2bed, False)
    .run(p.sort_bed2bed, False)
    .run(p.convertBedToFasta_bed2fa, False)

    .branch(True) # Plot nucleotide abundance
        .run(p.getNucleotideAbundanceTable_fa2csv, False)
        .run(p.addTreatment_csv2txt, True)
        .cat(p.mergeNucleotideAbundance, True)
    .stop()

    .branch(True) # Plot dinucleotide abundance
        .run(p.getDimerAbundanceTable_fa2csv, False)
        .run(p.addTreatment_csv2txt, True)
        .cat(p.mergeNucleotideAbundance, True, '_diNuc')
    .stop()

    .run(p.getDamageSites_fa2bed, False)
    .run(p.sort_bed2bed, False)


    .branch(False)
        .run(p.writeTotalMappedReads_bed2txt, False)
    .stop()

    .branch(True)
        .run(p.geneMap_bed2txt, False)
        .run(p.normalizeCounts_txt2txt, False)
        .run(p.addTSNTS_txt2txt, True)
        .branch(False)
            .run(p.makeScoreBed6_txt2bed, True)
        .stop()
        .cat(p.mergeGeneCounts, False)
    .stop()

   # TCR Analysis
    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False)
        .run(p.addTreatment_txt2txt, True, 'real')
    .stop()

    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False, {'random':True})
        .run(p.addTreatment_txt2txt, True, 'random')
    .stop()

    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False, {"slice":True, "n":4, "sliceNum":1, "keyword":'_Q1'})
        .run(p.addTreatment_txt2txt, True, 'Q1')
    .stop()

    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False, {"slice":True, "n":4, "sliceNum":2, "keyword":'_Q2'})
        .run(p.addTreatment_txt2txt, True, 'Q2')
    .stop()

    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False, {"slice":True, "n":4, "sliceNum":3, "keyword":'_Q3'})
        .run(p.addTreatment_txt2txt, True, 'Q3')
    .stop()

    .branch(True)
        .run(p.transcriptIntersect_bed2txt, False, {"slice":True, "n":4, "sliceNum":4, "keyword":'_Q4'})
        .run(p.addTreatment_txt2txt, True, 'Q4')
        .cat(p.mergeTCR, True)
    .stop()

    .run(p.toBg_bed2bg, False)
    .run(p.toBw_bg2bw, False)

    .stop()
)
