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
import numpy

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        OUTPUT_DIR = '0825'
        SAMPLE_STAT_FILE = 'samples.csv'    
        self.input = os.path.join(os.path.curdir, '..', 'hiSeq', 'dataDir', 'raw', self.input)
        # self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        self.outputDir = os.path.realpath(os.path.join(os.path.curdir, 'dataDir', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.attributes = sorted(sampleDictionary.keys())
        self.paths = referenceGenomePath()
        
        input1 = self.input
        input2 = generalUtils.in2out(self.input, '.1.fastq', '.2.fastq')
        self.saveInput([input1, input2])
        
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.treatment + '.txt',
            '-e ': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams

        if self.product == 'CPD':
            self.motifRegex = '\'.{4}(TT|tt|CT|ct|cT|Ct).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
            # self.motifRegex = '\'.{4}(TT|tt|Tt|tT).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        elif self.product == '_6-4':
            self.motifRegex = '\'.{4}(TT|tt|Tt|tT|TC|tc|Tc|tC).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        elif self.product == 'GG':
            self.motifRegex = '\'.{4}(G|g)(G|g).{4}\'' # Get GG only at the positions 5 and 6.
        elif self.product == 'NA':
            self.motifRegex = '\'.{10}\'' # Get GG only at the positions 5 and 6.
        else:
            raise ValueError("Unknown product type. Exiting...")


    def prettyOutput(self, customWord = ''):
        newOutputs = []
        for o in self.output:
            if 'Plus' in o:
                extraWord = '_Plus'
            elif 'Minus' in o:
                extraWord = '_Minus'
            else:
                extraWord = ''

            if '_TS' in o:
                extraWord += '_TS'
            elif '_NTS' in o:
                extraWord += '_NTS'

            extension = pipeTools.getExtension(o)
            newOutputs.append(os.path.join(os.path.dirname(o),self.treatment_title + extraWord + customWord + '.' + extension))
        return newOutputs

    def tailFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "tail --lines=+2 " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        # code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        print(code)
        os.system(code)

    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        # code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        print(code)
        os.system(code)

    def fullPath2wildcard(self, fullPath):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join(["*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def mergeNucleotideFrequencies(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geSiLe" + str(self.fragmentLengths[0]), "*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength']
        output = os.path.join(self.outputDir, '..', 'merged_NucleotideFrequencies.txt')
        self.tailFiles(wildcard, headers, output)
        return self

    def cutadapt_fastq2fastq(self):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)
        adapter = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
        codeList = [
            'cutadapt',
            '--discard-trimmed', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
            '-g', adapter, # The adapter is located at the 5' end
            '-o', self.output[0],
            '-p', self.output[1],
            self.input[0],
            self.input[1]
        ]
        self.execM(codeList)
        return self
    
    def bowtie_fastq2sam(self, referenceGenome = 'hg19'):
        noCpus = 8
        self.reference = self.paths.get(referenceGenome)
        output = [self.output[0].replace(".sam", self.reference['name'] + ".sam")]
        self.saveOutput(output)
        codeList = [
            'bowtie',
            '-t', self.reference['bowtie'],
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            '-m', 4, # Do not report the reads that are mapped on to more than 4 genomic locations
            '-X', 1000,
            '--seed', 123, # Randomization parameter in bowtie,
            '-p', noCpus,
            '-1', self.input[0],
            '-2', self.input[1],
            self.output,
            '&>',
            self.output[0] + '.log'
        ]
        self.execM(codeList)
        return self

    def convertToBam_sam2bam(self):
        codeList = [
            'samtools',
            'view',
            '-bf', '0x2', #	each segment properly aligned according to the aligner
            #'-Sb'
            '-o',
            self.output,    
            self.input
        ]
        self.execM(codeList)
        return self
 
    def convertToBed_bam2bedpe(self):
        codeList = [
            'bedtools',
            'bamtobed',
            '-bedpe',
            '-mate1',
            '-i', self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self
    
    def uniqueSort_bedpe2bedpe(self):
        codeList = [
            'sort',
            '-u',
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            # '-k4,4',
            # '-k5,5',
            # '-k6,6',
            #'-k9,9',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertBedpeToSingleFrame_bedpe2bed(self):
        codeList = [
            'bedpe2bed.py',
            self.input,
            '>', self.output
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
    
    def sortBed_bed2bed(self):
        codeList = [
            'sort',
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        self.finalBed = self.output[0]
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
    
    def getPyrimidineDimers_fa2bed(self):
        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', self.motifRegex 
        ]
        self.execM(codeList)
        return self

    def makeBed6_bed2bed(self):
        codeList = [
            'bed4tobed6.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self
    
    def transcriptStrandCount_bed2txt(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        self.saveOutput(self.prettyOutput('_transcriptCounts'))
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["transcripts"],
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
 
    def relativeCount_txt2txt(self):
        codeList = [
            'bedCount2relativeRatio.py',
            '-i', self.input,
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

    def LADcount_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["LAD"],
            '-b', self.input,
            '-wa',
            '-c',
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self
    
    def binnedCountsToPositions_txt2txt(self):
        # codeList = [
        #     'bedBinned2totalCounts.py',
        #     '-i', self.input,
        #     '-o', self.output,
        #     '-n', 200,
        #     '-reverseStrand', '"-"'
        # ]
        startPosition = "-5000"
        codeList = [
            'bedBinned2totalCounts.py',
            '-i', self.input,
            '-o', self.output,
            '-n', 100,
            '-reverseStrand', '"-"',
            '--mergeStrands',
            '-start', startPosition,
            '-winLength', 100,
            '--writePosition'
        ]
        self.execM(codeList)
        return self

def getArgs():
    parser = argparse.ArgumentParser(description='UV Damage-seq Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-n', required= True, help='input index')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-n', required= False, default="1", help='input index')
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
    .run(p.cutadapt_fastq2fastq, False)
    .run(p.bowtie_fastq2sam, False)
    .run(p.convertToBam_sam2bam, False)
    .run(p.convertToBed_bam2bedpe, False)
    .run(p.uniqueSort_bedpe2bedpe, False)
    .run(p.convertBedpeToSingleFrame_bedpe2bed, False)
    .run(p.slopBed_bed2bed, False)
    .run(p.convertToFixedRange_bed2bed, False)
    .run(p.sortBed_bed2bed, False)
    .run(p.convertBedToFasta_bed2fa, False)
    .run(p.getPyrimidineDimers_fa2bed, False)
    .run(p.sortBed_bed2bed, False)

    .branch(False)
        .run(p.makeBed6_bed2bed, True)
        .run(p.transcriptStrandCount_bed2txt, True)
        .run(p.normalizeCounts_txt2txt, True)
        .run(p.relativeCount_txt2txt, True)
    .stop()

    .branch(True and p.cell == "GM12878") # LAD
        .run(p.LADcount_bed2txt, True)
        .run(p.normalizeCounts_txt2txt, True)
        .run(p.binnedCountsToPositions_txt2txt, True)
        .run(p.addTreatment_txt2txt, True)
    .stop()


)
