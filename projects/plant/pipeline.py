import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob
import argparse
import argument
sys.path.append('..')
from referenceGenomePath import referenceGenomePath

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        OUTPUT_DIR = '0707'
        SAMPLE_STAT_FILE = 'samples.csv'    
        self.input = os.path.join(os.path.curdir, 'dataDir', 'raw', self.input)
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.attributes = sorted(sampleDictionary.keys())
        self.saveInput([self.input])
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.treatment + '.txt',
            '-e ': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams
        paths = referenceGenomePath()
        self.reference = paths.get('TAIR9')
        self.fragmentLengths = range(17, 32+1)
    
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

    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "tail --lines=+2 " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
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
        wildcard = wildcard.replace("geSiLe17", "*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength']
        output = os.path.join(self.outputDir, '..', 'merged_NucleotideFrequencies.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def cutadapt_fastq2fastq(self):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)
        adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'
        codeList = [
            'cutadapt',
            '-a', adapter,
            '-o', output,
            self.input
        ]
        self.execM(codeList)
        return self

    def bowtie_fastq2sam(self):
        noCpus = 1
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0]]
        codeList = [
            'bowtie',
            '-t', self.reference['bowtie'],
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

    def splitByStrand_bed2bed(self):
        strand = ['+', '-']
        output = [
            self.addExtraWord(self.output[0], '_Plus'), 
            self.addExtraWord(self.output[0], '_Minus')
        ]
        self.saveOutput(output)
        self.saveOutput(self.prettyOutput())
        codeList = [
            'grep',
            strand,
            self.input[0],
            '>',
            self.output
        ]
        self.execM(codeList)
        return self

    def convertToBedGraph_bed2bdg(self):
        if self.runFlag and self.runMode:
            scaleFactor = float(1000000)/self.internalRun(bed.bed(self.input[0]).getHitNum, [], self.runFlag, 'get hit number')
        else:
            scaleFactor = 1
        codeList = [
            'bedtools',
            'genomecov',
            '-i', self.input,
            '-g', self.reference['limits'],
            '-bg',
            '-scale', scaleFactor,
            '>', self.output
        ]
        self.execM(codeList)
        return self
   
    def toBigWig_bdg2bw(self):
        codeList = [
            'bedGraphToBigWig',
            self.input,
            self.reference['limits'],
            self.output
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

    def getNucleotideAbundanceTable_fa2csv(self):
        # nucleotideOrder = 'GCTA'
        headerList = []
        columns = self.list2attributes(self.attributes)
        for i in range(len(columns)):
            headerList.append(sorted(self.attributes)[i] + ',' + columns[i])
        codeList = [
            'fa2kmerAbundanceMeltedData.py',
            '-i', self.input,
            '-o', self.output,
            '-k', 1,
            '--percentage',
            '-c', ' '.join(headerList)
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

    def getSingleLength_bed2bed(self):
        self.addWordsToOutput(self.fragmentLengths)
        print(self.output)
        codeList = [
            'bed2getCertainIntervalLengths.py',
            '-i', self.input * len(self.fragmentLengths),
            '-o', self.output,
            '-l', self.fragmentLengths
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


    def addTSNTS_txt2txt(self, site):
        columns = self.list2attributes(self.attributes)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', [ ' '.join(columns) +  ' TS ' + site, 
                    ' '.join(columns) +  ' NTS ' + site]
        ]
        self.execM(codeList)
        columnHeaders = self.attributes
        return self
   
    def addLength_csv2csv(self, *nargs):
        columns = self.list2attributes(self.attributes) + list(nargs)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', self.fragmentLengths
        ]
        self.execM(codeList)
        columnHeaders = self.attributes
        return self

    def geneStrandMap_bed2bed(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
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

    def TSScount_bed2txt(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["TSS"],
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def TEScount_bed2txt(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["TES"],
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

    def mergeGeneCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S")
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes + ['TSNTS']
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts.txt')
        self.catFiles(wildcard, headers, output)
        return self
        
    def mergeTSSTES(self):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S").replace("TSS", "*").replace("TES", "*")
        headers = ['position', 'value'] + self.attributes + ['strand', 'site']
        output = os.path.join(self.outputDir, '..', 'merged_TSSTES.txt')
        self.catFiles(wildcard, headers, output)
        return self
   
    def binnedCountsToPositions_txt2txt(self, startPosition):
        codeList = [
            'bedBinned2totalCounts.py',
            '-i', self.input,
            '-o', self.output,
            '-n', 60,
            '-reverseStrand', '"-"',
            '--mergeStrands',
            '-start', startPosition,
            '-winLength', 100,
            '--writePosition'
        ]
        self.execM(codeList)
        return self

    def chromosomeCounts_bed2txt(self):
        codeList = [
            'bed2chromosomeRPKM.py',
            '-i', self.input,
            '-g', self.reference['limits'],
            '-o', self.output
        ]
        self.execM(codeList)
        return self

def getArgs():
    parser = argparse.ArgumentParser(description='XR-seq Mouse Organs Pipeline', prog="pipeline.py")
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
    .run(p.convertToBed_bam2bed, False)
    .run(p.uniqueSort_bed2bed, False)
    
    .branch()
        .run(p.lengthDistribution_bed2csv, False)
        .run(p.plotLengthDistribution_csv2pdf, False)
    .stop()

    .branch(False)
        .run(p.getSingleLength_bed2bed, False)
        .run(p.convertBedToFasta_bed2fa, False)
        .run(p.getNucleotideAbundanceTable_fa2csv, False)
        .run(p.addLength_csv2csv, False)
        .cat(p.mergeNucleotideFrequencies, False)
    .stop()
   
    .branch(False)
        .run(p.splitByStrand_bed2bed, False)
        
        .branch(False)
            .run(p.convertToBedGraph_bed2bdg, True)
            .run(p.toBigWig_bdg2bw, True)
        .stop()
    .stop()

    .branch(False)
        .run(p.geneStrandMap_bed2bed, True)
        .run(p.normalizeCounts_txt2txt, True)
        
        .branch()
            .run(p.addTreatmentAndStrand_txt2txt, True)
            .cat(p.mergeGeneCounts, False)
        .stop()
    .stop()    

    .branch(False)
        .run(p.TSScount_bed2txt, True)
        .run(p.normalizeCounts_txt2txt, True)
        
        .branch()
            .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"')
            .run(p.addTSNTS_txt2txt, True, "TSS")
            # .cat(p.mergeTSSTES, True)
        .stop()
    .stop()

    .branch(False)
        .run(p.TEScount_bed2txt, True)
        .run(p.normalizeCounts_txt2txt, True)
        
        .branch()
            .run(p.binnedCountsToPositions_txt2txt, True, -4000)
            .run(p.addTSNTS_txt2txt, True, "TES")
            .cat(p.mergeTSSTES, True)
        .stop()
    .stop() 

    .branch(True)
        .run(p.chromosomeCounts_bed2txt, True)
    .stop()

    

)
