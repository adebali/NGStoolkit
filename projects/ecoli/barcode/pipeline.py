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
# subsample

class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        SAMPLE_STAT_FILE = 'samples.csv'    
        OUTPUT_DIR = "0405"
        ADAPTER_SEQUENCE = "GGCTCAGTTCGTATGAGTGCCG"
        pipe.__init__(self, input, args)
        self.input = os.path.realpath(os.path.join(os.path.curdir, 'dataDir', 'raw', input))
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        input = generalUtils.file(os.path.realpath(os.path.join('dataDir', 'raw', input)))
        self.saveInput([self.input])

        self.ecoliReferenceRoot = '/nas02/home/a/d/adebali/ncbi/ecoli'
        self.referenceVersion = 'NC_000913.2'
        self.referenceRoot = os.path.join(self.ecoliReferenceRoot, self.referenceVersion)
        self.referenceBowtieIndex = os.path.join(self.referenceRoot, self.referenceVersion)
        self.referenceFasta = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.fa'))
        self.referenceGeneTTcontent = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.genes.AA_TT.bed'))
        self.referenceGeneNames = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.genes.txt'))
        self.referenceGenesRNA = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '_RNAseq_SRR1173967_PlusMinus.txt'))
        self.referenceChromosomeLimits = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.sizes'))
        self.referenceGenesBed = generalUtils.file(os.path.join(self.referenceRoot, self.referenceVersion + '.chr.genes.6col.bed'))

        self.TSSfile = os.path.join('genomeFiles', 'TSS_1h1k_noOverlap.bed')

        self.temp = os.path.join('/pine/scr/a/d/adebali', self.treatment_title + '.txt')

        self.totalMappedReads = None
        self.adapter = ADAPTER_SEQUENCE

    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "tail --lines=+2 " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        os.system(code)


    def fullPath2wildcard(self, fullPath):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join(["*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def mergeUniqueCounts(self):
        wildcard = os.path.join(self.outputDir, "*.unReCo.adTr.csv")
        headers = ['value'] + self.attributes + ['status']
        output = os.path.join(self.outputDir, '..', 'merged_UniqueCounts.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeGeneCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['chr', 'start', 'end', 'x1', 'x2', 'strand', 'TSrep', 'NTSrep'] + self.attributes + ['name', 'pTT', 'mTT', 'pRNA', 'mRNA']
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeTSS(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("_TS", "*")
        headers = ['position', 'value', 'name'] + self.attributes + ['strand']
        output = os.path.join(self.outputDir, '..', 'merged_TSS.txt')
        self.catFiles(wildcard, headers, output)
        return self
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


    def writeSeqToHeader_fastq2fastq(self):
        self.output = [os.path.join(self.outputDir,os.path.basename(self.output[0]))]
        self.saveOutput(self.output)
        codeList = [
            'fastq2addSeqToHeader.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def cutAdaptor_fastq2fastq(self):
        codeList = [
            'cutadapt',
            '-a', self.adapter, # The adapter is located at the 3' end
            '-o', self.output,
            self.input
        ]
        self.execM(codeList)
        return self
    
    def getLengthDistribution_fastq2csv(self):
        codeList = [
            # 'head', self.input, '|'
            "echo -e 'length\\tcount'", '>', self.output,
            '&&',
            "awk 'BEGIN { OFS = \"\t\" } NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'",
            self.input,
            '>>', self.output
        ]
        self.coreColumns = ['length', 'count']
        self.mergedFileName = "lengthDistribution"
        self.execM(codeList)
        return self

    def addTreatment_csv2csv(self, *nargs):
        columns = self.list2attributes(self.attributes) + list(nargs)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', ' '.join(columns)
        ]
        self.execM(codeList)
        columnHeaders = self.attributes

        if self.runFlag and self.args.get("subprogram") == "cat":
            # f = self.output[0]
            # basename = os.path.basename(f)
            # directory = os.path.dirname(f)
            # wildcard = ".".join(["*"] + basename.split(".")[1:])
            # fullWildcard = os.path.join(directory, wildcard)
            fullWildcard = self.fullPath2wildcard(self.output[0])
            outputFileName = self.mergedFileName
            outputFile = os.path.join(self.outputDir, "..", 'merged_' + outputFileName + '.txt')
            headers = self.coreColumns + columnHeaders
            # header = '\t'.join(self.coreColumns + columnHeaders)
            # code = "echo -e '" + header + "' >" + outputFile + " & " + "tail --lines=+2 " + fullWildcard + " | grep -v '^==' | grep -v -e '^$' >>" + outputFile
            self.catFiles(fullWildcard, headers, outputFile)
            # print(code)
            # os.system(code)
        return self

    def addGeneData_csv2csv(self):
        codeList = [
            'pasteLastColumn.py',
            '-i', self.input,
            self.referenceGeneNames,
            self.referenceGeneTTcontent,
            self.referenceGenesRNA,
            '-n', 2,
            '>', self.output,
        ]
        self.execM(codeList)
        return self

    def singleLine_fastq2csv(self):
        codeList = [
            'fastq2singleLine.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def fourLines_csv2fastq(self):
        codeList = [
            'fastq2singleLine.py',
            '-i', self.input,
            '-o', self.output,
            '--reverse'
        ]
        self.execM(codeList)
        return self

    def sortBySequence_csv2csv(self):
        codeList = [
            'sort',
            "-t$'\\t'",
            '-k2',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getAnalytics_csv2csv(self):
        codeList = [
            './barcodeAnalysis.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def reduceArtifacts_csv2csv(self):
        codeList = [
            './reduceArtifacts.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def uniqReadCount_csv2csv(self):
        codeList = [
            'cut',
            '-f', 2,
            self.input,
            "| uniq -c | tr -s ' ' | cut -d ' ' -f 2",
            '>', self.output
        ]
        self.execM(codeList)
        self.mergedFileName = "uniqueReadCount"
        return self

    def align_fastq2sam(self):
        codeList = [
            'bowtie',
            '-t', self.referenceBowtieIndex,
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            #'--all', 
            #'--strata',
            #'--best',
            '--seed 123', # Randomization parameter in bowtie,
            self.input,
            self.output
        ]
        self.execM(codeList)
        return self

    def toBam_sam2bam(self):
        codeList = [
            'samtools',
            'view',
            '-bS',
            '-o',
            self.output,
            self.input
        ]
        self.execM(codeList)
        return self

    def sortBam_bam2bam(self):
        codeList = [
            'samtools',
            'sort',
            '-T', self.temp, # temporary file
            '-o', self.output,
            self.input
        ]
        self.execM(codeList)
        return self

    def getIndex_xbam2bai(self):
        codeList = [
            'samtools',
            'index',
            self.input,
            self.output
        ]
        self.execM(codeList)
        return self

    def toBed_bam2bed(self):
        codeList = [
            'bedtools',
            'bamtobed',
            '-i', self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getFasta_bed2fa(self):
        codeList = [
            'bedtools',
            'getfasta',
            '-fi', self.referenceFasta,
            '-bed', self.input,
            '-fo', self.output,
            '-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
        ]
        self.execM(codeList)
        return self

    def get13mer_fastq2fastq(self):
        codeList = [
            'fastq2getCertainLengths.py',
            '-i', self.input,
            '-o', self.output,
            '-l', 13
        ]
        self.execM(codeList)
        return self

    def nucleotideFrequency_fa2csv(self):
        codeList = [
            'fa2kmerAbundanceMeltedData.py',
            '-i', self.input,
            '-o', self.output,
            '-k', 1,
            '--percentage'
        ]
        self.coreColumns = ['position', 'sequence', 'value']
        self.mergedFileName = "singleNucleotide"
        self.execM(codeList)
        return self

    def diNucleotideFrequency_fa2csv(self):
        codeList = [
            'fa2kmerAbundanceMeltedData.py',
            '-i', self.input,
            '-o', self.output,
            '-k', 2,
            '--percentage'
        ]
        self.coreColumns = ['position', 'sequence', 'value']
        self.mergedFileName = "diNucleotide"
        self.execM(codeList)
        return self

    def getTT_fa2bed(self):
        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', '\'.{7}(TT|Tt|tT).{4}\''
        ]
        self.execM(codeList)
        return self

    def sortAnalytics_csv2csv(self):
        codeList = [
            'sort',
            "-t$'\\t'",
            '-k9,9n',
            self.input,
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
    
    def toBam_bed2bam(self):
        codeList = [
            'bedToBam',
            '-i', self.input,
            '-g', self.referenceChromosomeLimits,
            '>', self.output
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
            '-g', self.referenceChromosomeLimits,
            '-bg',
            '-scale', scaleFactor,
            '>', self.output
        ]
        self.execM(codeList)
        return self


    def convertToWig_bed2wig(self):
        codeList = [
            'igvtools',
            'count',
            self.input,
            self.output,
            self.referenceChromosomeLimits
        ]
        self.execM(codeList)
        return self

    def toBigWig_bdg2bw(self):
        codeList = [
            'bedGraphToBigWig',
            self.input,
            self.referenceChromosomeLimits,
            self.output
        ]
        self.execM(codeList)
        return self

    def geneCount_bed2bdg(self):
        self.input.append(self.input[0])
        self.saveOutput([self.addExtraWord(self.output[0],'_TS'), self.addExtraWord(self.output[0],'_NTS')])
        strandedness = ['-S', '-s']
        codeList = [
            'bedtools',
            'coverage',
            '-counts',
            strandedness,
            '-a', self.referenceGenesBed,
            '-b', self.input,
            '>', self.output
        ]
        self.execM(codeList)

    def normalizeGeneCounts_bdg2bdg(self):
        try:
            totalMappedReads = int(open(self.totalMappedReadsFile, 'r').readline())
        except:
            totalMappedReads = 1000000
        codeList = [
            'bedCount2normalizedCount.py',
            '-i', self.input,
            '-o', self.output,
            '-c', 7,
            '-l', 1000,
            '-m', 1000000/float(totalMappedReads)
        ]
        self.execM(codeList)
        return self

    def mergeStrands_bdg2csv(self):
        self.saveOutput([self.output[0].replace("_TS","")])
        codeList = [
            'pasteLastColumn.py',
            '-inputs',
            self.input[0], self.input[1],
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def slopDamage_bed2bed(self):
        codeList = [
            'bedtools', 
            'slop',
            '-i', self.input,
            '-g', self.referenceChromosomeLimits,
            '-s',
            '-l', -7,
            '-r', -4,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def intersectTSS_bed2txt(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.TSSfile,
            '-b', self.input,
            '-wa',
            '-wb',
            strandParameters,
            '-f', 0.001,
            '>', self.output
        ]
        self.execM(codeList)
        self.frame = [-100, 1000]
        return self

    def getPositions_txt2txt(self):
        codeList = [
            'bedIntersect2parsedPosition.py',
            '-i', self.input,
            '--reverse',
            '--name',
            '-l', self.frame[1] - self.frame[0] + 1,
            # '-rpm', self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number'),
            # '-rpm', self.internalRun(bed.bed(self.finalBed).getHitNum, [], True, 'get hit number', True),
            # '-rpm', self.internalRun(bed.bed(self.input[0]).getHitNum, [], self.runFlag, 'get hit number', True),
            '-rpm', 1000000,
            '-nameNorm', 1000,
            '>', self.output
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

def getArgs():
    parser = argparse.ArgumentParser(description='XR-seq Barcode Pipeline', prog="pipeline.py")
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
    .run(p.writeSeqToHeader_fastq2fastq, True)

    .branch()
        .run(p.cutAdaptor_fastq2fastq, True)

        .branch()
            .run(p.getLengthDistribution_fastq2csv, True)
            .run(p.addTreatment_csv2csv, True)
        .stop()

        .run(p.get13mer_fastq2fastq, True)

        .branch()
            .run(p.align_fastq2sam, True)
            .run(p.toBam_sam2bam, True)
            .run(p.sortBam_bam2bam, True)

            .branch(False)
                .run(p.getIndex_xbam2bai, False)
            .stop()

            .run(p.toBed_bam2bed, True)
            .run(p.getFasta_bed2fa, True)

            .branch()
                .run(p.nucleotideFrequency_fa2csv, True)
                .run(p.addTreatment_csv2csv, True)
            .stop()

            .branch()
                .run(p.diNucleotideFrequency_fa2csv, True)
                .run(p.addTreatment_csv2csv, True)
            .stop()

            .run(p.getTT_fa2bed, True)
            
            .branch(True)
                .run(p.writeTotalMappedReads_bed2txt, True)
            .stop()
            
            .branch(True)
                .run(p.geneCount_bed2bdg, True)
                .run(p.normalizeGeneCounts_bdg2bdg, True)
                .run(p.mergeStrands_bdg2csv, True)
                .run(p.addTreatment_csv2csv, True)
                .run(p.addGeneData_csv2csv, True)
                .cat(p.mergeGeneCounts, True)
            .stop()

        .stop()


        .branch(False)
            .run(p.singleLine_fastq2csv, False)
            .run(p.sortBySequence_csv2csv, False)

            .branch(True)
                .run(p.uniqReadCount_csv2csv, False)
                .run(p.addTreatment_csv2csv, False, 'initial')
            .stop()
            
            .branch(True)
                .run(p.reduceArtifacts_csv2csv, False)
                
                .branch(False)
                    .run(p.uniqReadCount_csv2csv, False)
                    .run(p.addTreatment_csv2csv, False, 'post')
                    .cat(p.mergeUniqueCounts, False)
                .stop()

                .run(p.fourLines_csv2fastq, False)
                .run(p.align_fastq2sam, False)
                .run(p.toBam_sam2bam, False)
                .run(p.sortBam_bam2bam, False)
                .run(p.toBed_bam2bed, False)
                .run(p.getFasta_bed2fa, False)
                .run(p.getTT_fa2bed, False)

                .branch(False)
                    .run(p.writeTotalMappedReads_bed2txt, False)
                .stop()

                .branch(True)
                    .run(p.splitByStrand_bed2bed, False)

                    .branch(True)
                        .run(p.toBam_bed2bam, False)
                   
                        .branch(True)
                            .run(p.getIndex_xbam2bai, False)
                        .stop()
                    .stop()

                    .branch(False)
                        .run(p.convertToBedGraph_bed2bdg, True)
                        .run(p.toBigWig_bdg2bw, True)
                    .stop()
                .stop()

                .branch(False)
                    .run(p.slopDamage_bed2bed, False)
                    .run(p.intersectTSS_bed2txt, False)
                    .run(p.getPositions_txt2txt, False)
                    .run(p.addTreatmentAndStrand_txt2txt, False)
                    .cat(p.mergeTSS, False)
                .stop()

                .branch(True)
                    .run(p.geneCount_bed2bdg, False)
                    .run(p.normalizeGeneCounts_bdg2bdg, True)
                    .run(p.mergeStrands_bdg2csv, True)
                    .run(p.addTreatment_csv2csv, True)
                    .run(p.addGeneData_csv2csv, True)
                    .cat(p.mergeGeneCounts, False)
                .stop()

            .stop()

            .run(p.getAnalytics_csv2csv, False)
            .run(p.sortAnalytics_csv2csv, False)
        .stop()

        
    .stop()
)
# def f(x):
#     return (2e-17) * x**6 - (1e-13) * x**5 + (5e-10)*x**4 - (8e-07)*x**3 + (0.0008)*x**2 + (0.8324)*x + 20.877
# # y = 2E-17x6 - 1E-13x5 + 5E-10x4 - 8E-07x3 + 0.0008x2 + 0.8324x + 20.877
# y = 1.8501x - 542.29
# 240.49e0.0013x
