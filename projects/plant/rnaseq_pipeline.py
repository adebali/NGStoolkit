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
import numpy
import tempfile


SAMPLE_STAT_FILE = 'rnaSamples.csv'

class pipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        OUTPUT_DIR = 'RNA-seq'
        # SAMPLE_STAT_FILE = 'rnaSamples.csv'
        self.input = os.path.join(os.path.curdir, 'dataDir', 'RNA-seq', self.input)
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
            '--job-name=': 'RNA-seq',
            '-o ': 'log_' + self.treatment + '.txt',
            '-e ': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams
        self.paths = referenceGenomePath()
        
        self.fragmentLengths = range(10, 32+1)
        # self.fragmentLengths = range(10, 17+1)
    
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

    def fullPath2wildcard(self, fullPath, prefix = ''):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join([prefix + "*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def mergeNucleotideFrequencies(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geSiLe" + str(self.fragmentLengths[0]), "*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength']
        output = os.path.join(self.outputDir, '..', 'merged_NucleotideFrequencies' + extraWord + '.txt')
        self.tailFiles(wildcard, headers, output)
        return self

    def mergeChrNucleotideFrequencies(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geSiLe" + str(self.fragmentLengths[0]), "*").replace("seCh_1", "seCh_*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength', 'chromosome']
        output = os.path.join(self.outputDir, '..', 'merged_ChrNucleotideFrequencies.txt')
        self.tailFiles(wildcard, headers, output)
        return self

    def mergeNucleotide(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes)
        output = os.path.join(self.outputDir, '..', 'merged_nucleotide.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeLengthDist(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geSiLe" + str(self.fragmentLengths[0]), "*").replace("seCh_1", "seCh_*")
        headers = ['length', 'value'] + sorted(self.attributes) + ['fragmentLength', 'chromosome']
        output = os.path.join(self.outputDir, '..', 'merged_lengthDistributions.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeLength(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['length', 'value'] + sorted(self.attributes)
        output = os.path.join(self.outputDir, '..', 'merged_length.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeChrLengthDist(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geSiLe" + str(self.fragmentLengths[0]), "*").replace("seCh_1", "seCh_*")
        headers = ['length', 'value'] + sorted(self.attributes) + ['chromosome']
        output = os.path.join(self.outputDir, '..', 'merged_chrLengthDistributions.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeChromtainStates(self, outputName):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['chromosome', 'start', 'end', 'state', 'score', 'strand', 'count'] + self.attributes
        output = os.path.join(self.outputDir, '..', outputName)
        self.catFiles(wildcard, headers, output)
        return self

    def mergeGeneCounts(self, outputName):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "*")
        headers = ['chromosome', 'start', 'end', 'gene', 'score', 'strand', 'count'] + self.attributes  + ['strand']
        output = os.path.join(self.outputDir, '..', outputName)
        self.catFiles(wildcard, headers, output)
        return self

    def mergeDNaseCounts(self, outputName):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['position', 'count'] + self.attributes 
        output = os.path.join(self.outputDir, '..', outputName)
        self.catFiles(wildcard, headers, output)
        return self

    def cutadapt_fastq2fastq(self):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)
        adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG'
        # adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'
        codeList = [
            'cutadapt',
            '-a', adapter,
            '-o', output,
            self.input
        ]
        self.execM(codeList)
        return self

    def bowtie_fa2sam(self, referenceGenome = "TAIR10"):
        noCpus = 8
        self.reference = self.paths.get(referenceGenome)
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0].replace(".sam", self.reference['name'] + ".sam")]
        self.saveOutput(output)
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)        
        codeList = [
            'bowtie',
            '-t', self.reference['bowtie'],
            '-f', # FASTA input
            '-p',8,
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

    def sort_bed2bed(self):
        codeList = [
            'sort',
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

    def get27mer_bed2bed(self):
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

    def get27mer_bed2bed(self):
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

    def addValue_txt2txt(self, *nargs):
        columns = list(nargs)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', ' '.join(columns)
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
   
    def addLength_csv2csv(self, *nargs):
        columns = self.list2attributes(self.attributes) + list(nargs)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', self.fragmentLengths * (len(self.input) / len(self.fragmentLengths))
        ]
        self.execM(codeList)
        columnHeaders = self.attributes
        return self

    def addLengthChr_csv2csv(self, *nargs):
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', [str(a) + ' ' + str(b) for a,b in zip(list(numpy.repeat(self.fragmentLengths, len(self.chromosomes))),self.chromosomes * (len(self.input) / len(self.chromosomes)))] 
        ]
        self.execM(codeList)
        return self

    def addChr_csv2csv(self, *nargs):
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', self.chromosomes
        ]
        self.execM(codeList)
        return self


    def geneMap_bed2txt(self, geneReference):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference[geneReference],
            '-b', self.input,
            '-wa',
            '-c',
            # strandParameters,
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self
  
    def dnaseCount_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["dnase"],
            '-b', self.input,
            '-wa',
            '-c',
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def transcriptSiteCount_bed2txt(self, siteKey): 
        """
        docstring here
            :param self: 
            :param siteKey: 
        """   
               
        TSSreference = self.reference[siteKey]
        newOutput = [self.addExtraWord(self.output[0], '_TS' + '_' + siteKey), self.addExtraWord(self.output[0], '_NTS' + '_' + siteKey)]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', TSSreference,
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def transcriptIntersect_bed2txt(self, args={}):
        transcripts = self.reference[args.get('transcripts', 'transcripts2')]
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

        windowLength = 100
        oneSideFlankingLength = 10000

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


    def dnaseIntersect_bed2txt(self, args={}):
        bedAfile = self.reference['dnase']

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
            bedA = bed.bed(bedAfile)
            averageLength = bedA.getAverageLength()
            totalRecord = bedA.getHitNum()
            perNmappedReads = 1000000
        
        shuffleCode = [
            'bedtools',
            'shuffle',
            '-i', bedAfile,
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
                bedAfile
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
            '-ends', 'remove',
            '-scale',
            1/float(totalRecord),
            100/float(averageLength),
            perNmappedReads/float(totalMappedReads),
            '-o', self.output
        ]
        self.execM(codeList)
        
        _temp, temp = tempfile.mkstemp()

        windowLength = 30
        oneSideFlankingLength = 3000

        if random:
            prepareAbedCodeList = list(shuffleCode)
        else:
            prepareAbedCodeList = [
                'cat',
                bedAfile
            ]
        prepareAbedCodeList += [
            '|',
            'bed2updownstream.py',
            '--fixed',
            '-l', oneSideFlankingLength,
            '|',
            'bed2removeChromosomeEdges.py',
            '--fixed',
            '-l', oneSideFlankingLength,
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

    def chmmIntersect_bed2txt(self, args={}):
        bedAfile = self.reference[args.get('bedA', 'chmm36')]

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
            bedA = bed.bed(bedAfile)
            averageLength = bedA.getAverageLength()
            totalRecord = bedA.getHitNum()
            perNmappedReads = 1000000
        
        shuffleCode = [
            'bedtools',
            'shuffle',
            '-i', bedAfile,
            '-g', self.reference['limits'],
            '|',
            'sort',
            '-k1,1',
            '-k2,2n',
            '-k3,3n'
        ]


        codeList = [
            'cat',
            bedAfile
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
            # '>', 'output.1.txt',
            # '&&',
            # 'cat', 'output.1.txt'
            '|',
            'bedIntersectPositionCount.py',
            '-count', 7,
            '-cat', 4, 8,
            '-ends', 'remove',
            '-scale',
            1/float(totalRecord),
            100/float(averageLength),
            perNmappedReads/float(totalMappedReads),
            '-o', self.output
        ]
        self.execM(codeList)
        
        _temp, temp = tempfile.mkstemp()

        windowLength = 10
        oneSideFlankingLength = 1000

        # if random:
            # prepareAbedCodeList = list(shuffleCode)
        # else:
        prepareAbedCodeList = [
            'cat',
            bedAfile
        ]
        prepareAbedCodeList += [
            '|',
            'bed2updownstream.py',
            '--fixed',
            '-l', oneSideFlankingLength,
            '|',
            'bed2removeChromosomeEdges.py',
            '--fixed',
            '-l', oneSideFlankingLength,
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
            '-cat', 4, 8,
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

    def epigeneticIntersect_bed2txt(self, args={}):
        return self.chmmIntersect_bed2txt(args)

    def OORintersect_bed2txt(self, args={}):
        return self.chmmIntersect_bed2txt(args)

    def TEScount_bed2txt(self, TESkey):
        TESreference = self.reference[TESkey]
        newOutput = [self.addExtraWord(self.output[0], '_TS' + '_' + TESkey), self.addExtraWord(self.output[0], '_NTS' + '_' + TESkey)]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', TESreference,
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

    def mergeTranscriptCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S")
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes + ['TSNTS']
        output = os.path.join(self.outputDir, '..', 'merged_transcriptCounts.txt')
        self.catFiles(wildcard, headers, output)
        return self
        
    def mergeTSSTES(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S").replace("TSS", "*").replace("TES", "*").replace("Q4","*")
        headers = ['position', 'value'] + self.attributes + ['strand', 'site', 'quartile']
        output = os.path.join(self.outputDir, '..', 'merged_TSSTES' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeTCR(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_no_neighbor", "*")
        headers = ['pos', 'count', 'cat', 'label']
        output = os.path.join(self.outputDir, '..', 'merged_TCR' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeDNase(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_random", "*")
        headers = ['pos', 'count', 'cat', 'label']
        output = os.path.join(self.outputDir, '..', 'merged_DNaseIntersect' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeChmm(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_random", "*")
        headers = ['pos', 'count', 'cat1', 'cat2', 'label']
        output = os.path.join(self.outputDir, '..', 'merged_CHMMintersect' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeEpigenetic(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_random", "*")
        headers = ['pos', 'count', 'cat1', 'cat2', 'label'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_epigenetic' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeOOR(self, extraWord = '', prefix = ''):
        wildcard =  self.fullPath2wildcard(self.input[0], prefix).replace("_random", "*")
        headers = ['pos', 'count', 'cat1', 'cat2', 'label'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_OOR' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

   
    def binnedCountsToPositions_txt2txt(self, startPosition, numberOfWindows = 80, windowsLength = 100):
        codeList = [
            'bedBinned2totalCounts.py',
            '-i', self.input,
            '-o', self.output,
            '-n', numberOfWindows,
            '-reverseStrand', '"-"',
            '--mergeStrands',
            '-start', startPosition,
            '-winLength', windowsLength,
            '--writePosition',
            '--average'
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

    def separateChromosomes_bed2bed(self):
        chromosomes = []
        outputs = []
        outputStringList = []
        inputs = []
        filein = open(self.reference['limits'], 'r')
        for line in filein:
            chromosome = line.split('\t')[0].strip()
            chromosomes.append(chromosome)
        self.chromosomes = chromosomes
        for i in range(len(self.input)):
            outString = ''
            for chromosome in chromosomes:
                inputs.append(self.input[i])
                outputName = self.output[i].replace('.bed', '_' + chromosome + '.bed')
                outString += outputName + ' '
                outputs.append(outputName)
            else:
                outputStringList.append(outString)

        self.saveOutput(outputs)
        print(len(inputs))
        print(len(outputs))
        codeList = [
            'bed2chromosomeSeparatedBed.py',
            '-i', self.input,
            '-chr', ' '.join(chromosomes),
            '-o', outputStringList
        ]
        # print(codeList)
        self.execM(codeList)
        return self

    def geneCount_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["genes"],
            '-b', self.input,
            '-c',
            '-F', 0.49,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countChromatinStates_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["chromatinStates"],
            '-b', self.input,
            '-c',
            '-F', 0.49,
            '>', self.output 
        ]
        self.execM(codeList)
        return self
   
    def subtractTCR_bed2bed(self):
        codeList = [
            'bedtools',
            'subtract',
            '-a', self.input,
            '-b', self.reference["transcripts"],
            '-S', # Force opposite strand
            '-A', # Remove if any overlap
            '>', self.output 
        ]
        self.execM(codeList)
        return self
   
    def addArtificialGlobalRepair_bed2bed(self):
        codeList = [
            'python addArtificialGlobalRepair.py',
            '-i', self.input,
            '-o', self.output,
            '-transcribed', self.reference["transcripts"]
        ]
        self.execM(codeList)
        return self
   
    def getIntersectPositions_txt2txt(self):
        codeList = [
            'bedIntersect2positionTxt.py',
            '-i', self.input,
            '-flanking', 'True',
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def countPositions_txt2txt(self):
        codeList =[
            'bedIntersectPositionCount.py',
            '-i', self.input,
            '-count', 7,
            '-cat', 8,
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
    # SAMPLE_STAT_FILE = 'samples.csv'
    return sampleIO(SAMPLE_STAT_FILE, n, 'no', 'sample')

args = getArgs()
inputIndex = args.get("n")
input = getInputFromIndex(inputIndex)


###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    p = pipeline(input, args)
    (p

        .run(p.bowtie_fa2sam, True, 'TAIR10')
        .run(p.convertToBam_sam2bam, True)
        .run(p.convertToBed_bam2bed, True)
        .run(p.uniqueSort_bed2bed, True)

    # Annotated transcript counts
        .branch(True)
            .run(p.geneMap_bed2txt, True, "transcripts")
            .run(p.normalizeCounts_txt2txt, True)
            .run(p.addTreatment_txt2txt, True)
            .cat(p.mergeGeneCounts, True, 'merged_RNAcounts.txt')                     
        .stop()
    )
###########################################################
