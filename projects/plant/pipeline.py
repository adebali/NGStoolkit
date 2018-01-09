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

class pipeline(pipe):
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
        self.paths = referenceGenomePath()
        
        self.fragmentLengths = range(10, 32+1)
        self.reference = self.paths.get('TAIR10')
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

    def mergeNucleotideFrequenciesAll(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geAlLe" + str(self.fragmentLengths[0]), "*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength']
        output = os.path.join(self.outputDir, '..', 'merged_NucleotideFrequencies' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
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
    
    def mergeBinCounts(self, outputName):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes
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

    def bowtie_fastq2sam(self, referenceGenome = "TAIR9"):
        noCpus = 8
        self.reference = self.paths.get(referenceGenome)
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0].replace(".sam", self.reference['name'] + ".sam")]
        self.saveOutput(output)
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

    def getAllLengths_bed2bed(self):
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


    def geneStrandMap_bed2txt(self, geneReference):
        newOutput = [self.addExtraWord(self.output[0], '_' + geneReference + '_TS'), self.addExtraWord(self.output[0], '_' + geneReference + '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference[geneReference],
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self
  

    def binCount_bed2txt(self, binSize = 1000):
        newOutput = [self.addExtraWord(self.output[0], '_' + str(binSize))]
        self.saveOutput(newOutput)

        codeList = [
            'bed2makeWindows.py',
            '-i', self.reference['limits_bed'],
            '-w', binSize,
            '|',
            'bedtools',
            'intersect',
            '-a', 'stdin',
            '-b', self.input,
            '-wa',
            '-c',
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

    def GEO_bed2txt(self):
        GEOtarget = self.reference['geo']
        def getPath(var):
            return os.path.join(GEOtarget, os.path.basename(var)).replace('.txt', '.md5.txt')

        def joinOutput(var):
            return os.path.join(GEOtarget, os.path.basename(var))

        inputList = pipeTools.listOperation(getPath, self.input)
        outputList = pipeTools.listOperation(joinOutput, self.output)
        outputList2 = pipeTools.listOperation(getPath, self.output)

        codeList = [
            'cp',
            self.input,
            outputList,
            '&&',
            'md5sum',
            outputList,
            '>',
            outputList2       
        ]
        self.saveOutput(outputList)
        self.execM(codeList)
        return self

    def GEO_fastq2txt(self):
        self.GEO_bed2txt()
        return self

    def GEO_txt2txt(self):
        self.GEO_bed2txt()
        return self

    def adjustColumns_txt2txt(self):
        def getPath(var):
            return var.replace('.txt', '.md5.txt')
        outputList2 = pipeTools.listOperation(getPath, self.output)
        
        codeList = [
            'echo -e "chromosome\\tstart\\tend\\tname\\tscore\\tstrand"',
            '>', self.output,
            '&&',
            'awk \'{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$7"\\t"$6}\'',
            self.input,
            '>>',
            self.output,
            '&&',
            'md5sum',
            self.output,
            '>',
            outputList2    
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
if __name__ == "__main__":
    p = pipeline(input, args)
    (p
        .branch(False)
            .run(p.GEO_fastq2txt, True)
        .stop()
        .run(p.cutadapt_fastq2fastq, False)
    # TAIR 9

        .branch(False)
            .run(p.bowtie_fastq2sam, False, 'TAIR9')
            .run(p.convertToBam_sam2bam, False)
            .run(p.convertToBed_bam2bed, False)
            .run(p.uniqueSort_bed2bed, False)
            
            .branch(False)
                .run(p.lengthDistribution_bed2csv, False)
                .run(p.plotLengthDistribution_csv2pdf, False)
            .stop()

            # .branch(False)
            #     .run(p.separateChromosomes_bed2bed, False)

            #     .branch(False)
            #         .run(p.lengthDistribution_bed2csv, False)
            #         .run(p.addTreatment_txt2txt, False)
            #         .run(p.addChr_csv2csv, False)
            #         .cat(p.mergeChrLengthDist, False)                
            #     .stop()
            # .stop()

            # .branch(True)
            #     .run(p.get27mer_bed2bed, False)

            #     .branch(False)
            #         .run(p.separateChromosomes_bed2bed, False)

            #         .branch(False)
            #             .run(p.convertBedToFasta_bed2fa, False)
            #             .run(p.getNucleotideAbundanceTable_fa2csv, False)
            #             .run(p.addLengthChr_csv2csv, False)
            #             .cat(p.mergeChrNucleotideFrequencies, False)
            #         .stop()
            #     .stop()
                
            #     .branch(False)
            #         .run(p.convertBedToFasta_bed2fa, False)
            #         .run(p.getNucleotideAbundanceTable_fa2csv, False)
            #         .run(p.addLength_csv2csv, False)
            #         .cat(p.mergeNucleotideFrequencies, False)
            #     .stop()
                            
            #     .branch(True)
            #         .run(p.convertBedToFasta_bed2fa, True)
            #         .run(p.getNucleotideAbundanceTable_fa2csv, True, 2)
            #         .run(p.addLength_csv2csv, True)
            #         .cat(p.mergeNucleotideFrequencies, True, '_diNuc')
            #     .stop()
            # .stop()
        
            .branch(False)
                .run(p.splitByStrand_bed2bed, True)
                
                .branch(False)
                    .run(p.convertToBedGraph_bed2bdg, True)
                    .run(p.toBigWig_bdg2bw, True)
                .stop()
            .stop()

            # .branch(True)
            #     .run(p.geneStrandMap_bed2bed, True)
            #     .run(p.normalizeCounts_txt2txt, True)
                
            #     .branch()
            #         .run(p.addTreatmentAndStrand_txt2txt, True)
            #         .cat(p.mergeGeneCounts, True)
            #     .stop()
            # .stop()    

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 120, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS')
                    # .cat(p.mergeTSSTES, True)
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch(True)
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-4000"', 120, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES')
                    .cat(p.mergeTSSTES, False)
                .stop()
            .stop() 

            .branch(False)
                .run(p.chromosomeCounts_bed2txt, True)
            .stop()

        .stop()

    # TAIR10
        .branch(True)
            .run(p.bowtie_fastq2sam, False, 'TAIR10')
            .run(p.convertToBam_sam2bam, False)
            .run(p.convertToBed_bam2bed, False)
            .run(p.uniqueSort_bed2bed, False)
            


            .branch(False)
                .run(p.lengthDistribution_bed2csv, True)
                .run(p.addTreatment_txt2txt, True)
                .cat(p.mergeLength, True)                
            .stop()

            .branch(False)
                .run(p.get27mer_bed2bed, True)
                .run(p.convertBedToFasta_bed2fa, True)
                .run(p.getNucleotideAbundanceTable_fa2csv, True)
                .cat(p.mergeNucleotide, True)
            .stop()

            .branch(False)
                .run(p.getAllLengths_bed2bed, False)
                .run(p.convertBedToFasta_bed2fa, False)
                .run(p.getNucleotideAbundanceTable_fa2csv, False)
                .run(p.addLength_csv2csv, True)
                .cat(p.mergeNucleotideFrequenciesAll, True)
            .stop()


        # Get BigWig Files
            .branch(False)
                .run(p.splitByStrand_bed2bed, True)
                
                .branch(True)
                    .run(p.convertToBedGraph_bed2bdg, True)
                    .run(p.toBigWig_bdg2bw, True)
                .stop()
            .stop()

        # Annotated transcript counts
            .branch(False)
                .run(p.geneStrandMap_bed2txt, True, "transcripts")
                .run(p.normalizeCounts_txt2txt, True)
                .run(p.addTSNTS_txt2txt, True)
                .cat(p.mergeGeneCounts, True, 'merged_annotatedTranscriptCounts.txt')                     
            .stop()

        # Bin count
            .branch(False)
                .run(p.binCount_bed2txt, True, 1000)
                .run(p.normalizeCounts_txt2txt, True)
                .run(p.addTreatment_txt2txt, True)
                .cat(p.mergeBinCounts, True, 'merged_bin1000.txt')
            .stop()

        # Measured transcript counts
            .branch(False)
                .run(p.geneStrandMap_bed2txt, True, "transcripts2")
                .run(p.normalizeCounts_txt2txt, True)
                .run(p.addTSNTS_txt2txt, True)
                .cat(p.mergeGeneCounts, True, 'merged_measuredTranscriptCounts.txt')                     
            .stop()

        # Coding gene counts
            .branch(True)
                .run(p.geneStrandMap_bed2txt, False, "genes")
                .run(p.normalizeCounts_txt2txt, False)
                .run(p.addTSNTS_txt2txt, False)
                    
                    .branch(True)
                        .run(p.GEO_txt2txt, True)
                        .run(p.adjustColumns_txt2txt, True)
                    .stop()
                
                .cat(p.mergeGeneCounts, False, 'merged_genes.txt')                     
            .stop()

        # DNase counts
            .branch(False)
                .run(p.dnaseCount_bed2txt, True)
                .run(p.binnedCountsToPositions_txt2txt, True, '"-1000"', 200, 10)
                .run(p.addTreatment_txt2txt, True)            
                .cat(p.mergeDNaseCounts, False, 'merged_dnaseCounts.txt')                                 
            .stop()

        # CHMM Global Repair only
            .branch(False)
                .run(p.subtractTCR_bed2bed, False)            
                .run(p.addArtificialGlobalRepair_bed2bed, False)
                .run(p.countChromatinStates_bed2txt, False)
                .run(p.normalizeCounts_txt2txt, False)
                .run(p.addTreatment_txt2txt, False)
                .cat(p.mergeChromtainStates, False, 'merged_noTCRchromatinStates.txt')                     
            .stop()

        # CHMM all
            .branch(False)
                .run(p.countChromatinStates_bed2txt, True)
                .run(p.normalizeCounts_txt2txt, True)
                .run(p.addTreatment_txt2txt, True)
                .cat(p.mergeChromtainStates, True, 'merged_chromatinStates.txt')                     
            .stop()

        # TSS and TES with expression quartiles
            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, False, 'TSS')
                .run(p.normalizeCounts_txt2txt, False)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS', 'Q')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS_Q1')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS', 'Q1')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS_Q2')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS', 'Q2')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS_Q3')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS', 'Q3')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TSS_Q4')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-8000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TSS', 'Q4')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, False, 'TES')
                .run(p.normalizeCounts_txt2txt, False)
                
                .branch(False)
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES', 'Q')
                    # .cat(p.mergeTSSTES, True, '_TAIR10Nasc')
                .stop()
            .stop() 

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TES_Q1')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES', 'Q1')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TES_Q2')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES', 'Q2')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TES_Q3')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES', 'Q3')
                .stop()
            .stop()

            .branch(False)
                .run(p.transcriptSiteCount_bed2txt, True, 'TES_Q4')
                .run(p.normalizeCounts_txt2txt, True)
                
                .branch()
                    .run(p.binnedCountsToPositions_txt2txt, True, '"-2000"', 100, 100)
                    .run(p.addTSNTS_txt2txt, True, 'TES', 'Q4')
                    .cat(p.mergeTSSTES, False, '_TAIR10Nasc_Q')                
                .stop()
            .stop()
        
        # TCR Analysis

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False)
                .run(p.addValue_txt2txt, False, 'real')
            .stop()

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {'random':False})
                .run(p.addValue_txt2txt, False, 'random')
            .stop()

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {"slice":False, "n":4, "sliceNum":1, "keyword":'_Q1'})
                .run(p.addValue_txt2txt, False, 'Q1')
            .stop()

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {"slice":False, "n":4, "sliceNum":2, "keyword":'_Q2'})
                .run(p.addValue_txt2txt, False, 'Q2')
            .stop()

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {"slice":False, "n":4, "sliceNum":3, "keyword":'_Q3'})
                .run(p.addValue_txt2txt, False, 'Q3')
            .stop()

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {"slice":False, "n":4, "sliceNum":4, "keyword":'_Q4'})
                .run(p.addValue_txt2txt, False, 'Q4')
            .stop()

        # Non overlapping transcripts

            .branch(False)
                .run(p.transcriptIntersect_bed2txt, False, {'transcripts': 'transcriptsNoNeighbor5K', 'keyword': '_no_neighbor'})
                .run(p.addValue_txt2txt, False, 'no_neighbor')
                .cat(p.mergeTCR, False)
            .stop()

        # DNase Analysis

            .branch(False)
                .run(p.dnaseIntersect_bed2txt, True)
                .run(p.addValue_txt2txt, True, 'real')
            .stop()

            .branch(False)
                .run(p.dnaseIntersect_bed2txt, True, {'random': True})
                .run(p.addValue_txt2txt, True, 'random')
                .cat(p.mergeDNase, True)
            .stop()
                
        # CHMM36 Analysis
            .branch(False)
                .run(p.chmmIntersect_bed2txt, True)
                .run(p.addValue_txt2txt, True, 'real')
            .stop()

            .branch(False)
                .run(p.chmmIntersect_bed2txt, True, {'random': True, 'bedA': 'chmm36_shuffled'})
                .run(p.addValue_txt2txt, True, 'random')
                .cat(p.mergeChmm, False)
            .stop()
        
        # Epigenetic Markers
            .branch(False)
                .run(p.epigeneticIntersect_bed2txt, False, {'bedA': 'epigeneticMarkers'})
                .run(p.addValue_txt2txt, False, 'real')
                .run(p.addTreatment_txt2txt, False)
            .stop()

            .branch(False)
                .run(p.epigeneticIntersect_bed2txt, False, {'random': False, 'bedA': 'epigeneticMarkers_shuffled'})
                .run(p.addValue_txt2txt, False, 'random')
                .run(p.addTreatment_txt2txt, False)                
                .cat(p.mergeEpigenetic, False)
            .stop()

        # Origin of replication
            .branch(False)
                .run(p.OORintersect_bed2txt, True, {'bedA': 'originOfReplication'})
                .run(p.addValue_txt2txt, True, 'real')
                .run(p.addTreatment_txt2txt, True)
            .stop()

            .branch(False)
                .run(p.OORintersect_bed2txt, True, {'random': True, 'bedA': 'originOfReplication_shuffled'})
                .run(p.addValue_txt2txt, True, 'random')
                .run(p.addTreatment_txt2txt, True)                
                .cat(p.mergeOOR, False)
            .stop()

        # CCA1 potential gene counts
            .branch(False)
                .run(p.geneStrandMap_bed2txt, True, "CCA1potentialGenes")
                .run(p.normalizeCounts_txt2txt, True)
                .run(p.addTSNTS_txt2txt, True)
                .cat(p.mergeGeneCounts, False, 'merged_CCA1potGenes.txt')                     
            .stop()

        .stop()
    )
###########################################################
