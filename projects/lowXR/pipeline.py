import os    
import sys
import argument
import bed
import pipeTools
import argparse
import generalUtils
from glob import glob
sys.path.append('..')
from referenceGenomePath import referenceGenomePath
from commonPipelineMethods import commonPipeline

class pipeline(commonPipeline):
    def __init__(self, input, args = argument.args()):
        commonPipeline.__init__(self, input, args)
        OUTPUT_DIR = '1802'
        SAMPLE_STAT_FILE = 'samples.csv'    
        self.input = os.path.join(os.path.curdir, 'data', 'raw', self.input)
        self.outputDir = os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR)
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.sampleBase = os.path.join(self.outputDir, os.path.splitext(os.path.basename(self.input))[0])
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.attributes = sorted(sampleDictionary.keys())
        self.saveInput([self.input])
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.title + '.txt',
            '-e ': 'err_' + self.title + '.txt',
        }
        self.wmParams = self.defaultWmParams
        self.reportCount = 0
        self.paths = referenceGenomePath()
        self.referenceGenome = 'hg19'
        self.fragmentLengths = [10, 15, 20, 25]       

    def fullPath2wildcard(self, fullPath, prefix = ''):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join([prefix + "*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        # code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "cat " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        print(code)
        os.system(code)

    def cutadapt_fastq2fastq(self):
        self.saveOutput([os.path.join(self.outputDir, os.path.basename(self.output[0]))])
        codeList = [
            'cutadapt'
            ]
        if self.method == 'loXR':
            # fivePrimeAdaptor = 'GTTCAGAGTTCTACAGTCCGACGATCTTCGCTACCTTAG'
            fivePrimeAdaptor = 'NTCGCTACCTTAG'
            threePrimeAdaptor = 'GACCGTTATAGTTATGGAATTCTCGGGTGCCAAGG'
            linkedAdapter = fivePrimeAdaptor + '...' + threePrimeAdaptor
            codeList += [
                '-a', linkedAdapter
                # '-e', 0.2
            ]
        elif self.method == 'XR':
            threePrimeAdaptor = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG'
            codeList += [
                '-a', threePrimeAdaptor
            ]
        else:
            raise ValueError('The method ' + self.method +  ' is not defined.')

        codeList += [
            '-o', self.output,
            self.input
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

    def getDinucleotideAbundanceTable_fa2csv(self):
        return self.getNucleotideAbundanceTable_fa2csv(2)

    def mergeNucleotide(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes)
        output = os.path.join(self.outputDir, '..', 'merged_nucleotide.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeLength(self):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['length', 'value'] + sorted(self.attributes)
        output = os.path.join(self.outputDir, '..', 'merged_length.txt')
        self.catFiles(wildcard, headers, output)
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

    def splitByStrand_bed2bed(self):
        strand = ['"\\t\+$"', '"\\t\-$"']
        output = [
            self.addExtraWord(self.output[0], '_Plus'), 
            self.addExtraWord(self.output[0], '_Minus')
        ]
        self.saveOutput(output)
        self.saveOutput(self.prettyOutput())
        codeList = [
            'grep -P ',
            strand,
            self.input[0],
            '>',
            self.output
        ]
        self.execM(codeList)
        return self

    def prettyOutput(self, additionalWord = ''):
        newOutputs = []
        for o in self.output:
            if 'Plus' in o:
                extraWord = '_Plus'
            elif 'Minus' in o:
                extraWord = '_Minus'
            else:
                extraWord = ''
            extension = pipeTools.getExtension(o)
            newOutputs.append(os.path.join(os.path.dirname(o),self.title + extraWord + additionalWord + '.' + extension))
        return newOutputs

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
        self.saveOutput(self.prettyOutput())
        codeList = [
            'bedGraphToBigWig',
            self.input,
            self.reference['limits'],
            self.output
        ]
        self.execM(codeList)
        return self

    def bowtie2_fastq2sam(self, referenceGenome = "hg19"):
        noCpus = 8
        self.reference = self.paths.get(referenceGenome)
        self.mutateWmParams({'-n ': str(noCpus)})
        output = [self.output[0].replace(".sam", self.reference['name'] + ".sam")]
        self.saveOutput(output)
        codeList = [
            'bowtie2',
            '--quiet',
            '-N', 1,
            '-p', noCpus,
            '-x', self.reference['bowtie2'],
            '-U', self.input,
            '-S', self.output
        ]
        self.execM(codeList)
        return self
   
    def convertToBam_sam2bam(self):
        codeList = [
            'samtools',
            'view',
            '-b',
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
        # print('self.output')
        # print(self.output)
        # print(self.input)
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

    def getAllLengths_bed2bed(self):
        self.addWordsToOutput(self.fragmentLengths)
        codeList = [
            'bed2getCertainIntervalLengths.py',
            '-i', self.input * len(self.fragmentLengths),
            '-o', self.output,
            '-l', self.fragmentLengths
        ]
        self.execM(codeList)
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

    def mergeNucleotideFrequenciesAll(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0])
        wildcard = wildcard.replace("geAlLe" + str(self.fragmentLengths[0]), "*")
        headers = ['position', 'sequence', 'value'] + sorted(self.attributes) + ['fragmentLength']
        output = os.path.join(self.outputDir, '..', 'merged_NucleotideFrequencies' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def countLine(self, denominator, labelList):
        labelNoList = []
        reportNo = self.reportCount
        for lebel in labelList:
            reportNo += 1            
            labelNoList.append(reportNo)
        self.reportCount = reportNo
        def getCount(filename):
            return sum(1 for line in open(filename))/denominator
        
        def getReportName(filename):
            return os.path.splitext(filename)[0] + '_report.txt'

        output = pipeTools.listOperation(getReportName, self.output)
        if self.runMode and self.runFlag:
            countList = pipeTools.listOperation(getCount, self.output)
        else:
            countList = len(self.output) * [0]
        codeList = [
            'echo',
            '-e',
            labelNoList,
            '\'\t\'',            
            countList,
            '\'\t\'',
            labelList,
            '>',
            output
        ]
        self.execM(codeList)
        return self

    def combineReport(self):
        wildcard = self.sampleBase + '*_report.txt'
        print(wildcard)
        output = self.sampleBase + '_combinedReport.txt'
        codeList = [
            'cat',
            wildcard,
            '|',
            'addColumns.py',
            '-c', ' '.join(self.list2attributes(self.attributes)),
            '>',
            output
        ]
        print(codeList)
        self.execM(codeList)
        return self

    def mergeReport(self):
        wildcard = self.sampleBase.replace(os.path.basename(self.sampleBase), '*') + '_combinedReport.txt'
        headers = ['step', 'count', 'name' ] + sorted(self.attributes)
        output = os.path.join(self.outputDir, '..', 'merged_report.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeBinCounts(self, outputName = 'mergedBinCounts.txt'):
        wildcard = self.fullPath2wildcard(self.input[0])        
        headers = ['chromosome', 'start', 'end', 'gene', 'score', 'strand', 'count'] + sorted(self.attributes)
        output = os.path.join(outputName)
        self.catFiles(wildcard, headers, output)
        return self

    def getCertainLengths_bed2bed(self):
        self.readLength = [10, 29]
        fragmentLengths = list(range(self.readLength[0], self.readLength[1]))
        codeList = [
            'bed2getCertainIntervalLengths.py',
            '-i', self.input,
            '-o', self.output,
            '-l', ' '.join(str(e) for e in fragmentLengths)
        ]
        self.execM(codeList)
        return self

    def binMap_bed2txt(self):
        codeList = [
            'bedtools makewindows',
            '-g', self.reference['limits'],
            '-w', 10000,
            '|',
            'addColumns.py',
            '-c . . .',
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
    
    def normalizeCountsBasedOnOriginal_txt2txt(self):
        try:
            totalMappedReads = int(open(self.totalMappedReadsFile, 'r').readline())
        except:
            if self.runFlag == True and self.runMode == True:
                raise ValueError('no read count file')
            totalMappedReads = 1000000
       
        self.scaleFactor = float(1000000)/totalMappedReads
        
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
    
    def subsample_bed2bed(self, readNumber):
        codeList = [
            'subsample',
            '-n', readNumber,
            '--seed', 123,
            self.input,
            '>', self.output
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

    parser_cat = subparsers.add_parser('report', help='report help')
    parser_cat.add_argument('-n', required= False, default="1", help='input index')
    parser_cat.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
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

###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    args = getArgs()
    inputIndex = args.get("n")
    input = getInputFromIndex(inputIndex)

    p = pipeline(input, args)
    (p
        # .branch(False)
        #     .run(p.GEO_fastq2txt, False)
        # .stop()
        
        .run(p.cutadapt_fastq2fastq, False)
            .report(p.countLine, False, 4, ['After the adaptor cut'])

        .branch(True)
            .run(p.bowtie2_fastq2sam, False, 'hg19')
            .run(p.convertToBam_sam2bam, False)
            .run(p.convertToBed_bam2bed, False)
                .report(p.countLine, False, 1, ['Aligned reads'])
            .run(p.sort_bed2bed, False, {'unique': True})
                .report(p.countLine, False, 1, ['Unique reads'])

            .branch(True)
                .run(p.getCertainLengths_bed2bed, False)
                .run(p.sort_bed2bed, False)
                
                .run(p.subsample_bed2bed, True, 6500000)
                
                .branch(True)
                    .run(p.writeTotalMappedReads_bed2txt, True)
                .stop()
                
                .run(p.binMap_bed2txt, True)
                .run(p.normalizeCountsBasedOnOriginal_txt2txt, True)
                .run(p.addTreatment_txt2txt, True)
                .cat(p.mergeBinCounts, True, 'mergedBins.txt') 
            .stop()

            .branch(False)
                .run(p.writeTotalMappedReads_bed2txt, False)
            .stop()

            .branch(False)
                .run(p.lengthDistribution_bed2csv, False)
                .run(p.addTreatment_txt2txt, False)
                .cat(p.mergeLength, False)                
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

                .branch(False)
                    .run(p.getNucleotideAbundanceTable_fa2csv, True)
                    .run(p.addLength_csv2csv, True)
                    .cat(p.mergeNucleotideFrequenciesAll, True)
                .stop()

                .branch(False)
                    .run(p.getDinucleotideAbundanceTable_fa2csv, True)
                    .run(p.addLength_csv2csv, True)
                    .cat(p.mergeNucleotideFrequenciesAll, True, '_diNuc')
                .stop()

            .stop()

            .branch(False)
                .run(p.splitByStrand_bed2bed, False)

                # Get BigWig Files
                .branch(False)
                    .run(p.convertToBedGraph_bed2bdg, True)
                    .run(p.toBigWig_bdg2bw, True)
                .stop()
            .stop()

            .branch(False)
                .run(p.convertToBedGraph_bed2bdg, True)
                .run(p.toBigWig_bdg2bw, True)
            .stop()

        .stop()
        .report(p.combineReport, False)
        .cat(p.mergeReport, False)
    )
###########################################################

