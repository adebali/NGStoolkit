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
sys.path.append('..')
from parent import pipeline as parentpipe
from referenceGenomePath import referenceGenomePath
import numpy
import tempfile

class pipeline(parentpipe):
    def __init__(self, input, args = argument.args()):
        print(input)
        pipe.__init__(self, input, args)
        OUTPUT_DIR = '1801'
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
            '--job-name=': 'Damage-seq',
            '-o ': 'log_' + self.title + '.txt',
            '-e ': 'err_' + self.title + '.txt',
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
            newOutputs.append(os.path.join(os.path.dirname(o),self.title + extraWord + '.' + extension))
        return newOutputs

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


    def cutadapt_fastq2fastq(self):
        output = [os.path.join(self.outputDir, os.path.basename(self.output[0]))]
        adapter = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'
        codeList = [
            'cutadapt',
            '--discard-trimmed', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
            '-g', adapter, # The adapter is located at the 5' end
            '-o', output,
            self.input,
        ]
        self.execM(codeList)
        self.saveOutput(output)
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

    def plotNucleotideAbundance_csv2pdf(self):
        codeList = [
            'plotNucleotideAbundance.r',
            self.input,
            self.title
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

    # def getAllLengths_bed2bed(self):
    #     self.addWordsToOutput(self.fragmentLengths)
    #     print(self.output)
    #     codeList = [
    #         'bed2getCertainIntervalLengths.py',
    #         '-i', self.input * len(self.fragmentLengths),
    #         '-o', self.output,
    #         '-l', self.fragmentLengths
    #     ]
    #     self.execM(codeList)
    #     return self

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

    def addTreatmentAndPlusMinus_txt2txt(self):
        columns = self.list2attributes(self.attributes)
        columnStringList = [' '.join(columns + ['plus']), ' '.join(columns + ['minus'])]
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
        self.motifRegex = '\'.{4}(T|t|C|c)(T|t|C|c).{4}\'' # Get GG only at the positions 5 and 6.

        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', self.motifRegex 
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


def getArgs():
    parser = argparse.ArgumentParser(description='Damage-seq Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-n', required= True, help='input index')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-n', required= False, default="13", help='input index')
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

###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    args = getArgs()
    inputIndex = args.get("n")
    input = getInputFromIndex(inputIndex)

    p = pipeline(input, args)
    (p
        .branch(False)
            .run(p.GEO_fastq2txt, False)
        .stop()
        
        .run(p.cutadapt_fastq2fastq, False)

        .branch(True)
            .run(p.bowtie2_fastq2sam, False, 'hg19')
            .run(p.convertToBam_sam2bam, False)
            .run(p.convertToBed_bam2bed, False)
            .run(p.sort_bed2bed, False)
            .run(p.slopBed_bed2bed, False)
            .run(p.convertToFixedRange_bed2bed, False)
            .run(p.sort_bed2bed, False)
            .run(p.convertBedToFasta_bed2fa, False)

            .branch(False) # Plot nucleotide abundance
                .run(p.getNucleotideAbundanceTable_fa2csv, True)
                .run(p.addTreatment_csv2txt, True)
                .cat(p.mergeNucleotideAbundance, True)
            .stop()

            .branch(False) # Plot dinucleotide abundance
                .run(p.getDimerAbundanceTable_fa2csv, True)
                .run(p.addTreatment_csv2txt, True)
                .cat(p.mergeNucleotideAbundance, True, '_diNuc')
            .stop()

            .run(p.getDamageSites_fa2bed, False)
            .run(p.sort_bed2bed, False)

            .branch(True)
                .run(p.writeTotalMappedReads_bed2txt, False)
            .stop()

            .branch(True)
                .run(p.splitByStrand_bed2bed, False)
                
                .branch(True)
                    .run(p.intersect10K_bed2txt, True)
                    .run(p.normalizeCountsBasedOnOriginal_txt2txt, True)
                    .run(p.addTreatmentAndPlusMinus_txt2txt, True)

                        .branch(True)
                            .run(p.filterChrPos_txt2txt, True, {'chromosome': 'chr2', 'startGT': 113000000, 'startLT': 130000000})
                            .cat(p.merge10KCounts, True, '_2_113_130')
                        .stop()
                    .cat(p.merge10KCounts, True)
                .stop()

                # Get BigWig Files
                .branch(False)
                    .run(p.convertToBedGraph_bed2bdg, True)
                    .run(p.toBigWig_bdg2bw, True)
                .stop()
            .stop()

        .branch(True)

        # Annotated transcript counts
            # .branch(False)
            #     .run(p.geneStrandMap_bed2txt, True, "transcripts")
            #     .run(p.normalizeCounts_txt2txt, True)
            #     .run(p.addTSNTS_txt2txt, True)
            #     .cat(p.mergeGeneCounts, True, 'merged_annotatedTranscriptCounts.txt')                     
            # .stop()
        .stop()
    )
###########################################################
