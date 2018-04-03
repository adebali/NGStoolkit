import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
from glob import glob
import argparse
import argument
import json
import bed
sys.path.append('../../..')
from referenceGenomePath import referenceGenomePath
from sampleHelpers import samples

sampleClass = samples('sample.json')

# def id2attribute(id, inputJson="sample.json"):
#     def recursiveBase(master, d):
#         if d.get('base', False):
#             baseD = master[d['base']]
#             d_updatedWithBase = recursiveBase(master, baseD)
#             d_updatedWithBase.update(d)
#             return d_updatedWithBase
#         return d
#     sheet = json.load(open(inputJson))
#     sampleDict = sheet[id]
#     completeSampleDict = recursiveBase(sheet, sampleDict)
#     return completeSampleDict

# print(id2attribute('SRX997094'))

class pipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        # information = json.load(open(input))
        information = sampleClass.key2attributes(input)
        for key in information.keys():
            setattr(self, key, information[key])
        self.outputDir = 'data'
        self.saveInput(os.path.join(self.outputDir, self.SRA_id + '.txt'))
        self.paths = referenceGenomePath()
        self.reference = self.paths.get(self.genome)
        self.noCpus = 8
        self.categories = []
        for e in self.columns:
            self.categories.append(getattr(self, e))

    def fastqdump_txt2fastq(self):
        self.saveOutput([os.path.join(self.outputDir, self.SRA_id + '.fastq')])
        codeList = ['rm', '-f', self.output]
        for run in self.runs:
            codeList += ['&&', 'fastq-dump', '--stdout', run]
        codeList += ['>>', self.output]

        self.execM(codeList)
        return self

    def addTreatment_txt2txt(self):
        columns = self.categories
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', ' '.join(str(x) for x in self.categories)
        ]
        self.execM(codeList)
        return self

    def cutadapt_fastq2fastq(self):
        inputBaseNames = []
        for inputFile in self.input:
            inputBaseNames.append(os.path.basename(inputFile))
        if not (self.adapter3 or self.adapter5):
            codeList = ['ln', '-sf', inputBaseNames, self.output]
        else:
            if self.method == 'XR-seq' or self.method == 'CPD-seq':
                codeList = [
                    'cutadapt',
                    '-a', ' -a '.join(self.adapter3),
                    '-o', self.output,
                    self.input
                ]
            elif self.method == 'Damage-seq':
                codeList = [
                    'cutadapt',
                    '--discard-trimmed', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
                    '-g', ' -g '.join(self.adapter5), # The adapter is located at the 5' end
                    '-o', self.output,
                    self.input
                ]
            # elif self.method == 'CPD-seq':
            #     codeList = [
            #         'cutadapt',
            #         '-a', ' -a '.join(self.adapter3),
            #         '-o', self.output,
            #         self.input
            #     ]
        self.execM(codeList)
        return self

    def bowtie2_fastq2sam(self):
        codeList = [
            'bowtie2',
            '--quiet',
            '-N', 1,
            '-p', self.noCpus,
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
        if self.method == 'RNA-seq':
            self.saveInput([os.path.join(self.tophatDirectory, 'accepted_hits.bam')])
            self.output = [os.path.join(self.outputDir, self.title + '.bed')]
            self.saveOutput(self.output)
                    
        codeList = [
            'bedtools',
            'bamtobed',
            '-i', self.input,
            '>', self.output
        ]
        self.finalBed = self.output
        self.execM(codeList)
        return self

    def sort_bed2bed(self, args={}):
        '''sort bed (or a related) file based on the first 3 columns:
        chromosome(string) start(integer) end(integer)
        {'unique': True} dictionary input will remove the identical rows based on the three columns
        '''
        uniqueFlag = ''
        if args.get('unique', False):
            uniqueFlag = '-u'
        codeList = [
            'sort',
            uniqueFlag,
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        self.finalSortedBed = self.output
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

    def getDamageSites_fa2bed(self):
        if not self.method in ['Damage-seq', 'CPD-seq']:
            raise ValueError('This function cannot be applied to the method ' + self.method)
        if self.damage in ['UV_CPD', 'UV_64']:
            self.motifRegex = '\'.{4}(T|t|C|c)(T|t|C|c).{4}\''
        elif self.damage in ['Cisplatin']:
            self.motifRegex = '\'.{4}(G|g)(G|g).{4}\''
            

        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', self.motifRegex 
        ]
        self.execM(codeList)
        return self
    
    def getCertainLengths_bed2bed(self):
        fragmentLengths = list(range(self.readLength[0], self.readLength[1]))
        codeList = [
            'bed2getCertainIntervalLengths.py',
            '-i', self.input,
            '-o', self.output,
            '-l', ' '.join(str(e) for e in fragmentLengths)
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

    def geneStrandMap_bed2txt(self):
        if self.stranded:
            extensions = ['_TS', '_NTS']
            if self.method == "XR-seq" or self.method == "Damage-seq":
                strandParameters = ['-S', '-s']
            elif self.method == "CPD-seq":
                strandParameters = ['-s', '-S']
            newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
            self.saveOutput(newOutput)
            self.input = [self.input[0], self.input[0]]
        else:
            strandParameters = ['']

        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference['genes'],
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
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

    def tophat_fastq2bam(self):
        newOutput = [self.input[0] + '_1']

        inputBaseNames = []
        for inputFile in self.input:
            inputBaseNames.append(os.path.basename(inputFile))
        
        codeList = [
            'ln -fs',
            inputBaseNames,
            newOutput
        ]

        self.execM(codeList)
        self.saveInput(newOutput)

        self.tophatDirectory = os.path.join(self.outputDir, 'tophat', self.title)

        codeList = [
            'mkdir -p',
            self.tophatDirectory,
            '&&',
            'tophat',
            '-p', self.noCpus,
            '-o', self.tophatDirectory,
            self.reference['bowtie2'],
            self.input,
        ]
        self.execM(codeList)
        return self

    def mergeGeneCounts(self, outputFile = 'mergedGeneCounts.txt'):
        for f in self.input:
            codeList = ['cat ', f, '>>', outputFile]
            os.system(' '.join(codeList))
        return self

def getArgs():
    parser = argparse.ArgumentParser(description='Standard XR-seq and Damage-seq Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-e', required= True, type=int, help='experiment index')
    parser_run.add_argument('-n', required= True, type=int, help='input index')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-e', required= True, type=int, help='experiment index')
    parser_cat.add_argument('-n', required= False, type=int, default=1, help='input index')
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
# inputIndex = args.get("n")
# input = getInputFromIndex(inputIndex)

# input = "SRX997094"
# def experimentAndNo2id(experiment, no):

# input = "SRX220485"
input = sampleClass.experimentNoAndNo2id(args.get('e'), args.get('n'))

###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    p = pipeline(input, args)
    (p
        .run(p.fastqdump_txt2fastq, True)
        .run(p.cutadapt_fastq2fastq, True)
    )

    if p.method == 'RNA-seq':
        p.run(p.tophat_fastq2bam, True)

    else:
        p.run(p.bowtie2_fastq2sam, True)
        p.run(p.convertToBam_sam2bam, True)
    
    p.run(p.convertToBed_bam2bed, True)
    p.run(p.sort_bed2bed, True, {'unique': p.sortUnique})

    if p.method == 'XR-seq':
        (p
            .run(p.getCertainLengths_bed2bed, True)
            .run(p.sort_bed2bed, True)
        )  
    elif p.method == 'Damage-seq':
        (p
            .run(p.slopBed_bed2bed, True)
            .run(p.convertToFixedRange_bed2bed, True)
            .run(p.sort_bed2bed, True)
            .run(p.convertBedToFasta_bed2fa, True)
            .run(p.getDamageSites_fa2bed, True)
            .run(p.sort_bed2bed, True)
        )

    elif p.method == 'CPD-seq':
        (p
            .run(p.slopBed_bed2bed, True)
            .run(p.convertToFixedRange_bed2bed, True)
            .run(p.sort_bed2bed, True)
            .run(p.convertBedToFasta_bed2fa, True)
            .run(p.getDamageSites_fa2bed, True)
            .run(p.sort_bed2bed, True)
        )

    (p
        .branch(True)
            .run(p.writeTotalMappedReads_bed2txt, True)
        .stop()

        .branch(True)
            .run(p.geneStrandMap_bed2txt, True)
            .run(p.normalizeCountsBasedOnOriginal_txt2txt, True)
            .run(p.addTreatment_txt2txt, True)
            .cat(p.mergeGeneCounts, True)
        .stop()
    )
 
    if p.stranded:
        (p
            .branch(False)
                .run(p.splitByStrand_bed2bed, False)
                
                .branch(False)
                    .run(p.convertToBedGraph_bed2bdg, False)
                    .run(p.toBigWig_bdg2bw, False)
                .stop()
            .stop()
        )
    else:
        (p
            .branch(False)
                .run(p.convertToBedGraph_bed2bdg, False)
                .run(p.toBigWig_bdg2bw, False)
            .stop()
        )

###########################################################
