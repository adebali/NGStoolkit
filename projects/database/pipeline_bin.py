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
import slurm
from configHelpers import samples, paths

def getArgs():
    parser = argparse.ArgumentParser(description='Standard XR-seq and Damage-seq Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-samples', default='sample.json', required= False, help='file for sample configuration in JSON format')
    parser_run.add_argument('-reference', default='reference.json', required= False, help='file for reference paths in JSON format')
    parser_run.add_argument('-e', required= True, type=int, help='experiment index')
    parser_run.add_argument('-cpu', required= False, type=int, default=1, help='requested cpu number')
    parser_run.add_argument('-n', required= True, type=int, help='input index')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-samples', default='sample.json', required= False, help='file for sample configuration in JSON format')
    parser_cat.add_argument('-reference', default='reference.json', required= False, help='file for reference paths in JSON format')
    parser_cat.add_argument('-e', required= True, type=int, help='experiment index')
    parser_cat.add_argument('-cpu', required= False, type=int, default=1, help='requested cpu number')    
    parser_cat.add_argument('-n', required= False, type=int, default=1, help='input index')
    parser_cat.add_argument('--mock', required= False, default=True, action='store_true', help='mock flag')
    parser_cat.add_argument('--noPrint', required= False, default=True, action='store_true', help='prints no code when stated')

    parser_batch = subparsers.add_parser('batch', help='batch help')
    parser_batch.add_argument('-samples', default='sample.json', required= False, help='file for sample configuration in JSON format')
    parser_batch.add_argument('-reference', default='reference.json', required= False, help='file for reference paths in JSON format')
    parser_batch.add_argument('-e', required= True, type=int, help='experiment index')
    parser_batch.add_argument('-cpu', required= False, type=int, default=1, help='requested cpu number')    
    parser_batch.add_argument('-mem', required= False, type=int, default=16000, help='requested memory (MB)')    
    parser_batch.add_argument('-s', required= False, type=int, nargs='+', help='space separated sample list')
    parser_batch.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_batch.add_argument('--noPrint', required= False, default=True, action='store_true', help='prints no code when stated')

    args = parser.parse_args()
    return argument.args(args)

args = getArgs()

sampleClass = samples(args.get('samples'))
pathClass = paths(args.get('reference'))

class pipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        # information = json.load(open(input))
        information = sampleClass.key2attributes(input)
        for key in information.keys():
            setattr(self, key, information[key])
        self.outputDir = 'data'
        self.saveInput(os.path.join(self.outputDir, self.SRA_id + '.txt'))
        self.reference = pathClass.getGenome(self.genome)
        self.noCpus = args.get('cpu')
        self.attributes = self.list2attributes(self.columns)
        self.attributesString = list((str(x) for x in self.attributes))

    def fastqdump_txt2fastq(self):
        self.saveOutput([os.path.join(self.outputDir, self.SRA_id + '.fastq')])
        codeList = ['rm', '-f', self.output]
        for run in self.runs:
            codeList += ['&&', 'fastq-dump', '--stdout', run]
        codeList += ['>>', self.output]

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

    def addTreatment_txt2txt(self):
        columns = self.list2attributes(self.columns)
        columnStringList = ' '.join(str(e) for e in columns)
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', columnStringList
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
            codeList = ['cut -f 7', f, '|', 'paste -s', '>>', outputFile]
            os.system(' '.join(codeList))
        return self
    
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

    def mergeBinCounts(self, outputName = 'mergedBinCounts.txt'):
        wildcard = self.fullPath2wildcard(self.input[0])        
        headers = ['chromosome', 'start', 'end', 'gene', 'score', 'strand', 'count'] + self.columns
        output = os.path.join(outputName)
        self.catFiles(wildcard, headers, output)
        return self


###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__" and args.get("subprogram") == "batch":
    os.system('mkdir -p log')
    sampleClass = samples(args.get('samples'))
    if args.get('s'):
        sampleNoList = args.get('s')
    else:
        sampleNoList = sampleClass.experiment2sampleNoList(args.get('e'))

    pipelineParameters = {
        "--job-name=": "db_pipeline",
        "-n ": args.get('cpu'),
        "--mem=": args.get('mem'),
        "--time=": "1-00:00:00",
        "--output=": "./log/%A_%a.out",
        "--error=": "./log/%A_%a.err",
        "--array=": ','.join(str(x) for x in sampleNoList),
        "--mail-type=": "END,FAIL",      # notifications for job done & fail
        "--mail-user=": "oadebali@gmail.com" # send-to address
    }

    job = slurm.Slurm('python ' + os.path.basename(__file__) + ' run -e ' + str(args.get('e')) + ' -n $SLURM_ARRAY_TASK_ID' + ' -cpu ' + str(args.get('cpu')))
    job.assignParams(pipelineParameters)
    job.printScript()
    if not args.get('mock'):
        jobId = job.run()


if __name__ == "__main__" and (args.get("subprogram") in ["run", "cat"]):
    input = sampleClass.experimentNoAndNo2id(args.get('e'), args.get('n'))

    p = pipeline(input, args)
    (p
        .run(p.fastqdump_txt2fastq, False)
        .run(p.cutadapt_fastq2fastq, False)
    )

    if p.method == 'RNA-seq':
        p.run(p.tophat_fastq2bam, False)
    else:
        p.run(p.bowtie2_fastq2sam, False)
        p.run(p.convertToBam_sam2bam, False)
    
    p.run(p.convertToBed_bam2bed, False)
    p.run(p.sort_bed2bed, False, {'unique': p.sortUnique})

    if p.method == 'XR-seq':
        (p
            .run(p.getCertainLengths_bed2bed, False)
            .run(p.sort_bed2bed, False)
        )  
    elif p.method == 'Damage-seq':
        (p
            .run(p.slopBed_bed2bed, False)
            .run(p.convertToFixedRange_bed2bed, False)
            .run(p.sort_bed2bed, False)
            .run(p.convertBedToFasta_bed2fa, False)
            .run(p.getDamageSites_fa2bed, False)
            .run(p.sort_bed2bed, False)
        )

    elif p.method == 'CPD-seq':
        (p
            .run(p.slopBed_bed2bed, False)
            .run(p.convertToFixedRange_bed2bed, False)
            .run(p.sort_bed2bed, False)
            .run(p.convertBedToFasta_bed2fa, False)
            .run(p.getDamageSites_fa2bed, False)
            .run(p.sort_bed2bed, False)
        )

    (p
        .branch(False)
            .run(p.writeTotalMappedReads_bed2txt, False)
        .stop()

        .branch(True)
            .run(p.binMap_bed2txt, False)
            .run(p.normalizeCountsBasedOnOriginal_txt2txt, False)
            .run(p.addTreatment_txt2txt, True)
            .cat(p.mergeBinCounts, True, 'mergedBins.txt')            
        .stop()

        .branch(False)
            .run(p.geneStrandMap_bed2txt, False)
            .run(p.normalizeCountsBasedOnOriginal_txt2txt, False)
            .cat(p.mergeGeneCounts, False)
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
