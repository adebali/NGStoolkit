import os
import sys
import pipeTools
import generalUtils
import bed
import argparse
import argument
sys.path.append('../..')
from referenceGenomePath import referenceGenomePath
sys.path.append('..')
from parent import pipeline as parentpipe
import tempfile

class pipeline(parentpipe):
    def __init__(self, input, args = argument.args()):
        parentpipe.__init__(self, input, args)
        OUTPUT_DIR = '1803'
        SAMPLE_STAT_FILE = 'samples.csv'    
        self.input = os.path.join(os.path.curdir, 'data', 'input')
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.attributes = sorted(sampleDictionary.keys())
        # self.saveInput([self.input])
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.sample + '.txt',
            '-e ': 'err_' + self.sample + '.txt',
        }
        self.wmParams = self.defaultWmParams
        self.paths = referenceGenomePath()
        self.referenceGenome = 'hg19'
        self.reference = self.paths.get(self.referenceGenome)
        self.reference['whiteListBed'] = '../wgEncodeDukeMapabilityRegionsExcludable.bed'

    def dump_2fastq(self):
        output = [self.sample + '.fastq']
        self.saveOutput([os.path.join(os.path.realpath(self.outputDir), os.path.basename(output[0]))])
        codeList = [
            'fastq-dump',
            self.SRA_id,
            '-Z',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def align_fastq2bam(self):
        temp1_, temp1 = tempfile.mkstemp()
        self.picardScript = '/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar'

        codeList = [
        'cutadapt -f fastq -e 0.1 -q 20 -O 3 -a CTGTCTCTTATACACATCT', self.input,
        '> /dev/stdout ',
        '| bwa mem -M ',
        self.reference['bwa'],
        '/dev/stdin',
        '| samtools view -S -u -F2304 -q 10 -L ', self.reference['whiteListBed'],
        '-',
        '| samtools sort - ', temp1, '&&',
        'java -Xmx5g -jar', self.picardScript,
        'MarkDuplicates I=', temp1, 'O=', self.output, 'M=' + self.sample + '.md.txt'
        ]
        self.execM(codeList)
        return self

    def cutadapt_fastq2fastq(self):
        self.execM(['cutadapt -f fastq -e 0.1 -q 20 -O 3 -a CTGTCTCTTATACACATCT', self.input, '>', self.output])
        return self

    def bwa_fastq2bam(self):
        self.execM(['bwa mem -M ', self.reference['bwa'], self.input, '>', self.output])
        return self

    def whitelist_bam2bam(self):
        self.execM(['samtools view -S -u -F2304 -q 10 -L ', self.reference['whiteListBed'], self.input, '>', self.output])
        return self

    def sort_bam2bam(self):
        self.execM(['samtools sort', self.input, '>', self.output])
        return self

    def markDuplicates_bam2bam(self):
        self.picardScript = '/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar'        
        self.execM(['java -Xmx5g -jar', self.picardScript, 'MarkDuplicates I=', self.input, 'O=', self.output, 'M=' + self.sample + '.md.txt'])
        return self

    def filter_bam2bam(self):
        self.execM(['samtools view -h -F1024', self.input, '| grep -v -P \'\\tchrM\\t\' | samtools view -b - >', self.output])
        return self
    
    def callPeaks_bam2txt(self):
        self.execM(['macs2 callpeak --keep-dup all -t', self.input, '-n', self.sample])
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

###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    args = getArgs()
    inputIndex = args.get("n")
    input = getInputFromIndex(inputIndex)

    p = pipeline(input, args)
    (p
        .run(p.dump_2fastq, False)
        
        .branch(False)
            .run(p.align_fastq2bam, True)
        .stop()

        .branch(True)
            .run(p.cutadapt_fastq2fastq, False)
            .run(p.bwa_fastq2bam, False)
            .run(p.whitelist_bam2bam, False)
            .run(p.sort_bam2bam, False)
            .run(p.markDuplicates_bam2bam, False)

            .branch(True)
                .run(p.filter_bam2bam, True)
                .run(p.callPeaks_bam2txt, True)
            .stop()
        .stop()
    )