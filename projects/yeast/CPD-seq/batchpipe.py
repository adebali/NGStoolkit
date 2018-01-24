import os
from pipe import pipe
import pipeTools
from glob import glob
import slurm
import argparse
import bed

defaultDir = 'dataDir/1801'
defaultExtension = '.wig.gz'
outDir = defaultDir

parser = argparse.ArgumentParser(description='pipe')
parser.add_argument('-i', default=defaultDir, help='input directory')
parser.add_argument('-o', required=False, default=outDir, help='output directory')
parser.add_argument('-extension', required=False, default=defaultExtension, help='extension eg .gz')
parser.add_argument('-mode', required=False, default='slurm', choices=['slurm', 'regular'], help='mode')
parser.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')


args = parser.parse_args()
files = glob(os.path.join(os.path.realpath(args.i), '*' + args.extension))
outputDir = os.path.realpath(args.o)
os.system('mkdir -p ' + outputDir)
os.system('rm err_*')
os.system('rm log_*')

class pipeline(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        # if args.mock == False:
            # self.runMode = True
        self.input = input
        self.outputDir = outputDir
        self.wmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'pipe',
            '-o ': 'log_' + os.path.basename(self.input).split('.')[0] + '.txt',
            '-e ': 'err_' + os.path.basename(self.input).split('.')[0] + '.txt'
        }
        if args.mode == 'slurm':
            self.execM = self.execMwm

    def toBed_gz2bed(self):
        codeList = [
            'zcat',
            self.input,
            '|',
            'wig2bed',
            '>', self.output
        ] 
        self.execM(codeList)
        return self

    def expand_bed2bed(self):
        codeList = [
            'bedScore2expandedBed.py',
            '-i',
            self.input,
            '-placeHolder', 1,
            '>', self.output
        ] 
        self.execM(codeList)
        self.finalBed = self.output
        return self

    def count_bed2txt(self):
        codeList = [
            'bedtools intersect',
            '-a', '/proj/sancarlb/users/ogun/seq/S288C/S288C_reference_genome_R64-2-1_20150113/genes.bed',
            '-b',
            self.input,
            '-wa',
            '-c',
            '>',
            self.output
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

for file in files:
    p = pipeline(file)
    (p
        .run(p.toBed_gz2bed, False)
        .run(p.expand_bed2bed, False)
        .run(p.count_bed2txt, False)
        .run(p.normalizeCounts_txt2txt, True)
    )

