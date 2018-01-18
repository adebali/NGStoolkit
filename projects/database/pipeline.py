import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
from glob import glob
import argparse
import argument
sys.path.append('..')
from referenceGenomePath import referenceGenomePath



class pipeline(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        self.SRA_id = self.input
        self.saveInput(self.SRA_id + '.txt')
       
    def fastqdump_txt2fastq(self):
        codeList = [
            'fastq-dump',
            self.SRA_id
        ]
        self.saveOutput(self.SRA_id + '.fastq')
        self.execM(codeList)
        return self


def getArgs():
    parser = argparse.ArgumentParser(description='Standard XR-seq and Damage-seq Pipeline', prog="pipeline.py")
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

input = "SRX997094"

###########################################################
#  Pipeline
###########################################################
if __name__ == "__main__":
    p = pipeline(input)
    (p
        .run(p.fastqdump_txt2fastq, True)
    )
###########################################################
