import os    
import sys
import argument
import bed
import pipeTools
sys.path.append('..')
from commonPipelineMethods import commonPipeline

class pipeline(commonPipeline):
    def __init__(self, input, args = argument.args()):
        commonPipeline.__init__(self, input, args)
    
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