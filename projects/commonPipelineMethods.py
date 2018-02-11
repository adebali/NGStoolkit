import os    
from pipe import pipe
import argument
import pipeTools

class commonPipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
    
    def sort_bed2bed(self, args={}):
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