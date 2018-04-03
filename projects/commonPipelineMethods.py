import os    
from pipe import pipe
import argument
import pipeTools

class commonPipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
    
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

    def sort_bdg2bdg(self, args={}):
        return self.sort_bed2bed()