import os    
import sys
import argument
import bed
import pipeTools
sys.path.append('..')
from commonPipelineMethods import commonPipeline
import tempfile

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

    def convertToBedGraphTwoPolars_bed2bdg(self):
        if self.runFlag and self.runMode:
            scaleFactor = float(1000000)/self.internalRun(bed.bed(self.input[0]).getHitNum, [], self.runFlag, 'get hit number')
        else:
            scaleFactor = 1
        scaleFactorList = [scaleFactor * 1, scaleFactor * -1]
        codeList = [
            'bedtools',
            'genomecov',
            '-i', self.input,
            '-g', self.reference['limits'],
            '-bg',
            '-scale', scaleFactorList,
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
    
    def bg4Intersect_bed2txt(self, args={}):
        bedAfileName = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'BG4/GSE76688_rhh_25cyc_BG4_12082015_peaks_summit1k_bed6.bed')
        keyword = args.get('keyword', '')
        output = [self.addExtraWord(self.output[0], keyword)]
        self.saveOutput(output)
        
        sliceFlag = args.get('slice', False)
        if sliceFlag:
            numberOfSlices = args.get('n')
            currentSlice = args.get('sliceNum')
            keyword = args.get('keyword')
            sliceTemp_, sliceTemp = tempfile.mkstemp()
            codeList = [
                'bed2sliceByScore.py',
                '-i', bedAfileName,
                '-n', numberOfSlices,
                '-slice', currentSlice,
                '>', sliceTemp
            ]
            self.execM(codeList)
            bedAfileName = sliceTemp

        random = args.get('random', False)
        if random:
            output = [
                self.addExtraWord(self.output[0], '_random'),
            ]
            self.saveOutput(output)

        averageLength = totalRecord = totalMappedReads = perNmappedReads= 1
        if self.runMode == True and self.runFlag == True:
            inputBed = bed.bed(self.input[0])
            totalMappedReads = inputBed.getHitNum()
            bedA = bed.bed(bedAfileName)
            averageLength = bedA.getAverageLength()
            totalRecord = bedA.getHitNum()
            perNmappedReads = 1000000
        
        shuffleCode = [
            'bedtools',
            'shuffle',
            '-i', bedAfileName,
            '-g', self.reference['limits'],
            '|',
            'sort',
            '-k1,1',
            '-k2,2n',
            '-k3,3n'
        ]
        if random:
            codeList = list(shuffleCode)
        else:
            codeList = [
                'cat',
                bedAfileName
            ]

        codeList += [
            '|',
            'bedtools',
            'intersect',
            '-a', 'stdin',
            '-b', self.input,
            '-wa',
            '-wb',
            '-F', 0.50,
            '|',
            'bedIntersect2positionTxt.py',
            '|',
            'bedIntersectPositionCount.py',
            '-count', 7,
            '-cat', 8,
            '-ends', 'double',
            '-scale',
            1/float(totalRecord),
            100/float(averageLength),
            perNmappedReads/float(totalMappedReads),
            '-o', self.output
        ]
        self.execM(codeList)
        
        _temp, temp = tempfile.mkstemp()

        windowLength = 20
        oneSideFlankingLength = 2000

        if random:
            prepareAbedCodeList = list(shuffleCode)
        else:
            prepareAbedCodeList = [
                'cat',
                bedAfileName
            ]
        prepareAbedCodeList += [
            '|',
            'bed2updownstream.py',
            '--fixed',
            '-l', oneSideFlankingLength,
            '|',
            'bed2removeChromosomeEdges.py',
            '--fixed',
            '-l', 10000,
            '-g', self.reference['limits'],
            '>', temp
        ]
        self.execM(prepareAbedCodeList)

        totalRecord = 1
        if self.runMode == True and self.runFlag == True:
            bedA = bed.bed(temp)
            totalRecord = bedA.getHitNum()/2

        flankingCodeList = [
            'bedtools',
            'intersect',
            # '-a', 'stdin',
            '-a', temp,
            '-b', self.input,
            '-wa',
            '-wb',
            '-F', 0.50,
            '|',
            'bedIntersect2positionTxt.py',
            '--flanking',
            '--fixed',
            '-w', windowLength,
            '|',
            'bedIntersectPositionCount.py',
            '-count', 7,
            '-cat', 8,
            '-ends', 'remove',
            '-scale', 
            1/float(totalRecord), 
            1/float(windowLength),
            perNmappedReads/float(totalMappedReads),
            '>>', self.output,
            '&&',
            'rm', temp
        ]
        self.execM(flankingCodeList)
        return self

    def addTreatment_txt2txt(self):
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
    
    def mergeBG4(self, extraWord=''):
        wildcard = self.fullPath2wildcard(self.input[0])
        headers = ['position', 'count', 'strand'] + self.attributes
        output = os.path.join(self.outputDir, '..', 'merged_BG4' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self    