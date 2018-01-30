import os    
from pipe import pipe
import argument
import bed
import pipeTools

class pipeline(pipe):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
    
    def intersect10K_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["10Kbed"],
            '-b', self.input,
            '-wa',
            '-c',
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def filterChrPos_txt2txt(self, args):
        codeList = [
            'awk', 
            '\'$1 == "' + args['chromosome'] + '" && ' +
            '$2 >= ' + str(args['startGT']) + ' && ' +
            '$2 <= ' + str(args['startLT']) + '\'',
            self.input,
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

    def merge10KCounts(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_Plus", "_*s")
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes + ['strandName']
        output = os.path.join(self.outputDir, '..', 'merged_10KCounts' + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
        return self

    def mergeLeadLag(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_Plus", "_*s")
        headers = ['position', 'count'] + self.attributes + ['strandName']
        output = os.path.join(self.outputDir, '..', 'merged_leadLag' + self.leadLagExtraWord + extraWord + '.txt')
        self.catFiles(wildcard, headers, output)
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
        codeList = [
            'bedGraphToBigWig',
            self.input,
            self.reference['limits'],
            self.output
        ]
        self.execM(codeList)
        return self

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
        self.finalBed = self.output[0]
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

    def leadLag_bed2txt(self, args={}):
        zone = args.get('zone', 'RIZ')
        distance = args.get('distance', 1000000)
        scoreCutoff = args.get('score', 900)
        self.leadLagExtraWord = '_' + zone + '_' + str(round(float(distance)/1000000,4)) + 'MB' + '_c' + str(scoreCutoff)
        self.saveOutput(pipeTools.listOperation(self.addExtraWord, self.output, self.leadLagExtraWord))
        # self.saveOutput(self.lineBa)
        if self.runFlag and self.runMode:
            self.scaleFactor = float(1000000)/self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number')
        else:
            self.scaleFactor = 1
        codeList = [
            'python ../LeadLagStrandProfiles.py',
            '-i', self.input,
            '-a', self.reference['HeLaS3'][zone],
            '-s', self.scaleFactor,
            '-g', self.reference['limits'],
            '-d', distance,
            '-cutoff', scoreCutoff,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    # remove the close ones (neighborDistance)
    # get only 900 to 1000 (awk command: colNo, sign, value)
    # center and flanking (exists)
    # separate into windows (makeWindows)
    # intersect separate strands (bedtools)
    # count bins (custom or exists)