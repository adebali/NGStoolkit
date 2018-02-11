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

    def intersectRepChmm_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.reference["HeLaS3"]["replicationChmm"],
            '-b', self.input,
            '-wa',
            '-c',
            '-F', 0.50,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def repChmmTotal_txt2txt(self):
        codeList = [
            'Rscript',
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 'repChmmTotal.R'),
            self.input,
            self.output
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

    def mergeRepChmm(self, extraWord = ''):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_Plus", "_*s")
        headers = ['RD', 'CH', 'length', 'RPM', 'RPKM'] + self.attributes + ['strandName']
        output = os.path.join(self.outputDir, '..', 'merged_repChmm' + extraWord + '.txt')
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
        self.saveOutput(self.prettyOutput())
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

    def G1bPeaks_bed2txt(self, args={}):
        zone = args.get('zone', 'RIZ')
        distance = args.get('distance', 100000)
        scoreCutoff = args.get('score', 500)
        self.leadLagExtraWord = '_' + zone + '_' + str(round(float(distance)/1000000,4)) + 'MB' + '_c' + str(scoreCutoff)
        self.saveOutput(pipeTools.listOperation(self.addExtraWord, self.output, self.leadLagExtraWord))
        # self.saveOutput(self.lineBa)
        if self.runFlag and self.runMode:
            self.scaleFactor = float(1000000)/self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number')
        else:
            self.scaleFactor = 1
        codeList = [
            'python ../G1b_peaks.py',
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

    def leadLagSelected_bed2txt(self, args={}):
        zone = args.get('zone', 'selectedRIZ')
        distance = args.get('distance', 100000)
        self.leadLagExtraWord = '_' + zone + '_' + str(round(float(distance)/1000000,4)) + 'MB'
        self.saveOutput(pipeTools.listOperation(self.addExtraWord, self.output, self.leadLagExtraWord))
        # self.saveOutput(self.lineBa)
        if self.runFlag and self.runMode:
            self.scaleFactor = float(1000000)/self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number')
        else:
            self.scaleFactor = 1
        codeList = [
            'python ../LeadLagStrandProfiles_selected.py',
            '-i', self.input,
            '-a', self.reference['HeLaS3'][zone],
            '-s', self.scaleFactor,
            '-g', self.reference['limits'],
            '-d', distance,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

    def replicationDomain_bed2txt(self, args={}):
        zone = args.get('zone', 'replicationDomains')
        keyword = args.get('keyword', 'DTZ')
        distance = args.get('distance', 100000)
        self.leadLagExtraWord = '_' + zone + '_' + keyword + '_' + str(round(float(distance)/1000000,4)) + 'MB'
        self.saveOutput(pipeTools.listOperation(self.addExtraWord, self.output, self.leadLagExtraWord))
        # self.saveOutput(self.lineBa)
        if self.runFlag and self.runMode:
            self.scaleFactor = float(1000000)/self.internalRun(bed.bed(self.finalBed).getHitNum, [], self.runFlag, 'get hit number')
        else:
            self.scaleFactor = 1
        codeList = [
            'python ../replicationDomains.py',
            '-i', self.input,
            '-a', self.reference['HeLaS3'][zone],
            '-s', self.scaleFactor,
            '-keyword', keyword,
            '-g', self.reference['limits'],
            '-d', distance,
            '-o', self.output
        ]
        self.execM(codeList)
        return self


    def strandAsymmetry_txt2bdg(self):
        import generalUtils
        def strandAssymmetry(x, y, arg=None):
            xl = x.strip().split('\t')
            yl = y.strip().split('\t')
            xval = float(xl[6])
            yval = float(yl[6])
            
            if xval + yval > 0:
                value = (xval - yval)/(xval + yval)
            else:
                return False
            xl[3] = str(value)
            return '\t'.join(xl[0:4])

        self.saveOutput(self.prettyOutput('_asymmetry'))
        self.internalRun(
            generalUtils.lineBasedTwoFilesOperation,
            [self.input[0], self.input[1], self.output[0].replace('_Plus', ''), strandAssymmetry, [None]], self.runFlag, 
            'strand Assymmetry'
            )
        return self

    def strandRatio_txt2bdg(self):
        import generalUtils
        def strandRatio(x, y, arg=None):
            xl = x.strip().split('\t')
            yl = y.strip().split('\t')
            xval = float(xl[6])
            yval = float(yl[6])
            
            if yval > 0 and xval>0:
                value = float(xval)/yval
            else:
                value = '.'
            xl[3] = str(value)
            return '\t'.join(xl[0:4])

        input1 = self.input[0]
        input2 = self.input[1]
        self.input = [self.input[0]]
        self.saveOutput([self.prettyOutput('_strandRatio')[0].replace('_Plus', '')])
        output = self.output[0]
        self.internalRun(
            generalUtils.lineBasedTwoFilesOperation,
            [input1, input2, output, strandRatio, [None]], self.runFlag, 
            'strand ratio'
            )
        
        return self

    # remove the close ones (neighborDistance)
    # get only 900 to 1000 (awk command: colNo, sign, value)
    # center and flanking (exists)
    # separate into windows (makeWindows)
    # intersect separate strands (bedtools)
    # count bins (custom or exists)