import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob

kHomeDir = os.path.join(os.path.expanduser("~"), "ogun")
kChromHmmDir = os.path.join(kHomeDir, "chromHMM")
chromHMMfiles = generalUtils.sorted_nicely(glob(os.path.join(kHomeDir, 'chromHMM', '*_chromHMM')))

class myPipe(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        outputDirName = '0101'
        self.outputDir = os.path.join(os.path.dirname(self.input), outputDirName)
        os.system('mkdir -p ' + self.outputDir)
        fileBaseName = os.path.basename(self.input)
        samples = generalUtils.table2dictionary(generalUtils.file('dataDir/samples_minReadCount.csv'), 'sample')
        sampleDictionary = samples[fileBaseName][0]
        self.treatment = sampleDictionary['treatment_title']
        self.minimumReadCount = round(int(sampleDictionary['minReadCount']) - 500000)
        self.group = sampleDictionary['group']
        self.sampleMinPyrCount = sampleDictionary['minPyrHitNo']
        #if '-group' in sys.argv:
        #    if self.group != int(sys.argv[sys.arg.index('-group') + 1]):
        #        self.runMode = False
        self.dnaseBed = '/nas02/home/a/d/adebali/ucsc/wgEncodeUwDnaseNHDFAd.fdr01peaks.hg19.bed'
        self.extra_name = 'ncl'
        self.bowtie_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19'
        self.fasta_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19.fa'
        self.chromosome_limits = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19.chrSizes.bed'
        #if runDir != False:
        #    self.outputDirectory = os.path.join(os.path.dirname(self.input), runDir)
        #    os.system('mkdir -p ' + self.outputDirectory)
        input1 = self.input
        input2 = generalUtils.in2out(self.input, '.1.fastq', '.2.fastq')
        self.saveInput([input1, input2])


    def cutadapt_fastq2fastq(self, runFlag=True):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput(output)
        adapter = 'GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT'

        codeList = [
            'cutadapt',
            '--discard-trimmed', # Remove the reads containing the adapter sequence. If a read contains adapter, it is because of the artifact, no damage is expected on this read.
            '-g', adapter, # The adapter is located at the 5' end
            '-o', self.output[0],
            '-p', self.output[1],
            self.input[0],
            self.input[1]
        ]
        self.execM(codeList)
        return self


    def bowtie_fastq2sam(self, runFlag=True):
        output = [pipeTools.in2out(self.input[0], '.fastq', '.bow.sam')]
        self.saveOutput(output)
        codeList = [
            'bowtie',
            '-t', self.bowtie_reference,
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            '-m', 4, # Do not report the reads that are mapped on to more than 4 genomic locations
            '-X', 1000,
            '--seed', 123, # Randomization parameter in bowtie,
            '-p', 4,
            '-1', self.input[0],
            '-2', self.input[1],
            self.output
        ]
        self.execM(codeList)
        return self

    def slopBed_bed2bed(self, runFlag=True):
        slopB = 6
        self.addWordToOutput('b6')
        codeList = [
            'bedtools',
            'slop',
            '-i', self.input,
            '-g', self.chromosome_limits, # Chomosomal lengths, needed in case region shift goes beyond the chromosomal limits.
            '-b', slopB,
            '-s', # Take strand into account
            '>', self.output
        ]
        self.execM(codeList)
        return self


    def convertBedToFasta_bed2fa(self):
        codeList = [
            'bedtools',
            'getfasta',
            '-fi', self.fasta_reference,
            '-bed', self.input,
            '-fo', self.output,
            '-s' # Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented
        ]
        self.execM(codeList)
        return self

    def getPyrimidineDimers_fa2bed(self):
        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', '\'.{4}[T|C][T|C].{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        ]
        self.execM(codeList)
        return self

    def sampleFromBed_bed2bed(self):
        codeList = [
            'sort', '-R', self.input,
            '|',
            '-head', self.sampleMinPyrCount,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getNucleotideAbundanceTable_fa2csv(self):
        nucleotideOrder = 'TCGA'
        codeList = [
            'fa2nucleotideAbundanceTable.py',
            '-i', self.input,
            '-o', self.output,
            '-n', nucleotideOrder,
            '--percentage'
        ]
        self.execM(codeList)
        return self

    def plotNucleotideAbundance_csv2pdf(self):
        codeList = [
            'plotNucleotideAbundance.r',
            self.input,
            self.treatment
        ]
        self.execM(codeList)
        return self

    def getDimerAbundanceTable_fa2csv(self):
        codeList = [
            'fa2kmerAbundanceTable.py',
            '-i', self.input,
            '-o', self.output,
            '-k', 2,
            '--percentage'
        ]
        self.execM(codeList)
        return self

    def plotDinucleotideAbundance_csv2pdf(self):
        codeList = [
            'plotNucleotideFreqLine.r',
            self.input,
            self.treatment
        ]
        self.execM(codeList)
        return self

    def convertToBam_sam2bam(self, runFlag=True):
        codeList = [
            'samtools',
            'view',
            '-bf', '0x2', #	each segment properly aligned according to the aligner
            #'-Sb'
            '-o',
            self.output,
            self.input
        ]
        self.execM(codeList)
        return self

    def convertToBed_bam2bedpe(self, runFlag=True):
        codeList = [
            'bedtools',
            'bamtobed',
            '-bedpe',
            '-mate1',
            '-i', self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def removeRedundancy_bedpe2bedpe(self, runFlag=True):
        codeList = [
            'bed2nonRedundantBed.py',
            '-i', self.input,
            '-o', self.output,
            '-c', '1,2,3,4,5,6,7,9' # columns
        ]
        self.execM(codeList)
        return self

    def uniqueSort_bedpe2bedpe(self, runFlag=True):
        codeList = [
            'sort',
            '-u',
            '-k1,1',
            '-k2,2n',
            # '-k3,3n',
            # '-k4,4',
            # '-k5,5',
            # '-k6,6',
            #'-k9,9',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def sortBed_bed2bed(self, runFlag=True):
        codeList = [
            'sort',
            '-k1,1',
            '-k2,2n',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self


    def convertBedpeToSingleFrame_bedpe2bed(self, runFlag=True):
        codeList = [
            'bedpe2bed.py',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertToFixedRange_bed2bed(self, runFlag=True):
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

    def convertToBedGraph_bed2bdg(self, runFlag=True):
        codeList = [
            'bedtools',
            'genomecov',
            '-i', self.input,
            '-g', self.chromosome_limits,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertToWig_bed2wig(self, runFlag=True):
        codeList = [
            'igvtools',
            'count',
            self.input,
            self.output,
            self.chromosome_limits
        ]
        self.execM(codeList)
        return self

    def convertToBigBed_bed2bb(self, runFlag=True):
        codeList = [
            'bedToBigBed',
            self.input,
            self.chromosome_limits,
            self.output
        ]
        self.execM(codeList)
        return self

    def convertToBigWig_bdg2bw(self, runFlag=True):
        codeList = [
            'wigToBigWig',
            '-clip',
            self.input,
            self.chromosome_limits,
            self.output
        ]
        self.execM(codeList)
        return self

    def addHeader_bdg2bdg(self, runFlag=True):
        codeList = [
            'echo "track type=bedGraph name=' + self.treatment + '"', '>', self.output,
            '&&',
            'tail -n +2', self.input, '>>', self.output
        ]
        self.execM(codeList)
        return self

    def coverageChromHMM_bed2bed(self):
        input = self.input[0]
        self.totalMappedReads = self.internalRun(bed.bed(input).getHitNum, [], self.runFlag)
        outputs = []
        os.system("mkdir -p " + os.path.join(os.path.dirname(input), "chromHMM"))
        for chromHMMfile in chromHMMfiles:
            regionName = chromHMMfile.split('/')[-1].split('_NHFL_')[0]
            os.system("mkdir -p " + os.path.join(os.path.dirname(input), "chromHMM", self.treatment))
            outputDir = os.path.join(os.path.dirname(input), "chromHMM", self.treatment)
            output = os.path.join(outputDir, regionName + '.bed')
            codeList = [
                'bedtools',
                'coverage',
                '-counts',
                '-a', chromHMMfile, # Primary genomic locations
                '-b', self.input, # Read file, to get the read counts overlapping with the primary genomic locations
                '>', output
            ]
            self.execM(codeList)
            outputs.append(output)
        self.saveOutput(outputs)
        return self

    def normalizeCoverageChromHMM_bed2bed(self):
        multiplicationFactor = float(10000000)/self.totalMappedReads if self.totalMappedReads else 1 # Normalizes counts by total mapped read number
        self.addWordToOutput('n1K')
        codeList = [
            'bedCount2normalizedCount.py',
            '-i', self.input,
            '-c', 7,
            '-o', self.output,
            '-l', 1000,
            '-m', multiplicationFactor
        ]
        self.execM(codeList)
        return self

    def appendNull_bed2bed(self):
        nullValue = "NA"
        maxEntry = 115306
        codeList = [
            'bedCount2appendedVoid.py',
            '-i', self.input,
            '-o', self.output,
            '-n', maxEntry,
            '-null', 'NaN'
        ]
        self.execM(codeList)
        return self

    def getCountColumn_bed2txt(self):
        columnNo = 7
        codeList = [
            'cat', self.input,
            '|',
            'cut', '-f', columnNo,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def mergeChromHMM_txt2txt(self):
        self.saveOutput(os.path.join(self.outputDir, self.treatment + "_mergedChromHMM.txt"))
        print(self.input)
        codeList = ['paste'] + self.input + ['>' + self.output[0]]
        self.execM(codeList)
        return self

    def plotChromHMM_txt2pdf(self):
        codeList = [
            "plotChromHMM.r",
            self.input,
            self.treatment
        ]
        self.execM(codeList)
        return self

    def dnaseClosest_bed2bed(self):
        codeList = [
            'bedtools',
            'closest',
            '-a', self.dnaseBed,
            '-b', self.input,
            '-k', 1,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def retrieveOnlyTheClosest_bed2bed(self):
        codeList = [
            'sort',
            '-u',
            '-k1,1',
            '-k2,1n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getDistanceFromClosest_bed2txt(self):
        codeList = [
            'bedClosest2distance.py',
            '-i', self.input,
            '-c1', 2,
            '-c2', 12,
            '-o', self.output
        ]
        self.execM(codeList)
        return self



p = myPipe(sys.argv[1])
(p
    .run(p.cutadapt_fastq2fastq, True)
    .run(p.bowtie_fastq2sam, True)
    .run(p.convertToBam_sam2bam, True)
    .run(p.convertToBed_bam2bedpe, True)
    .run(p.uniqueSort_bedpe2bedpe, True)
    .run(p.convertBedpeToSingleFrame_bedpe2bed, True)
    .run(p.slopBed_bed2bed, True)
    .run(p.convertToFixedRange_bed2bed, True)
    .run(p.sortBed_bed2bed, True)
    .run(p.convertBedToFasta_bed2fa, True)
        
        .branch()
        .run(p.getNucleotideAbundanceTable_fa2csv, False)
        .run(p.plotNucleotideAbundance_csv2pdf, False)
        .stop()

        .branch()
        .run(p.getDimerAbundanceTable_fa2csv, False)
        .run(p.plotDinucleotideAbundance_csv2pdf, False)
        .stop()

    .run(p.getPyrimidineDimers_fa2bed, True)
    .run(p.sampleFromBed_bed2bed, True)
        
        .branch()
        .run(p.coverageChromHMM_bed2bed, False)
        .run(p.normalizeCoverageChromHMM_bed2bed, False)
        .run(p.appendNull_bed2bed, False)
        .run(p.getCountColumn_bed2txt, False)
        .run(p.mergeChromHMM_txt2txt, False)
        .run(p.plotChromHMM_txt2pdf, False)
        .stop()
        
        .branch(False)
        .run(p.dnaseClosest_bed2bed, False)
        .run(p.retrieveOnlyTheClosest_bed2bed, False)
        .run(p.getDistanceFromClosest_bed2txt, False)
        .stop()

    .run(p.convertToBigBed_bed2bb, False)
)