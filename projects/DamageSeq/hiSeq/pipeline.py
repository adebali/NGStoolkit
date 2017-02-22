import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob
import argparse

# Required tools in path
# bowtie
# bedtools
# samtools
# subsample

SAMPLE_STAT_FILE = 'dataDir/samples.csv'

parser = argparse.ArgumentParser(description='DamageSeq Pipeline')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-n', required= False, help='output')
parser.add_argument('--mock', required= False, action='store_true', help='mock flag')
parser.add_argument('--noPrint', required= False, action='store_true', help='prints no code when stated')
args = parser.parse_args()
inputFile = args.i
inputIndex = args.n

if inputFile == None and inputIndex == None:
    raise ValueError("No input or index was given! Exiting ...")
    sys.exit()
elif inputFile != None:
    fileBaseName = os.path.basename(inputFile)
    samples = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'sample')
    sampleDictionary = samples[fileBaseName][0]
    input = generalUtils.file(inputFile)
else:
    indexes = generalUtils.table2dictionary(generalUtils.file(SAMPLE_STAT_FILE), 'no')
    sampleDictionary = indexes[inputIndex][0]
    input = generalUtils.file(os.path.realpath(os.path.join('dataDir', 'raw', sampleDictionary['sample'])))


class myPipe(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        outputDirName = '0106'
        self.totalMappedReads = None
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', outputDirName))
        os.system('mkdir -p ' + self.outputDir)
        self.treatment = sampleDictionary['treatment_title']
        # self.minimumReadCount = round(int(sampleDictionary['minReadCount']) - 500000)
        self.group = sampleDictionary['group']
        self.cell = sampleDictionary['cell']
        # self.sampleMinPyrCount = sampleDictionary['minPyrHitNo']
        self.sampleMinPyrCount = None
        self.product = sampleDictionary['product']
        if int(self.group) == 1:
            self.sampleMinPyrCount = 7700464
        if int(self.group) == 2:
            self.sampleMinPyrCount = 1833986

        # if self.product == 'CPD' or self.product == '_6-4':
        #     self.motifRegex = '\'.{4}(T|C|t|c)(T|C|t|c).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        # elif self.product == 'GG':
        #     self.motifRegex = '\'.{4}(G|g)(G|g).{4}\'' # Get GG only at the positions 5 and 6.
        # else:
        #     raise ValueError("Unknown product type. Exiting...")

        if self.product == 'CPD':
            self.motifRegex = '\'.{4}(TT|tt|CT|ct|cT|Ct).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        elif self.product == '_6-4':
            self.motifRegex = '\'.{4}(TT|tt|TC|tc|Tc|tC).{4}\'' # Get sequences pyrimidine dimer at the positions 5 and 6.
        elif self.product == 'GG':
            self.motifRegex = '\'.{4}(G|g)(G|g).{4}\'' # Get GG only at the positions 5 and 6.
        elif self.product == 'NA':
            self.motifRegex = '\'.{10}\'' # Get GG only at the positions 5 and 6.
        else:
            raise ValueError("Unknown product type. Exiting...")

        ## Transcription factor binding sites
        if self.cell == "GM12878":
            # self.txnBedList = ['/proj/sancarlb/users/ogun/ENCODE/Txn/wgEncodeRegTfbsClusteredWithCellsV3.GM12878.sorted.1K.sorted.bed']
            self.txnBedList = ['/proj/sancarlb/users/ogun/ENCODE/Txn/wgEncodeRegTfbsClusteredWithCellsV3.GM12878.sorted.1K.sorted.top100K.bed']
            self.dnaseBedList = ['/proj/sancarlb/users/ogun/ENCODE/dnase/GM12878.bed']
        elif self.cell == "NHF1":
            self.txnBedList = ['/proj/sancarlb/users/ogun/ENCODE/Txn/wgEncodeRegTfbsClusteredWithCellsV3.BJ.sorted.1K.sorted.bed']
            self.dnaseBedList = ['/proj/sancarlb/users/ogun/ENCODE/dnase/BJ.bed']
        else:
            raise ValueError("No such a cell is defined: " + self.cell)

        #if '-group' in sys.argv:
        #    if self.group != int(sys.argv[sys.arg.index('-group') + 1]):
        #        self.runMode = False

        self.dnaseBed = '/nas02/home/a/d/adebali/ucsc/wgEncodeUwDnaseNHDFAd.fdr01peaks.hg19.bed'
        self.extra_name = 'ncl'
        self.bowtie_reference = None
        self.fasta_reference = None
        self.chromosome_limits = None

        self.nucleosomeBed = None
        # self.nucleosomeBedRep1Top = '/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep1.t1M.bed'
        # self.nucleosomeBedRep1Top = '/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep1.all.bed'
        self.nucleosomeBedRep1Top = None

        self.nucleosomeBedList = None

        # self.nucleosomeBedRep1Top = '/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep1.bed'
        self.nucleosomeTopBed = None
        self.nucleosomeLowBed = None
        self.void147bed = None

        self.TSSbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min10000_TSS_win100.bed'
        self.TESbedFile = '/proj/sancarlb/users/ogun/seq/hg19/annotations/expression/hg19_expression_noHeader_sgt300_noNeighIn6000_min10000_TES_win100.bed'

        if self.cell == 'GM12878':
            self.chmmCellLine = self.cell
            self.fifteenChromatinStates = '/proj/sancarlb/users/ogun/ENCODE/chromatinStates/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed'
        elif self.cell == 'NHF1':
            self.chmmCellLine = 'NHLF'
            self.fifteenChromatinStates = '/proj/sancarlb/users/ogun/ENCODE/chromatinStates/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed'
        # self.fifteenChromatinStates = os.path.join('/nas/longleaf/home/adebali/ogun/chromHMM/statebyline', self.chmmCellLine, 'states.bed')
        #if runDir != False:
        #    self.outputDirectory = os.path.join(os.path.dirname(self.input), runDir)
        #    os.system('mkdir -p ' + self.outputDirectory)
        input1 = self.input
        input2 = generalUtils.in2out(self.input, '.1.fastq', '.2.fastq')
        self.saveInput([input1, input2])

    def changeDefaultValues(self, key):
        if key == "hg19":
            self.referenceNickname = 'hg19'
            output = pipeTools.listOperation(self.addExtraWord, self.output, '.hg19')
            self.saveOutput(output)
            self.bowtie_reference = '/proj/seq/data/HG19_UCSC/Sequence/BowtieIndex/genome'
            self.fasta_reference = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.fa')
            self.chromosome_limits = generalUtils.file('/nas/longleaf/home/adebali/ogun/seq/hg19/HG19_UCSC/genome.chrSizes.bed')
            kHomeDir = os.path.join(os.path.expanduser("~"), "ogun")
            kChromHmmDir = os.path.join(kHomeDir, "chromHMM")
            if p.cell == "NHF1":
                self.chromHMMfiles = generalUtils.sorted_nicely(glob(os.path.join(kHomeDir, 'chromHMM', 'w1000', '*_chromHMM')))
                self.maxChromHMMhitNumber = 19166369
                self.nucleosomeBedList = None
                self.CTCF = None
                self.allBoundTFBSlist = None
                self.stat3 = None
            elif p.cell == "GM12878":
                self.chromHMMfiles = None
                self.maxChromHMMhitNumber = None
                # self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/nucleosomeFromSheera.bed']
                self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/nucleosome/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.100K.1Kfl.bed']
                self.CTCF = ["/proj/sancarlb/users/ogun/ENCODE/Txn/wgEncodeRegTfbsClusteredWithCellsV3.CTCF.GM12878.bed.fa.motif.sorted.slop20.bed"]
                self.allBoundTFBSlist = [
                    "/nas/longleaf/home/adebali/ogun/ENCODE/Txn/motifIntersectedChip-seq.hittingMotifs.motifIntervals.bed",
                    "/nas/longleaf/home/adebali/ogun/ENCODE/Txn/motifIntersectedChip-seq.hittingMotifs.motifIntervals.upstream.bed",
                    "/nas/longleaf/home/adebali/ogun/ENCODE/Txn/motifIntersectedChip-seq.hittingMotifs.motifIntervals.downstream.bed",
                ]
                self.STAT3 = ["/nas/longleaf/home/adebali/ogun/ENCODE/Txn/STAT3_ENCFF001VFM.bed"]
        elif key == "hg19_nucleosome":
            self.referenceNickname = 'hg19nuc'
            output = pipeTools.listOperation(self.addExtraWord, self.output, '.hg19nuc')
            self.saveOutput(output)
            self.bowtie_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19'
            self.fasta_reference = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19_2.fa'
            self.chromosome_limits = '/nas02/home/a/d/adebali/ogun/ENCODE/hg19/nucleosome/hg19.chrSizes.bed'
            # self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep' + str(i + 1) + '.bed' for i in range(1, 9)]
            # self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/top100K/Nucleosome_Gm12878_rep' + str(i + 1) + '.bed' for i in range(1, 9)]
            self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep1.sga.bed.1K.intersect.selected.146.bed']
            self.nucleosomeBedList = ['/proj/sancarlb/users/ogun/ENCODE/Nucleosome_Gm12878_rep1.sga.500flanking.bed']
            self.void147bed = generalUtils.file('/proj/sancarlb/users/ogun/ENCODE/hg19/nucleosome/hg19.147.bed')
        return self
    
    def cutadapt_fastq2fastq(self):
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


    def bowtie_fastq2sam(self):
        output = [self.output[0]]
        self.saveOutput([self.addExtraWord(output[0], '.' + self.referenceNickname)])
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



    def slopBed_bed2bed(self):
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
            '-r', self.motifRegex 
        ]
        self.execM(codeList)
        return self

    def sampleFromBed_bed2bed(self):
        codeList = [
            'subsample',
            '-n', self.sampleMinPyrCount,
            '--seed', 123,
            self.input,
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

    def convertToBam_sam2bam(self):
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

    def convertToBed_bam2bedpe(self):
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

    def removeRedundancy_bedpe2bedpe(self):
        codeList = [
            'bed2nonRedundantBed.py',
            '-i', self.input,
            '-o', self.output,
            '-c', '1,2,3,4,5,6,7,9' # columns
        ]
        self.execM(codeList)
        return self

    def uniqueSort_bedpe2bedpe(self):
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

    def sortBed_bed2bed(self):
        codeList = [
            'sort',
            '-k1,1',
            '-k2,2n',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self


    def convertBedpeToSingleFrame_bedpe2bed(self):
        codeList = [
            'bedpe2bed.py',
            self.input,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertToFixedRange_bed2bed(self):
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

    def convertToBedGraph_bed2bdg(self):
        codeList = [
            'bedtools',
            'genomecov',
            '-i', self.input,
            '-g', self.chromosome_limits,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def convertToWig_bed2wig(self):
        codeList = [
            'igvtools',
            'count',
            self.input,
            self.output,
            self.chromosome_limits
        ]
        self.execM(codeList)
        return self

    def convertToBigBed_bed2bb(self):
        codeList = [
            'bedToBigBed',
            self.input,
            self.chromosome_limits,
            self.output
        ]
        self.execM(codeList)
        return self

    def convertToBigWig_bdg2bw(self):
        codeList = [
            'wigToBigWig',
            '-clip',
            self.input,
            self.chromosome_limits,
            self.output
        ]
        self.execM(codeList)
        return self

    def addHeader_bdg2bdg(self):
        codeList = [
            'echo "track type=bedGraph name=' + self.treatment + '"', '>', self.output,
            '&&',
            'tail -n +2', self.input, '>>', self.output
        ]
        self.execM(codeList)
        return self

    def coverageChromHMM_bed2bed(self):
        input = self.input[0]
        self.totalMappedReads = self.internalRun(bed.bed(input).getHitNum, [], self.runFlag, 'get hit number')
        outputs = []
        os.system("mkdir -p " + os.path.join(os.path.dirname(input), "chromHMM"))
        for chromHMMfile in self.chromHMMfiles:
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
            # '-l', 1000,
            '-l', 1000,
            '-m', multiplicationFactor
        ]
        self.execM(codeList)
        return self

    def appendNull_bed2bed(self):
        nullValue = "NA"
        maxEntry = self.maxChromHMMhitNumber
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
        self.saveOutput([os.path.join(self.outputDir, self.treatment + "_mergedChromHMM.txt")])
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

    def countNucleosomeDamages_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.nucleosomeBed,
            '-b', self.input,
            '-c',
            '-f', float(6)/147,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countChmm_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.fifteenChromatinStates,
            '-b', self.input,
            '-c',
            '-F', 0.49,
            '|',
            'cut',
            '-f', '1-4,9-10',
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countTSS_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.TSSbedFile,
            '-b', self.input,
            '-c',
            '-F', 0.5,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def countTES_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.TESbedFile,
            '-b', self.input,
            '-c',
            '-F', 0.5,
            '>', self.output 
        ]
        self.execM(codeList)
        return self

    def normalizeChmm_bed2bed(self):
        input = self.input[0]
        self.totalMappedReads = self.internalRun(bed.bed(input).getHitNum, [], self.runFlag and (not self.totalMappedReads), 'get hit number')
        multiplicationFactor = float(1000000)/self.totalMappedReads if self.totalMappedReads else 1 # Normalizes counts by total mapped read number
        codeList = [
            'bedCount2normalizedCount.py',
            '-i', self.input,
            '-c', 5,
            '-o', self.output,
            # '-l', 1000,
            '-l', 200,
            '-m', multiplicationFactor
        ]
        self.execM(codeList)
        return self


    def removeSitesWoDamage_bed2bed(self):
        def keepLine(line):
            keepFlag = line.strip().split('\t')[-1] != '0'
            return keepFlag
        self.internalRun(generalUtils.lineBasedFiltering, [self.input[0], self.output[0], keepLine], self.runFlag, 'remove sites without damage')
        self.internalRun(generalUtils.lineBasedFiltering, [self.input[1], self.output[1], keepLine], self.runFlag, 'remove sites without damage')
        return self

    def getClosestDamage_bed2bed(self):
        codeList = [
            'bedtools',
            'closest',
            '-a', self.nucleosomeBed,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageOnAll_bed2bed(self):
        self.input = self.input * (9-1)
        output = pipeTools.listOperation(self.addExtraWord, self.output, 'rep')
        self.saveOutput(output)
        newOutputList = []
        for i in range(1,9):
            newOutputList += pipeTools.listOperation(self.addExtraWord, self.output, str(i + 1))
        self.saveOutput(newOutputList)

        newNucleosomeList =[]
        for i in range(len(self.nucleosomeBedList)):
            for j in range(2):
                newNucleosomeList.append(self.nucleosomeBedList[i])

        codeList = [
            'bedtools',
            'closest',
            '-a', newNucleosomeList,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageTxn_bed2bed(self):
        txnList =[]
        for i in range(len(self.txnBedList)):
            for j in range(2):
                txnList.append(self.txnBedList[i])

        codeList = [
            'bedtools',
            'closest',
            '-a', txnList,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageLoop_bed2bed(self):
        output = [self.addExtraWord(self.output[0], '_Plus'), self.addExtraWord(self.output[0], '_Minus')]
        output = pipeTools.listOperation(self.addExtraWord, self.output, 'COPY')
        self.saveOutput(output)

        codeList = [
            'bedtools',
            'closest',
            '-a', self.nucleosomeBedRep1Top,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageTopNucleosomes_bed2bed(self):
        codeList = [
            'bedtools',
            'closest',
            '-a', self.nucleosomeTopBed,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageLowNucleosomes_bed2bed(self):
        codeList = [
            'bedtools',
            'closest',
            '-a', self.nucleosomeLowBed,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def getClosestDamageVoid_bed2bed(self):
        codeList = [
            'bedtools',
            'closest',
            '-a', self.void147bed,
            '-b', self.input,
            '-k', 1,
            '-D', 'a',
            '>', self.output
        ]
        self.execM(codeList)
        return self


    def adjustExactLocation_bed2bed(self):
        def adjustLocation(line):
            # NC_000001.10	10077	10223	+	2	NC_000001.10	10079	10089	+	0
            ll = line.strip().split('\t')
            nucleosomeStart = int(ll[1])
            nucleosomeEnd = int(ll[2])
            nucleosomeStrand = ll[3]
            damageStart = int(ll[6])
            damageEnd = int(ll[7])
            damageStrand = ll[8]
            distance = int(ll[9])
            damagePosition = (((damageEnd + damageStart) / 2 - nucleosomeStart))
            keepFlag = distance == 0 and damagePosition >= 0 and damagePosition <= 146
            if keepFlag:
                newLine = line.strip() + '\t' + str(damagePosition)
            else:
                newLine = None
            return newLine
        for i in range(len(self.input)):
            self.internalRun(generalUtils.lineBasedFileOperation, [self.input[i], self.output[i], adjustLocation, []], self.runFlag, 'adjustExactLocation_bed2bed')


    def adjustExactTxnLocation_bed2bed(self):
        def adjustLocation(line):
            # NC_000001.10	10077	10223	+	2	NC_000001.10	10079	10089	+	0
            #chrX	192763	193106	CTCF	193	GM12878,AG04449,AoAF,HBMEC,HCPEpiC,HEEpiC,HPAF	chrX	195157	195167	+	2052
            ll = line.strip().split('\t')
            txnStart = int(ll[1])
            txnEnd = int(ll[2])
            TF = ll[3]
            TFscore = int(ll[4])
            cells = ll[5]
            damageStart = int(ll[7])
            damageEnd = int(ll[8])
            damageStrand = ll[9]
            distance = int(ll[10])
            damagePosition = (((damageEnd + damageStart) / 2 - txnStart))
            keepFlag = distance == 0 and damagePosition >= 0 and damagePosition <= 1000
            if keepFlag:
                newLine = line.strip() + '\t' + str(damagePosition)
            else:
                newLine = None
            return newLine
        for i in range(len(self.input)):
            self.internalRun(generalUtils.lineBasedFileOperation, [self.input[i], self.output[i], adjustLocation, []], self.runFlag, 'adjustTxnExactLocation_bed2bed')


    def nucleosomeCountPlot_bed2txt(self):
        def treatmentHeader(input):
            if 'Plus' in input:
                return '"' + self.treatment + '_+"'
            elif 'Minus' in input:
                return '"' + self.treatment + '_-"'
            else:
                raise ValueError('No "Plus" or "Minus" in Input')
        def getDamageCount(input, runFlag):
            if runFlag == False:
                return 1
            else:
                bedObject = bed.bed(input)
                return float(1000000) / bedObject.getHitNum()
        # multiplicationFactor = float(1000000) / int(self.sampleMinPyrCount)
        codeList = [
            'rm', '-f', self.output,
            '&&',
            'cat', self.input,
            '|',
            'cut', '-f', 11,
            # '>', self.output
            '|',
            './plotNucleosomeDamage.R',
            self.output,
            pipeTools.listOperation(getDamageCount, self.input, self.runFlag and self.runMode),
            pipeTools.listOperation(treatmentHeader, self.input),
            73
        ]
        self.execM(codeList)
        return self 
    def TxnCountPlot_bed2txt(self):
        def treatmentHeader(input):
            if 'Plus' in input:
                return '"' + self.treatment + '_+"'
            elif 'Minus' in input:
                return '"' + self.treatment + '_-"'
            else:
                raise ValueError('No "Plus" or "Minus" in Input')
        def getDamageCount(input, runFlag):
            if runFlag == False:
                return 1
            else:
                bedObject = bed.bed(input)
                return float(1000000) / bedObject.getHitNum()
        # multiplicationFactor = float(1000000) / int(self.sampleMinPyrCount)
        codeList = [
            'rm', '-f', self.output,
            '&&',
            'cat', self.input,
            '|',
            'cut', '-f', 12,
            # '>', self.output
            '|',
            './plotNucleosomeDamage.R',
            self.output,
            pipeTools.listOperation(getDamageCount, self.input, self.runFlag and self.runMode),
            pipeTools.listOperation(treatmentHeader, self.input),
            500
        ]
        self.execM(codeList)
        return self 
    def separateStrands_bed2bed(self):

        strand = ['+', '-']
        output = [
            self.addExtraWord(self.output[0], '_Plus'), 
            self.addExtraWord(self.output[0], '_Minus')
        ]
        self.saveOutput(output)
        codeList = [
            'grep',
            strand,
            self.input[0],
            '>',
            self.output
        ]
        self.execM(codeList)
        return self

    def intersectWithTFBS_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.txnBedList[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-1000,1000]
        self.execM(codeList)
        return self
    def intersectWithDNase_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.dnaseBedList[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-1000,1000]
        self.execM(codeList)
        return self

    def intersectWithNucleosome_bed2txt(self):
        self.addWordToOutput('_1K_flanking')
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.nucleosomeBedList[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-573,573]
        self.execM(codeList)
        return self

    def parseIntersectResults_txt2txt(self):
        def parseIntersectOutput(line):
            # chr1	236749	238749	RUNX3	258	GM12878	chr1	238050	238060	-
            ll = line.split('\t')
            chromosomeA = ll[0]
            startA = int(ll[1])
            endA = int(ll[2])
            midA = (startA + endA) / 2
            startB = int(ll[7])
            endB = int(ll[8])
            midB = (startB + endB) / 2            
            damgePosition = midB - midA
            return str(damgePosition)
        for i in range(len(self.input)):
            self.internalRun(generalUtils.lineBasedFileOperation, [self.input[i], self.output[i], parseIntersectOutput, []], self.runFlag, 'parseIntersectResults_txt2txt')

    def intersectWithCTCF_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.CTCF[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-26,26]
        self.execM(codeList)
        return self

    def intersectWithBoundTFBS_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.allBoundTFBSlist[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [0,22]
        self.execM(codeList)
        return self

    def intersectWithBoundTFBS_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.allBoundTFBSlist[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-11,11]
        self.execM(codeList)
        return self

    def intersectWithBoundTFBSUpstream_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.allBoundTFBSlist[1],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-20,20]
        self.execM(codeList)
        return self

    def intersectWithBoundTFBSDownstream_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.allBoundTFBSlist[2],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-20,20]
        self.execM(codeList)
        return self

    def intersectWithSTAT3_bed2txt(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-wb',
            '-F', 0.5,
            '-a', self.STAT3[0],
            '-b', self.input,
            '>', self.output
        ]
        self.frame = [-1000,1000]
        self.execM(codeList)
        return self

    def binnedCountsToPositions_bed2txt(self):
        codeList = [
            'bedBinned2totalCounts.py',
            '-i', self.input,
            '-o', self.output,
            '-n', 200,
            '-reverseStrand', '"-"'
        ]
        self.execM(codeList)
        return self

    def intersectToPositions_txt2txt(self):
        codeList = [
            'cat', self.input,
            '|',
            'bedIntersect2parsedPosition.py',
            '|',
            'list2countTable.py',
            '|',
            'countTable2filledGaps.py',
            '-s', self.frame[0],
            '-e', self.frame[1],
            '|',
            'cut -f2'
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def listOfNumbersToCounts_txt2txt(self):
        codeList = [
            'list2countTable.py',
            '-i', self.input,
            '-o', self.output
        ]
        self.execM(codeList)
        return self

###########################################################
#  Pipeline
###########################################################
p = myPipe(input)
(p
    .run(p.cutadapt_fastq2fastq, False)

    .branch(False)
        .changeDefaultValues("hg19")
        .run(p.bowtie_fastq2sam, False)
        .run(p.convertToBam_sam2bam, False)
        .run(p.convertToBed_bam2bedpe, False)
        .run(p.uniqueSort_bedpe2bedpe, False)
        .run(p.convertBedpeToSingleFrame_bedpe2bed, False)
        .run(p.slopBed_bed2bed, False)
        .run(p.convertToFixedRange_bed2bed, False)
        .run(p.sortBed_bed2bed, False)
        .run(p.convertBedToFasta_bed2fa, False)
            
        .branch(False) # Plot nucleotide abundance
            .run(p.getNucleotideAbundanceTable_fa2csv, False)
            .run(p.plotNucleotideAbundance_csv2pdf, False)
        .stop()

        .branch(False) # Plot dinucleotide abundance
            .run(p.getDimerAbundanceTable_fa2csv, False)
            .run(p.plotDinucleotideAbundance_csv2pdf, False)
        .stop()

        .run(p.getPyrimidineDimers_fa2bed, False)
        .run(p.sortBed_bed2bed, False)
        
        # .branch(False) # ChromHMM analysis for NHF1 cells
        #     .run(p.coverageChromHMM_bed2bed, True)
        #     .run(p.normalizeCoverageChromHMM_bed2bed, True)
        #     .run(p.appendNull_bed2bed, True)
        #     .run(p.getCountColumn_bed2txt, True)
        #     .run(p.mergeChromHMM_txt2txt, True)
        #     .run(p.plotChromHMM_txt2pdf, True)
        # .stop()

        # .branch(False and p.cell == "NHF1") # CHMM
        #     .run(p.countChmm_bed2bed, False)
        #     .run(p.normalizeChmm_bed2bed, False)
        # .stop()

        .branch(False) # CHMM
            .run(p.countChmm_bed2bed, False)
            # .run(p.normalizeChmm_bed2bed, True)
        .stop()


        # .branch(False) # DNase analysis
            # .run(p.dnaseClosest_bed2bed, False)
            # .run(p.retrieveOnlyTheClosest_bed2bed, False)
            # .run(p.getDistanceFromClosest_bed2txt, False)
        # .stop()

        .branch(True)
            .run(p.separateStrands_bed2bed, False)

            .branch(False) # TSS
                .run(p.countTSS_bed2bed, False)
                .run(p.binnedCountsToPositions_bed2txt, False)
            .stop()

            .branch(False) # TES
                .run(p.countTES_bed2bed, False)
                .run(p.binnedCountsToPositions_bed2txt, False)
            .stop()


            .branch(False) # TFBS analysis
                .run(p.intersectWithTFBS_bed2txt, False)
                .run(p.intersectToPositions_txt2txt, False)
            .stop()

            .branch(False) # DNase analysis
                .run(p.intersectWithDNase_bed2txt, True)
                .run(p.intersectToPositions_txt2txt, True)
            .stop()

            # .branch(False and p.cell == "GM12878")
            #     .run(p.intersectWithNucleosome_bed2txt, False)
            #     .run(p.intersectToPositions_txt2txt, False)
            # .stop()

            # .branch(False and p.cell == "GM12878")
            #     .run(p.intersectWithCTCF_bed2txt, True)
            #     .run(p.intersectToPositions_txt2txt, True)
            # .stop()

            # .branch(False and p.cell == "GM12878")
            #     .run(p.intersectWithBoundTFBS_bed2txt, False)
            #     .run(p.intersectToPositions_txt2txt, False)
            # .stop()
            
            # .branch(False and p.cell == "GM12878")
            #     .run(p.intersectWithBoundTFBSUpstream_bed2txt, False)
            #     .run(p.intersectToPositions_txt2txt, False)
            # .stop()

            # .branch(False and p.cell == "GM12878")
            #     .run(p.intersectWithBoundTFBSDownstream_bed2txt, False)
            #     .run(p.intersectToPositions_txt2txt, False)
            # .stop()

            # .branch(True and p.cell == "GM12878")
            #     .run(p.intersectWithSTAT3_bed2txt, True)
            #     .run(p.intersectToPositions_txt2txt, True)
            # .stop()
        .stop()
    .stop()

    .branch(True and p.cell == "GM12878") # GM12878 Nucleosome Analysis
        .changeDefaultValues("hg19_nucleosome")
        .run(p.bowtie_fastq2sam, False)
        .run(p.convertToBam_sam2bam, False)
        .run(p.convertToBed_bam2bedpe, False)
        .run(p.uniqueSort_bedpe2bedpe, False)
        .run(p.convertBedpeToSingleFrame_bedpe2bed, False)
        .run(p.slopBed_bed2bed, False)
        .run(p.convertToFixedRange_bed2bed, False)
        .run(p.sortBed_bed2bed, False)
        .run(p.convertBedToFasta_bed2fa, False)
        .run(p.getPyrimidineDimers_fa2bed, False)
        .run(p.sortBed_bed2bed, False)

        .branch(True and p.cell == "GM12878")
            .run(p.separateStrands_bed2bed, False)
                
            .branch(True and p.cell == "GM12878")
                .run(p.intersectWithNucleosome_bed2txt, False)
                .run(p.intersectToPositions_txt2txt, True)
            .stop()


            # .branch(True and p.cell == "GM12878")
            #     .run(p.getClosestDamageOnAll_bed2bed, True)
            #     .run(p.adjustExactLocation_bed2bed, True)
            #     .run(p.nucleosomeCountPlot_bed2txt, True)
            # .stop()

            # .branch(False and p.cell == "GM12878")
            #     .run(p.getClosestDamageVoid_bed2bed, False)
            #     .run(p.adjustExactLocation_bed2bed, False)
            #     .run(p.nucleosomeCountPlot_bed2txt, False)
            # .stop()
        .stop()
    .stop()


    # .run(p.convertToBigBed_bed2bb, False)
)
