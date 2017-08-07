import os
import sys
from pipe import pipe
import pipeTools
import generalUtils
import bed
from glob import glob
import argparse
import argument


class myPipe(pipe):
    def __init__(self, input, args = argument.args()):
        SAMPLE_STAT_FILE = 'samples.csv'  
        OUTPUT_DIR = '0125'
        pipe.__init__(self, input, args)
        self.input = os.path.realpath(os.path.join(os.path.curdir, 'dataDir', 'raw', input))
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self.attributes = sorted(sampleDictionary.keys())
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.motifRegex = '\'.{4}(G|g)(G|g).{4}\'' # Get GG only at the positions 5 and 6.
        self.referenceNickname = 'mm10'
        self.bowtie_reference = '/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome'
        self.fasta_reference = generalUtils.file('/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa')
        self.chromosome_limits = generalUtils.file('/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai')
        input1 = self.input
        input2 = generalUtils.in2out(self.input, '.1.fastq', '.2.fastq')
        self.saveInput([input1, input2])
        self.genesBedFile = "/proj/sancarlb/users/ogun/seq/mm10/geneList.bed"

        self.mm10 = {
            "name": "",
            "bowtie": "/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome",
            "fasta": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa",
            "limits": "/proj/seq/data/MM10_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai",
            "genes": "/proj/sancarlb/users/ogun/seq/mm10/geneList_yy_6.bed",
            "chmm": {
                "liver": "NA",
                "spleen": "NA",
                "kidney": "NA",
                "testes": "NA"
            }
        }
        self.mm9 = {
            "name": "_mm9",
            "bowtie": "/proj/seq/data/MM9_UCSC/Sequence/BowtieIndex/genome",
            "fasta": "/proj/seq/data/MM9_UCSC/Sequence/WholeGenomeFasta/genome.fa",
            "limits": "/proj/seq/data/MM9_UCSC/Sequence/WholeGenomeFasta/genome.fa.fai",
            "genes": "NA",
            "chmm": {
                "liver": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/liver_cStates_HMM",
                "spleen": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/spleen_cStates_HMM",
                "kidney": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/kidney_cStates_HMM",
                "testes": "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9/testes_cStates_HMM"
            }
        }
    
    def catFiles(self, wildcard, headers, output):
        code = "echo -e '" + '\t'.join(headers) + "' >" + output + " & " + "tail --lines=+2 " + wildcard + " | grep -v '^==' | grep -v -e '^$' >>" + output
        print(code)
        os.system(code)

    def fullPath2wildcard(self, fullPath):
        basename = os.path.basename(fullPath)
        directory = os.path.dirname(fullPath)
        wildcard = ".".join(["*"] + basename.split(".")[1:])
        fullWildcard = os.path.join(directory, wildcard)
        return fullWildcard

    def mergeGeneCounts(self):
        wildcard = self.fullPath2wildcard(self.input[0]).replace("_TS", "_*S")
        headers = ['chr', 'start', 'end', 'name', 'score', 'strand', 'count'] + self.attributes + ['TSNTS']
        output = os.path.join(self.outputDir, '..', 'merged_geneCounts.txt')
        self.catFiles(wildcard, headers, output)
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

    def bowtie_fastq2sam(self, referenceDictionary):
        output = pipeTools.listOperation(pipeTools.changeDir, self.output, self.outputDir)
        self.saveOutput([self.addExtraWord(output[0], referenceDictionary["name"])])
        codeList = [
            'bowtie',
            '-t', referenceDictionary["bowtie"],
            '-q', # FASTAQ input (default)
            '--nomaqround', # Do NOT round MAC
            '--phred33-quals', # Depends on the sequencing platform
            '-S', # Output in SAM format
            '-m', 4, # Do not report the reads that are mapped on to more than 4 genomic locations
            '-X', 1000,
            '--seed', 123, # Randomization parameter in bowtie,
            '-p', 8,
            '-1', self.input[0],
            '-2', self.input[1],
            self.output
        ]
        self.referenceDictionary = referenceDictionary
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

    def getDamageSites_fa2bed(self):
        codeList = [
            'fa2bedByChoosingReadMotifs.py',
            '-i', self.input,
            '-o', self.output,
            '-r', self.motifRegex 
        ]
        self.execM(codeList)
        return self

    # def sampleFromBed_bed2bed(self):
    #     codeList = [
    #         'subsample',
    #         '-n', self.sampleMinPyrCount,
    #         '--seed', 123,
    #         self.input,
    #         '>', self.output
    #     ]
    #     self.execM(codeList)
    #     return self

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
            '-k3,3n',
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
        self.finalBed = self.output[0] 
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

    def geneMap_bed2bed(self):
        codeList = [
            'bedtools',
            'intersect',
            '-wa',
            '-a', self.genesBedFile,
            '-b', self.input,
            '-c',
            '-F', 0.5,
            '>', self.output
        ]
        self.execM(codeList)
        return self

    def geneStrandMap_bed2bed(self):
        newOutput = [self.addExtraWord(self.output[0], '_TS'), self.addExtraWord(self.output[0], '_NTS')]
        self.saveOutput(newOutput)
        strandParameters = ['-S', '-s'] #[different strand, same strand]
        self.input = [self.input[0], self.input[0]]
        codeList = [
            'bedtools',
            'intersect',
            '-a', self.referenceDictionary["genes"],
            '-b', self.input,
            '-wa',
            '-c',
            strandParameters,
            '-F', 0.50,
            '>', self.output
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

    def categorizeTSvsNTS_bed2txt(self):
        codeList = [
            'categorize_TSvsNTS.py',
            '-plus', self.input[0],
            '-minus', self.input[1],
            '-s', 4,
            '-o', self.output[0]
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

    def uniqueSort_bed2bed(self):
        codeList = [
            'sort',
            '-u',
            '-k1,1',
            '-k2,2n',
            '-k3,3n',
            self.input,
            '>', self.output
        ]
        self.finalBed = self.output[0] 
        self.execM(codeList)
        return self

    def addTreatmentAndStrand_txt2txt(self):
        columns = self.list2attributes(self.attributes)
        columnStringList = [' '.join(columns + ['TS']), ' '.join(columns + ['NTS'])]
        codeList = [
            'addColumns.py',
            '-i', self.input,
            '-o', self.output,
            '-c', columnStringList
        ]
        self.execM(codeList)
        return self

def getArgs():
    parser = argparse.ArgumentParser(description='XR-seq ZT Pipeline', prog="pipeline.py")
    parser.add_argument('--outputCheck', required= False, default=False, action='store_true', help='checkOutput flag')
    
    subparsers = parser.add_subparsers(help='pipeline help', dest="subprogram")

    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('-n', required= True, help='output')
    parser_run.add_argument('--mock', required= False, default=False, action='store_true', help='mock flag')
    parser_run.add_argument('--noPrint', required= False, default=False, action='store_true', help='prints no code when stated')

    parser_cat = subparsers.add_parser('cat', help='cat help')
    parser_cat.add_argument('-n', required= False, default="1", help='output')
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
###########################################################
#  Pipeline
###########################################################
p = myPipe(input, args)
(p
    .run(p.cutadapt_fastq2fastq, False)
    .run(p.bowtie_fastq2sam, False, p.mm10)
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

    .run(p.getDamageSites_fa2bed, False)
    .run(p.uniqueSort_bed2bed, False)

    .branch(True)
        .run(p.geneStrandMap_bed2bed, False)
        .run(p.normalizeCounts_txt2txt, False)
        
        .branch()
            .run(p.addTreatmentAndStrand_txt2txt, True)
            .cat(p.mergeGeneCounts, True)
        .stop()
    .stop()    

    .run(p.separateStrands_bed2bed, False)

    .branch(False)
        .run(p.geneMap_bed2bed, False)
        .run(p.categorizeTSvsNTS_bed2txt, False)
    .stop()
)

