from pipeline import *

class seqContextPipeline(pipeline):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
        OUTPUT_DIR = '0707'
        SAMPLE_STAT_FILE = 'seqSamples.csv'    
        # self.input = os.path.join('dataDir', OUTPUT_DIR, 'Hetzel_sorted_noEdge10K_s10K_TT.bed')
        self.input = os.path.join('dataDir', OUTPUT_DIR, 'genome_TT.bed')
        self.outputDir = os.path.realpath(os.path.join(os.path.dirname(self.input), '..', OUTPUT_DIR))
        os.system('mkdir -p ' + self.outputDir)
        sampleDictionary = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'sample')[input][0]
        self = pipeTools.assignProperty(self, sampleDictionary)
        self.attributes = sorted(sampleDictionary.keys())
        self.saveInput([self.input])
        self.defaultWmParams = {
            '--mem=': 32000,
            '-n ': 1,
            '-t ': '24:00:00',
            '--job-name=': 'XR-seq',
            '-o ': 'log_' + self.treatment + '.txt',
            '-e ': 'err_' + self.treatment + '.txt',
        }
        self.wmParams = self.defaultWmParams
        self.paths = referenceGenomePath()
        self.reference = self.paths.get('TAIR10')

def getInputFromIndex(n):
    SAMPLE_STAT_FILE = 'seqSamples.csv'
    return sampleIO(SAMPLE_STAT_FILE, n, 'no', 'sample')

args = getArgs()
inputIndex = args.get("n")
input = getInputFromIndex(inputIndex)


if __name__ == "__main__":
     p = seqContextPipeline(input, args)
     (p
        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True)
            .run(p.addValue_txt2txt, True, 'real')
        .stop()

        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True, {'random':True})
            .run(p.addValue_txt2txt, True, 'random')
        .stop()

        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True, {"slice": True, "n": 4, "sliceNum": 1, "keyword": '_Q1'})
            .run(p.addValue_txt2txt, True, 'Q1')
        .stop()

        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True, {"slice": True, "n": 4, "sliceNum": 2, "keyword": '_Q2'})
            .run(p.addValue_txt2txt, True, 'Q2')
        .stop()

        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True, {"slice": True, "n": 4, "sliceNum": 3, "keyword": '_Q3'})
            .run(p.addValue_txt2txt, True, 'Q3')
        .stop()

        .branch(False)
            .run(p.transcriptIntersect_bed2txt, True, {"slice": True, "n": 4, "sliceNum": 4, "keyword": '_Q4'})
            .run(p.addValue_txt2txt, True, 'Q4')
            .cat(p.mergeTCR, True, '_seq', 'genome_TT')
        .stop()

        .branch(True)
            .run(p.dnaseIntersect_bed2txt, True)
            .run(p.addValue_txt2txt, True, 'DNase')
        .stop()

        .branch(True)
            .run(p.dnaseIntersect_bed2txt, True, {'random': True})
            .run(p.addValue_txt2txt, True, 'DNase')
            .cat(p.mergeDNase, True, '_seq', 'genome_TT')
        .stop()
     )
        
