import os
import sys
import json
import copy

class samples():
    def __init__(self, inputJson):
        self.input = inputJson
        self.sampleDict = json.load(open(self.input))
        self.qualityTest()
        self.completeSampleDict = self.completeSamples()

    def qualityTest(self):
        experiments = self.filterDictionary(self.sampleDict, 'isExperiment', True)
        experimentNoList = []
        for key in experiments.keys():
            if experiments[key].get('experimentNo', False):
                experimentNoList.append(experiments[key].get('experimentNo'))
            else:
                raise ValueError('Experiment must have a field of "experimentNo"')
        if len(experimentNoList) != len(set(experimentNoList)):
            print(experimentNoList)
            raise ValueError('Duplicated experiment no was found!')

    def key2attributes(self, key):
        def recursiveBase(d):
            if d.get('base', False):
                baseD = copy.copy(self.sampleDict[d['base']])
                d_updatedWithBase = copy.copy(recursiveBase(baseD))
                d_updatedWithBase.update(d)
                return d_updatedWithBase
            return d
        singleSampleDict = self.sampleDict[key]
        completeSampleDict = recursiveBase(singleSampleDict)
        return completeSampleDict

    def completeSamples(self):
        completeDict = {}
        for key in self.sampleDict.keys():
            sample = self.sampleDict[key]
            if sample.get('template', False) != True:
                completedSample = self.key2attributes(key)
                completeDict[key] = completedSample
        return completeDict

    def filterDictionary(self, dictionary, key, value):
        return {k : v for k,v in dictionary.iteritems() if key in v.keys() and v[key] == value}
        
    def getSamplesByKey(self, key, value):
        return self.filterDictionary(self.completeSampleDict, key, value)
        
    def getSamplesByExperimentName(self, experimentName):
        return self.getSamplesByKey('experiment', experimentName)

    def getSamplesByExperimentNo(self, experimentNo):
        return self.getSamplesByKey('experimentNo', experimentNo)

    def experimentNoAndNo2id(self, experimentNo, no):
        filteredByExperiment = self.getSamplesByExperimentNo(experimentNo)
        filteredByNo = self.filterDictionary(filteredByExperiment, 'no', no)
        if len(filteredByNo.keys()) != 1:
            raise ValueError('Expected 1 but got ' + str(len(filteredByNo.keys())) + ' samples!')
        return filteredByNo.keys()[0]

    def experiment2sampleNoList(self, experimentNo):
        noList = []
        experimentSamples = self.getSamplesByExperimentNo(experimentNo)
        for key in experimentSamples.keys():
            noList.append(experimentSamples[key]['no'])
        return sorted(noList)

from pprint import pprint

class paths():
    def __init__(self, inputJson):
        self.input = inputJson
        self.pathDict = json.load(open(self.input))

    def genome2fullPaths(self, genome):
        completeDict = {}
        for key in self.pathDict[genome].keys():
            completeDict[key] = os.path.join(self.pathDict[self.pathDict[genome][key]['base']], self.pathDict[genome][key]['location'])
        return completeDict
    
    def getGenome(self, genome):
        return self.genome2fullPaths(genome)

    def getPath(self, genome, key):
        return self.getGenome(genome)[key]


def test_paths():
    pathClass = paths('reference.json')
    assert pathClass.getPath('hg19', 'bowtie2') == '/proj/seq/data/HG19_UCSC/Sequence/Bowtie2Index/genome'

def test_samples():
    sampleClass = samples('sample.json')