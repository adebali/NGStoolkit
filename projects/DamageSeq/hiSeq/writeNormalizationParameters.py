# Global modules
import os
import sys
import json
from glob import glob

# Local modules
import bed
import multiUtils

# Constants
kHome = os.path.expanduser("~")
kDataDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dataDir')
kChromHMMdir = os.path.join(kHome, "ogun", "chromHMM")
kSequenceDir = kDataDir
kNormalizationParameters = os.path.join(kSequenceDir, "normalizationParameters.json")


# chromHMM key
chromHMMfiles_unsorted = glob(os.path.join(kChromHMMdir, "*_chromHMM"))
chromHMMfiles = sorted(chromHMMfiles_unsorted, key=lambda name: int(name.split('chromHMM/')[1].split('_')[0]))
chromHMMlengths = multiUtils.filesClassFunctionField2list(chromHMMfiles, bed.bed, 'getTotalRegionLength', 'totalRegionLength')
###########################################################

bedFileList = [
	"N6-0HB19N.cutadapt.1.bowtie_hg19.uSorted.sf.slopB3.frL10.bed",
	"N6-20MB9N.cutadapt.1.bowtie_hg19.uSorted.sf.slopB3.frL10.bed",
	"N6-1HB20N.cutadapt.1.bowtie_hg19.uSorted.sf.slopB3.frL10.bed",
	"N6-2HB12N.cutadapt.1.bowtie_hg19.uSorted.sf.slopB3.frL10.bed",
	"N6-4HB21N.cutadapt.1.bowtie_hg19.uSorted.sf.slopB3.frL10.bed"
]
# Mapped read files
mappedReadFiles = []
for file in bedFileList:
	mappedReadFiles.append(os.path.join(kSequenceDir, file))

print(mappedReadFiles)
readNumbers = multiUtils.filesClassFunctionField2list(mappedReadFiles, bed.bed, 'getHitNum', 'readNum')
print(readNumbers)
# Actual coverage files key
fileList = []
wildCardList = []
for file in bedFileList:
    wildCardList.append(os.path.join(kSequenceDir, file))

for wildcard in wildCardList:
    fileList += sorted(glob(wildcard), key=lambda name: int(name.split('coverage_')[1].split('_')[0]))
###########################################################

## only for local machine
newFileList = []
for file in fileList:
    newFileList.append(file.split('/')[-1])

fileList = list(newFileList)

# Write to a JSON file
myDict = {
    # "chromHMM": chromHMMlengths,
    "sampleReadNumbers": readNumbers,
    "fileList": fileList
}

with open(kNormalizationParameters, 'w') as f:
     json.dump(myDict, f, indent=4)

with open(kNormalizationParameters, 'a') as f:
    f.write('\n')
###########################################################
