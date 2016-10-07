# Global modules
import os
import sys
import json
from glob import glob

# Local modules
import bedUtils
import fastqUtils
import multiUtils

# Constants
kHome = os.path.expanduser("~")
kDataDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dataDir")
kChromHMMdir = os.path.join(kHome, "ogun", "chromHMM")
kSequenceDir = os.path.join(kDataDir, "160831_UNC23_0048_000000000-ARR90")
kNormalizationParameters = os.path.join(kSequenceDir, "normalizationParameters.json")


# chromHMM key
chromHMMfiles_unsorted = glob(os.path.join(kChromHMMdir, "*_chromHMM"))
chromHMMfiles = sorted(chromHMMfiles_unsorted, key=lambda name: int(name.split('chromHMM/')[1].split('_')[0]))
chromHMMlengths = multiUtils.filesFunctionField2list(chromHMMfiles, bedUtils.getTotalRegionLength, 'totalRegionLength')
###########################################################

# Sample Read Numbers key
# readFiles = []
# readFiles += glob(os.path.join(kSequenceDir, "NCPD-A*001.fastq"))
# readFiles += glob(os.path.join(kSequenceDir, "NCPD0h*001.fastq"))
# readFiles += glob(os.path.join(kSequenceDir, "NCPD1h*001.fastq"))
# readFiles += glob(os.path.join(kSequenceDir, "NCPD8h*001.fastq"))
# readFiles += glob(os.path.join(kSequenceDir, "NCPD1d*001.fastq"))
# readFiles += glob(os.path.join(kSequenceDir, "NCPD2d*001.fastq"))
#
# readNumbers = multiUtils.filesFunctionField2list(readFiles, fastqUtils.getSequenceCount, 'readNum')
###########################################################

# Mapped read files
mappedReadFiles = []
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD-A*sorted.NR.slop-l10.bed"))
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD0h*sorted.NR.slop-l10.bed"))
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD1h*sorted.NR.slop-l10.bed"))
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD8h*sorted.NR.slop-l10.bed"))
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD1d*sorted.NR.slop-l10.bed"))
mappedReadFiles += glob(os.path.join(kSequenceDir, "NCPD2d*sorted.NR.slop-l10.bed"))
readNumbers = multiUtils.filesFunctionField2list(mappedReadFiles, bedUtils.getHitNum, 'readNum')


# Actual coverage files key
fileList = []
wildCardList = [
    os.path.join(kSequenceDir, "NCPD-A*coverage*.norm.bed"),
    os.path.join(kSequenceDir, "NCPD0h*coverage*.norm.bed"),
    os.path.join(kSequenceDir, "NCPD1h*coverage*.norm.bed"),
    os.path.join(kSequenceDir, "NCPD8h*coverage*.norm.bed"),
    os.path.join(kSequenceDir, "NCPD1d*coverage*.norm.bed"),
    os.path.join(kSequenceDir, "NCPD2d*coverage*.norm.bed")
]
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
