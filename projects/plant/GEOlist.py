import os
import sys
sys.path.append('..')
from referenceGenomePath import referenceGenomePath
from glob import glob
import generalUtils

SAMPLE_STAT_FILE = 'samples.csv'    
sampleNumbers = generalUtils.table2dictionary(SAMPLE_STAT_FILE, 'no')


paths = referenceGenomePath()
reference = paths.get('TAIR10')
geoPath = reference['geo']


sampleNoList = range(3,18+1)

for i in sampleNoList:
    sampleName = sampleNumbers[str(i)][0]['sample'].split('.')[0]
    extension = 'fastq'
    filename = glob(os.path.join(geoPath, sampleName + '*' + extension))[0]    
    print(os.path.basename(filename) + '\t' + open(glob(os.path.join(geoPath, filename.replace(extension, '')) + 'GEO*')[0]).readline().strip().split(' ')[0])
    
for i in sampleNoList:
    sampleName = sampleNumbers[str(i)][0]['sample'].split('.')[0]
    extension = 'bed'
    filename = glob(os.path.join(geoPath, sampleName + '*' + extension))[0]
    print(os.path.basename(filename) + '\t' + open(glob(os.path.join(geoPath, filename.replace(extension, '')) + 'GEO*')[0]).readline().strip().split(' ')[0])
    