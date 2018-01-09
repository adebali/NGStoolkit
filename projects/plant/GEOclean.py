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

extensions = ['*.md5.txt', '*.GEO.txt']

for extension in extensions:
    os.system('rm ' + os.path.join(geoPath, extension))