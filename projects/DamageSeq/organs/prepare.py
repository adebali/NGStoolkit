import os
import sys
from glob import glob
import re

gunzipFiles = []
# gunzipFiles += sorted(glob('/proj/sancarlb/HTSF/161010_UNC17-D00216_0412_BCA3T0ANXX/*.fastq.gz'))
# gunzipFiles += sorted(glob('/proj/sancarlb/HTSF/161102_UNC32-K00270_0028_AHFTWJBBXX/*.fastq.gz'))
gunzipFiles += sorted(glob('/proj/sancarlb/HTSF/170223_UNC18-D00493_0399_ACAJD6ANXX/*.fastq.gz'))
outDir = '/proj/sancarlb/users/ogun/DamageSeq/organs/raw/'
fastqFiles =[]

def run(code):
    print(code)
    if '--run' in sys.argv:
        os.system(code)
    
startedFileList = []


for gunzipFile in gunzipFiles:
    fastqFile = gunzipFile[:-3]
    fastqBasename = os.path.basename(fastqFile)
    currentFastq = outDir + fastqBasename
    fastqFiles.append(currentFastq)

    fields = re.split("[_\-]", fastqBasename)
    if not fields[2]:
        fields[2] = "0"
    name = "_".join(fields[0:3])
    index = fields[3]
    # Lnumber = fields[2]
    Lnumber = fields[5]
    # R1R2 = fields[3]
    R1R2 = fields[6]

    if R1R2 == 'R1':
        output = outDir + name + '.1.fastq'
    elif R1R2 == 'R2':
        output = outDir + name + '.2.fastq'
    else:
        sys.exit('Unexpected R1R2 field')
    
    if output in startedFileList:
        writeMode = '>>'
    else:
        writeMode = '>'

    startedFileList.append(output) 

    code = 'gunzip -c ' + gunzipFile + ' ' + writeMode + output
    run(code)