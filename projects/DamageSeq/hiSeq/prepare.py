import os
import sys
from glob import glob

gunzipFiles = []
# gunzipFiles += sorted(glob('/proj/sancarlb/HTSF/161010_UNC17-D00216_0412_BCA3T0ANXX/*.fastq.gz'))
gunzipFiles += sorted(glob('/proj/sancarlb/HTSF/161102_UNC32-K00270_0028_AHFTWJBBXX/*.fastq.gz'))
outDir = '/proj/sancarlb/users/ogun/DamageSeq/hiSeq/'
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

    fields = fastqBasename.split('_')
    name = fields[0]
    index = fields[1]
    Lnumber = fields[2]
    R1R2 = fields[3]

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