import os
import sys
from glob import glob
import re

# files = sorted(glob("dataDir/*50win.nGene.nRead.bed")) + sorted(glob("../rnaSeq/dataDir/RNA_*50win.nGene.nRead.bed"))
#files = sorted(glob("dataDir/*2900K.bowtie_ecoli.sorted.bed.TT.*.geneCov.bed")) + sorted(glob("../rnaSeq/dataDir/RNA_*genesCov.nGene.nRead.bed"))
files = sorted(glob("dataDir/*2900K.bowtie_ecoli.sorted.bed.TT.*.transcript.Cov.nRead.bed"))
output = 'dataDir/transcriptsMergedCounts.csv'
tempFile = 'dataDir/temp.csv'
headerFile = 'dataDir/headers.csv'
totalColumn = 5
unusedCol = 4
start = True
metaDataFiles = ["../rnaSeq/dataDir/FPKM.merged.noOpTr.sorted.bed"]

def run(code):
    print(code)
    if '--run' in sys.argv:
        if '--ask' in sys.argv:
            runFlag = raw_input("Shall I run this? (y/n)")
            if runFlag == 'y':
                os.system(code)
        else:
            os.system(code)

for file in files:
    basename = os.path.basename(file)
    sample = re.sub(r"[^A-Za-z]+", '', basename.split('.')[0].split('_')[0])
    if '.Minus.' in file:
        ext = '_neg'
    elif '.Plus.' in file:
        ext = '_pos'
    else:
        # ext = ''
        sys.exit('Something wrong: no .Minus. or .Plus. pattern')
    sample += ext
    if start:
        run('cat ' + file + ' >' + output)
        header = 'chromosome\tstart\tend\tstrand\t' + sample + '\t'
        start = False 
        writeMode = '>'
    else:
        run('cp ' + output + ' ' + tempFile)
        run('paste ' + tempFile + ' ' + file + ' | cut -f 1-' + str(totalColumn) + ',' + str(totalColumn + unusedCol + 1) + ' >' + output)
        run('rm ' + tempFile)
        totalColumn += 1
        header = sample + '\t'
        writeMode = '>>'
    run('echo -n \'' + header + '\' ' + writeMode + headerFile)

for file in metaDataFiles:
    basename = os.path.basename(file)
    sample = re.sub(r"[^A-Za-z]+", '', basename.split('.')[0].split('_')[0])
    run('cp ' + output + ' ' + tempFile)
    run('paste ' + tempFile + ' ' + file + ' >' + output)
    run('rm ' + tempFile)
    totalColumn += 1
    header = sample + '\t'
    writeMode = '>>'
    run('echo -n \'' + header + '\' ' + writeMode + headerFile)

run('echo -n \'\n\' ' + writeMode + headerFile)
run('cat ' + output + ' | cut -f 1-12,17 >>' + tempFile)
run('cat ' + headerFile + ' '  + tempFile + ' >' + output)
run('rm ' + headerFile)