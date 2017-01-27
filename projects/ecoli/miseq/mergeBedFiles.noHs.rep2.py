import os
import sys
from glob import glob
import re

# files = sorted(glob("dataDir/*50win.nGene.nRead.bed")) + sorted(glob("../rnaSeq/dataDir/RNA_*50win.nGene.nRead.bed"))
files = sorted([fn for fn in glob("dataDir/STL*cutadapt13.2900K.bowtie_ecoli.sorted.bed.TT.noHS.*.gene.Cov.bed") if not os.path.basename(fn).startswith('STLWT4')]) + sorted(glob("../rnaSeq/dataDir/RNA_*genesCov.nGene.nRead.bed"))
output = 'dataDir/genesMergedCounts_noHs_rep2.csv'
tempFile = 'dataDir/temp.csv'
headerFile = 'dataDir/headers.csv'
totalColumn = 5
unusedCol = 4
start = True
metaDataFiles = ["/nas02/home/a/d/adebali/ncbi/ecoli/NC_000913.2/NC_000913.2.genes.txt"]

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
    sample = re.sub(r"[^A-Za-z]+", '', basename.split('.')[0].split('_')[0]).replace("STL",'')
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
run('cp ' + output + ' ' + tempFile)
run('cat ' + headerFile + ' '  + tempFile + ' >' + output)
run('rm ' + headerFile)