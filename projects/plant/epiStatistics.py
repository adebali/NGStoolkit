import bed
from glob import glob

files = glob('/nas/longleaf/home/adebali/ogun/seq/TAIR10/Liu2017/epigeneticMarks/*bed')

for file in sorted(files):
    bedObject = bed.bed(file)
    print(file.split('/')[-1] + '\t' + str(bedObject.getHitNum()) + '\t' + str(bedObject.getAverageLength()))
