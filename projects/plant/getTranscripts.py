import os
import sys
from fasta import fasta

directory = "/nas/longleaf/home/adebali/ogun/seq/TAIR9"
os.chdir(directory)
# os.system("wget -O transcripts.fa http://www.plantgdb.org/download/Download/xGDB/AtGDB/ATtranscriptTAIR9")

myFasta = fasta("transcripts.fa")
transcripts = open("transcripts.bed", "w")
minLength = 4000
transcripts_gtCutoff = open("transcripts_gt" + str(minLength)[0] + "K.bed", "w")

chromosomeDict = {
    "chr1": "1",
    "chr2": "2",
    "chr3": "3",
    "chr4": "4",
    "chr5": "5",
    "chrM": "Mt",
    "chrC": "Pt"
}

strandDict = {"FORWARD": "+", "REVERSE": "-"}

nameDict = {}

for sequence in myFasta.stream():
    header = sequence['h']
    name = header.split("|")[0].strip().replace(">","")
    versionlessName = name.split('.')[0]
    location = header.split("|")[-1].strip()
    position = location.split(" ")[0]
    strand = strandDict[location.split(" ")[1]]
    ll = position.split(":")
    chromosome = chromosomeDict[ll[0].strip()]
    interval = ll[1]
    interval_list = interval.split("-")
    start = int(interval_list[0]) - 1
    end = int(interval_list[1]) - 1
    currentDict = nameDict.get(versionlessName, {"chromosome": chromosome, "start": start, "end": end, "strand": strand})
    nameDict[versionlessName] = currentDict
    nameDict[versionlessName]["start"] = min(start, currentDict["start"])
    nameDict[versionlessName]["end"] = min(end, currentDict["end"])


for name in nameDict.keys():
    lineDict = nameDict[name]
    bedLine = "\t".join([lineDict['chromosome'], str(lineDict['start']), str(lineDict['end']), name, ".", lineDict['strand']]) + "\n"
    transcripts.write(bedLine)
    if lineDict['end'] - lineDict['start'] >= minLength:
        transcripts_gtCutoff.write(bedLine)
transcripts.close()
transcripts_gtCutoff.close()
