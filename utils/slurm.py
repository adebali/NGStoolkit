import subprocess
import random 

slurmKeyword = '#SBATCH '
firstLine = '#!/bin/bash\n'
defaultFileName = 'temp.sh'
expectedOutputStringInitial = 'Submitted batch job '

def dictionary2slurmLines(dictionary):
    lines = '\n'
    for key in dictionary.keys():
        lines += slurmKeyword + str(key) + str(dictionary[key]) + '\n'
    return lines

class Slurm:
    def __init__(self, code):
        self.code = code
        self.script = firstLine
        self.jobIndex = 0

    def assignParams(self, dictionary):
        # memoryLine = slurmKeyword + '--mem ' + str(memory) + '000' + '\n'
        # timeLine = slurmKeyword + '-t ' + str(t) + ':00:00' + '\n'
        # cpuLine = slurmKeyword + '-n ' + str(n) + '\n'
        self.parameters = dictionary
        if '-o ' in dictionary.keys():
            dictionary['-o '] = dictionary['-o '] + '.' + str(self.jobIndex)
        if '-e ' in dictionary.keys():
            dictionary['-e '] = dictionary['-e '] + '.' + str(self.jobIndex)
        self.script += dictionary2slurmLines(dictionary)
    
    def assignDefaultParameters(self):
        d = {
            "--mem=": 32000,
            "-t ": "2-00:00:00",
            "--mail-type=": "END,FAIL",
            "-e ": "error.err",
            "-o ": "log.out"
        }
        self.assignParams(d)

    def qrun(self):
        self.assignDefaultParameters()
        self.run()

    def setDependencies(self, dependencies):
        if dependencies:
            if dependencies[0] != None:
                line = ''
                for d in dependencies:
                    line += slurmKeyword + '--dependency=afterok:' + str(d) + '\n'
                self.script += line

    def write_(self, fileName):
        out = open(fileName, 'w')
        out.write(self.script)
        out.close()

    def addJobIndex(self, n):
        self.jobIndex = n

    def execute_(self, fileName):
        output = subprocess.check_output('sbatch ' + fileName, shell=True)
        return output

    def rmFile_(self, fileName):
        subprocess.check_output('rm ' + fileName, shell=True)

    def checkOutput_(self, executeOutput):
        if not executeOutput.startswith(expectedOutputStringInitial):
            raise ValueError(executeOutput + ' does not start with the expected string ' + expectedOutputStringInitial)

    def assignJobId_(self, output):
        self.checkOutput_(output)
        self.jobId = int(output.strip().split(expectedOutputStringInitial)[1])

    def getJobId_(self):
        return self.jobId

    def getRandomFileName_(self, N=6):
        return 'sb_' + ''.join(random.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(N)) + '.sh'

    def run(self):
        fileName = self.getRandomFileName_()
        self.script += self.code
        self.write_(fileName)
        output = self.execute_(fileName)
        self.rmFile_(fileName)
        self.assignJobId_(output)
        return self.getJobId_()

    def printCode(self):
        print(self.code)

    def printScript(self):
        print(self.script + self.code)
