import unittest
import os
import slurm
import datetime
import subprocess

def assignProperty(theClass, dictionary):
    for key in dictionary.keys():
        setattr(theClass, key, dictionary[key])
    return theClass

def is_binary(filename):
    """Return true if the given filename is binary.
    @raise EnvironmentError: if the file does not exist or cannot be accessed.
    @attention: found @ http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text on 6/08/2010
    @author: Trent Mick <TrentM@ActiveState.com>
    @author: Jorge Orpinel <jorge@orpinel.com>"""
    fin = open(filename, 'rb')
    try:
        CHUNKSIZE = 1024
        while 1:
            chunk = fin.read(CHUNKSIZE)
            if '\0' in chunk: # found null byte
                return True
            if len(chunk) < CHUNKSIZE:
                break # done
    finally:
        fin.close()
    return False

def getExtension(fileName):
    return fileName.strip().split('.')[-1]


def countHits(fileName):
    extension = getExtension(fileName)
    ignore = False
    if extension == 'bed' or extension == 'bedpe' or extension == 'sam' or extension == 'csv' or extension == 'txt':
        pattern = "^"
    elif extension == 'fastq':
        pattern = "^+"
    elif extension == "fa":
        pattern = "^>"
    else:
        ignore = True

    if not ignore:
        count = subprocess.check_output('grep -c "' + pattern + '" ' + fileName, shell=True)
    else:
        count = "NA"
    return count

def countLines(fileName):
    """Return cout lines of a file if the file does exist (if not returns 0) and is not binary (if so reutrns -1)."""
    if not os.path.isfile(fileName):
        return 0
    if is_binary(fileName):
        return -1
    divider = 1
    if fileName.endswith('.fastq'):
        divider = 4
    elif fileName.endswith('.fa'):
        divider = 2
    return int(subprocess.check_output('wc -l ' + fileName, shell=True).split(" ")[0])/divider

def replaceLast(source_string, replace_what, replace_with):
    head, sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail

def in2out(input, oldExtension, newExtension):
    if input.endswith(oldExtension):
        output = replaceLast(input, oldExtension, newExtension)
    else:
        output = input + '.' + getExtension(newExtension)
    return output

def listOperation(function, theList, *args):
    newList = []
    for e in theList:
        newList.append(function(e, *args))
    return newList

def funIn2out(functionName, input, extraWord = '', abbreviationLength = 3):
    fl = functionName.split('_')
    method = fl[0]
    extensions = fl[1]
    if extensions.count('2') != 1:
        raise ValueError('2 must be used and it has to be used only once. Eg: functionName_fa2csv')
    el = extensions.split('2')
    extensionIn = el[0]
    extensionOut = el[1]
    string = ''
    start = True
    for letter in fl[0]:
        if letter.isupper() or start == True:
            wordStart = 1
            string += letter
        elif wordStart < abbreviationLength:
            string += letter
            wordStart += 1
        start = False
    string += extraWord
    return in2out(input, '.' + extensionIn, '.' + string + '.' + extensionOut)

def changeDir(fullPath, dir):
    return os.path.join(dir, os.path.basename(fullPath))

def list2gappedString(theList):
    string = ''
    for e in theList:
        string += str(e) + ' '
    return string[:-1]

def list2allStringList(theList):
    newList = []
    for e in theList:
        newList.append(str(e))
    return newList

# def run(codeList, runFlag, runMode, printFlag, jobIndex, outputCheckMode):
def run(codeList, pipelineObject):
    def TrueFalse2binary(trueOrFalse):
        if trueOrFalse:
            return 1
        else:
            return 0
    code = list2gappedString(codeList)
    # code = code.replace('"', '\\"')
    # code = code.replace('(', '\\(')
    # code = code.replace(')', '\\)')
    allStringList = list2allStringList(codeList)
    branchSpace = '  '
    if pipelineObject.runFlag:
        if pipelineObject.printFlag and (not pipelineObject.outputCheckMode):
            i = 0
            print('\033[94m' + str(pipelineObject.jobIndex) + ' =>\t' + pipelineObject.currentBranch * branchSpace + code + '\033[0m')
        if pipelineObject.runMode:
            failedHere = os.system(code)
            if failedHere:
                raise ValueError('we cannot execute the code: ' + code)
        # else:
        # 	print("gave up running the code, because the command is not given in the default 'RUN' mode.")
    else:
        if pipelineObject.printFlag and (not pipelineObject.outputCheckMode):
            print(str(pipelineObject.jobIndex) + ' =X\t' + pipelineObject.currentBranch * branchSpace + code)

    if pipelineObject.outputCheckMode and (not pipelineObject.runMode) and pipelineObject.printFlag:
        i = 0
        for inputFile in pipelineObject.input:
            i += 1
            print(str(pipelineObject.jobIndex) + '.' + str(i) + ' ' + inputFile + ' ' + str(countLines(inputFile))) 


def runWm(codeList, runFlag, runMode, printFlag, wmParams, dependencies, jobIndex):
    code = list2gappedString(codeList)
    allStringList = list2allStringList(codeList)
    if runFlag:
        if printFlag:
            print(str(jobIndex) + ' -->\t' + code)
        if runMode:
            slurmObject = slurm.Slurm(code)
            slurmObject.addJobIndex(jobIndex)
            slurmObject.assignParams(wmParams)
            slurmObject.setDependencies(dependencies)
            return slurmObject.run()
            if failedHere:
                raise ValueError('we cannot execute the code: ' + code)
        # else:
        # 	print("gave up running the code, because the command is not given in the default 'RUN' mode.")
    else:
        if printFlag:
            print(str(jobIndex) + ' X\t' + code)
    return None
    # return dependencies

def codeList2multiCodeList(codeList):
    def getParallelJobNumber(codeList):
        maxListLength = 1
        for e in codeList:
            if type(e).__name__ == 'list':
                maxListLength = max(maxListLength, len(e))
        return maxListLength

    n = getParallelJobNumber(codeList)
    parallelJobLists = []
    for i in range(n):
        parallelJobLists.append([])
        for e in codeList:
            if type(e).__name__ != 'list':
                parallelJobLists[i].append(e)
            else:
                parallelJobLists[i].append(e[i])
    return parallelJobLists

def makeCommandsFromListArrays(l):
    '''Converts a list of lists andd strings and spits out the combinations
    eg. input = ['bash', ['A','B'], ['C', 1]] 
        output = makeCommandsFromListArrays(input)
        print(output)
        ['bash A C', 'bash A 1', 'bash B C', 'bash B 1']
    '''
    def makeNonlistElementListE(l):
        ''' Makes nonlist element such as string or int a list element containing a single item'''
        newL = []
        for item in l:
            if type(item) != list:
                newL.append([str(item)])
            else:
                newL.append(item)
        return newL

    tuple2list = lambda x: list(x)
    joinList = lambda x: ' '.join(str(v) for v in x)
    commandList = list(map(joinList, list(map(tuple2list, list(itertools.product(*makeNonlistElementListE(l)))))))
    return commandList

# def execM(multiCodeList, runFlag, runMode, printFlag, jobIndex, outputCheckMode):
def execM(multiCodeList, pipelineObject):
    # pipelineObject.jobIndex -= 1
    parallelJobLists = codeList2multiCodeList(multiCodeList)
    for codeList in parallelJobLists:
        # pipelineObject.jobIndex += 1
        # run(codeList, runFlag, runMode, printFlag, jobIndex, outputCheckMode)
        run(codeList, pipelineObject)
    return len(parallelJobLists)

def execL(listArray, pipelineObject):
    parallelJobLists = makeCommandsFromListArrays(listArray)
    for codeList in parallelJobLists:
        run(codeList, pipelineObject)
    return len(parallelJobLists)

def execMwm(multiCodeList, runFlag, runMode, printFlag, wmParams, dependencies, jobIndex):
    jobIndex -= 1    
    parallelJobLists = codeList2multiCodeList(multiCodeList)
    jobIdList = []
    for codeList in parallelJobLists:
        jobIndex += 1        
        jobIdList.append(runWm(codeList, runFlag, runMode, printFlag, wmParams, dependencies, jobIndex))
    return jobIdList

if __name__ == "__main__":
    import unittest

    class unitTests(unittest.TestCase):
        def test_codeList2multiCodeList(self):
            codeList = ['code.sh', '-m', 'm1', '-i', ['in1', 'in2', 'in3'], '-o', ['out1', 'out2', 'out3']]
            parallelJobLists = codeList2multiCodeList(codeList)
            results = [
                ['code.sh', '-m', 'm1', '-i', 'in1', '-o', 'out1'],
                ['code.sh', '-m', 'm1', '-i', 'in2', '-o', 'out2'],
                ['code.sh', '-m', 'm1', '-i', 'in3', '-o', 'out3']
            ]
            self.assertEqual(results, parallelJobLists)

        def test_list2gappedString(self):
            theList = ['a', 'b', 1, 89]
            expectedResult= 'a b 1 89'
            self.assertEqual(list2gappedString(theList), expectedResult)

        def test_in2out(self):
            input = '~/dir1/dir2/somefile.otherNames.extension1'
            expectedResult= '~/dir1/dir2/somefile.otherNames.extension2'
            self.assertEqual(in2out(input, 'extension1', 'extension2'), expectedResult)

        def test_replaceLast(self):
            source_string = 'someNameReplacehereOtherNameReplaceherename2Replaceherename3'
            replace_what = 'Replacehere'
            replace_with = 'ThisIsReplaced'
            expectedResult = 'someNameReplacehereOtherNameReplaceherename2ThisIsReplacedname3'
            self.assertEqual(replaceLast(source_string, replace_what, replace_with), expectedResult)

        def test_funIn2out(self):		
            self.assertEqual(funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa'), 'input.method1.getNucAbu.csv')
            self.assertEqual(funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 1), 'input.method1.gNA.csv')
            self.assertEqual(funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 2), 'input.method1.geNuAb.csv')
            self.assertEqual(funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', '', 4), 'input.method1.getNuclAbun.csv')
            self.assertEqual(funIn2out('getNucleotideAbundance_fa2csv', 'input.method1.fa', 'Extra'), 'input.method1.getNucAbuExtra.csv')

        def test_listOperation(self):
            def basicSumFunction(e, *args):
                result = e
                for arg in args:
                    print(arg)
                    result += arg
                return result

            theList = [3,4,5,6]
            self.assertEqual(listOperation(basicSumFunction, theList, 1,2,3), [9,10,11,12])

    unittest.main()