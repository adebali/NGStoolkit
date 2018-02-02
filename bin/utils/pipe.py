import pipeTools 
import sys
import os
import argument
import itertools
kMaximumIndentation = 10


class void():
    def run(self, *args):
        return self
    def branch(self, *args):
        return self
    def stop(self, *args):
        return self
    def changeDefaultValues(self, *args):
        return self

class pipe(object):

    def __init__(self, input, args = argument.args()):
        self.input = input
        self.args = args
        self.latestInput = input
        self.inputLevels = kMaximumIndentation * [None]
        self.outputLevels = kMaximumIndentation * [None]
        self.inputLevels[0] = input
        self.currentBranch = 0
        self.justBranched = False
        self.justStopped = False
        self.branchRunFlag = True
        self.branchRunFlagLevels = kMaximumIndentation * [True]
        self.branchDependencyLevels = kMaximumIndentation * [None]
        self.i = 0
        self.jobIndex = 0
        self.subJobIndex = 0
        self.runMode = True
        self.printFlag = True
        self.wmParams = {}
        self.defaultWmParams = {}
        self.outputCheckMode = False
        self.outputCheckModeCount = False
        if '--mock' in sys.argv or self.args.get("mock"):
            self.runMode = False
        if '--noPrint' in sys.argv or self.args.get("noPrint"):
            self.printFlag = False
        if '--outputCheck' in sys.argv or self.args.get("outputCheck"):
            self.runMode = False
            self.outputCheckMode = True
            if '--count' in sys.argv or args.get("count"):
                self.outputCheckModeCount = True

    def list2attributes(self, l):
        newList = []
        for e in l:
            newList.append(self.__getattribute__(e))
        return newList

    def run(self, function, runFlag=True, *nargs):
        self.prepare_(function.__name__)
        if runFlag == "Skip":
            return self
        if runFlag == True and self.branchRunFlag == True:
            runFlag = True
        else:
            runFlag = False
        self.runFlag = runFlag
        # print(function.__name__)
        if self.justBranched == True:
            branchPlotString = self.currentBranch * '   ' + '-|* '
        else:
            branchPlotString = self.currentBranch * '   ' + ' |* '
        if self.justStopped == True:
            stopPlotString = '\n\n'
        else:
            stopPlotString = ''
        self.i += 1
        function(*nargs)
        self.finalize_()
        return self

    def cat(self, function, runFlag=True, *nargs):
        if runFlag == True and self.branchRunFlag == True and self.args.get("subprogram") == "cat":
            function(*nargs)
        return self
    
    def runSingle(self, function, runFlag=True):
        self.output = self.input
        if runFlag == True and self.branchRunFlag == True:
            runFlag = True
        else:
            runFlag = False
        self.runFlag = runFlag
        if self.justBranched == True:
            branchPlotString = self.currentBranch * '   ' + '-|* '
        else:
            branchPlotString = self.currentBranch * '   ' + ' |* '
        if self.justStopped == True:
            stopPlotString = '\n\n'
        else:
            stopPlotString = ''
        self.i += 1
        function()
        return self

    def branch(self, branchRunFlag=True):
        previousBranchFlag = self.branchRunFlagLevels[self.currentBranch]
        self.currentBranch += 1
        if previousBranchFlag == True:
            self.branchRunFlag = branchRunFlag
        else:
            self.branchRunFlag = previousBranchFlag
        self.branchRunFlagLevels[self.currentBranch] = self.branchRunFlag
        self.justBranched = True
        return self

    def stop(self):
        self.branchRunFlagLevels[self.currentBranch] = True
        self.currentBranch -= 1
        self.branchRunFlag = self.branchRunFlagLevels[self.currentBranch]
        self.justStopped = True
        return self

    def saveInput(self, nextInput):
        self.inputLevels[self.currentBranch] = nextInput
        self.latestInput = nextInput
        self.input = nextInput

    def saveOutput(self, output):
        self.outputLevels[self.currentBranch] = output
        self.latestOutput = output
        self.output = output
    
    def addExtraWord(self, fileName, word):
        fl = fileName.split('.')
        return '.'.join(fl[:-1]) + word + '.' + fl[-1]

    def addWordToOutput(self, word):
        newList = []
        for o in self.outputLevels[self.currentBranch]:
            newList.append(self.addExtraWord(o, word))
        self.outputLevels[self.currentBranch] = newList
        self.output = self.outputLevels[self.currentBranch]

    def addWordsToOutput(self, words):
        newList = []
        for word in words:
            for output in self.output:
                newList.append(self.addExtraWord(output,str(word)))
        self.outputLevels[self.currentBranch] = newList
        self.output = self.outputLevels[self.currentBranch]

    def getInput_(self):
        if self.justBranched:
            input = self.inputLevels[self.currentBranch - 1]
        else:
            input = self.inputLevels[self.currentBranch]
        if input == None:
            print("Problem here: " + str(self.inputLevels) )
        return input

    def getOutput_(self):
        return self.outputLevels[self.currentBranch]

    def prepare_(self, functionName):
        def input2output(input):
            return pipeTools.funIn2out(functionName, input, '', 2)
        input = self.getInput_()
        self.latestInput = input
        self.input = input
        output = self.listOp_(input2output, input)
        self.latestOutput = output
        self.output = output
        # print(self.inputLevels)
        # print(output)
        self.saveOutput(output)
        # print(self.inputLevels)
        
        self.justBranched = False

    def finalize_(self):
        nextInput = self.outputLevels[self.currentBranch]
        self.saveInput(nextInput)
        self.justStopped = False
        # self.inputLevels[self.currentBranch] = nextInput


    def listOp_(self, operation, var):
        if type(var).__name__ != 'list':
            newVar = operation(var)
            return newVar
        else:
            newList = []
            for e in var:
                newList.append(operation(e))
            return newList

    def execM(self, codeList):
        # jobNumber = pipeTools.execM(codeList, self.runFlag, self.runMode, self.printFlag, self.jobIndex, self.outputCheckMode)
        jobNumber = pipeTools.execM(codeList, self)
        self.jobIndex += jobNumber

    def execMwm(self, codeList):
        dependencyIndex = self.currentBranch
        jobIds = pipeTools.execMwm(codeList, self)
        self.jobIndex += len(jobIds)        
        self.branchDependencyLevels[self.currentBranch] = jobIds
        self.branchDependencyLevels[self.currentBranch + 1] = jobIds
        # self.wmParams = self.defaultWmParams
        
    def getDependencies(self):
        return self.branchDependencyLevels[self.currentBranch]

    def internalRun(self, function, arguments, runFlag=True, operationName=False, force=False):
        branchSpace = '  '
        if operationName:
            functionName = operationName
        else:
            functionName = function.__name__
        
        if runFlag:
            if self.printFlag:
                print('\033[94m' + str(self.jobIndex) + ' =>\t' + self.currentBranch * branchSpace + functionName + '\033[0m')
            if self.runMode or force:
                return function(*arguments)
        else:
            if self.printFlag:
                print(str(self.jobIndex) + ' =X\t' + self.currentBranch * branchSpace + functionName)
            return 1

    def mutateWmParams(self, dictionary):
        for key in dictionary.keys():
            self.wmParams[key] = dictionary[key]

    def exit(self):
        voidClass = void()
        return voidClass

    # def forEach(self, items, command):
    #     for item in items:
    #         self.currentKeyItem = item

    #         return self