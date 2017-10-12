import os

def f (string, *args):
    for i in args[0]:
        print(i)
        if len(args) > 1:
            for j in args[1]:
                string += j
                return f(string, args[1:])
        elif args == []:
            return string 

# a = f('', ['A','B'], ['C', 'D'], ['E', 'F', 'G', 'H'])
# print(a)

import itertools
d = [{'A':'_A','B':'_B'}, {'C':'_C', 'D':'_D'}, {'E':'_E', 'F':'_F', 'G':'_G', 'H':'_H'}]


def makeCommandsFromListArrays(l, separator= ' '):
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
    joinList = lambda x: separator.join(str(v) for v in x)
    commandList = list(map(joinList, list(map(tuple2list, list(itertools.product(*makeNonlistElementListE(l)))))))
    return commandList

def ListArrayAnnotationToCommands(l, annotation, outputNo, processWordAndOutputExtension):
    '''annotation = ['input', ['_C', '_1'], ['_E', '_F', '_G', '_H'], 'output']
    '''
    newList = []
    i = 0
    for item in annotation:
        i += 1
        if type(item) == list:
            newList.append(item)
        elif item == 'input':
            inputIndex = i
        elif item == 'output':
            outputIndex = i

    annotationList = list(itertools.product(*newList))
    annotationStringList = makeCommandsFromListArrays(annotationList, '')
    for outputFile in l[outputIndex - 1]:
        baseFileName = '.'.join(outputFile.split('.')[:-1])


    for inputFile in l[inputIndex - 1]:
        baseFileName = '.'.join(outputFile.split('.')[:-1])
        

    renameOutput(inputFile, normalOutput)
    


a = ['bash', ['A','B'], ['C', 1], 'X1', ['E', 'F', 'G', 'H']]
print(makeCommandsFromListArrays(a))

print(makeCommandsFromListArrays(['bash', 'A','B', 'C', 1, 'X1']))
input = ['bash', ['A','B'], ['C', 1]] 
output = makeCommandsFromListArrays(input)
print(output)

['bash', ['A.txt','B.txt'], ['C', 1], 'X1', '-o', '$output', ['E', 'F', 'G', 'H']], '_process1.bed', 2, [['_C', '_1'], ['_E', '_F', '_G', '_H']]


listOfArrays, keyword, inputIndex, annotationList

inputList = listOfCommand[inputIndex - 1]

listOfCommand = ['bash', 'A', 'C', 'X1', '-o', '$output', 'E']
inputFile = listOfCommand[inputIndex - 1]
baseFileName = '.'.join(inputFile.split('.')[:-1])
outputFile = baseFileName + ''.join(listOfCommand[]


'bash A.txt C X1 -o $output E'



[
    'bash A.txt C X1 -o A_C_E_process1.bed E',
    'bash A.txt C X1 -o A_C_F_process1.bed F',
    'bash A.txt C X1 -o A_C_G_process1.bed G',
    'bash A.txt C X1 -o A_C_H_process1.bed H',
    'bash A.txt 1 X1 -o A_1_E_process1.bed E',
    'bash A.txt 1 X1 -o A_1_F_process1.bed F',
    'bash A.txt 1 X1 -o A_1_G_process1.bed G',
    'bash A.txt 1 X1 -o A_1_H_process1.bed H',
    'bash B.txt C X1 -o B_C_E_process1.bed E',
    'bash B.txt C X1 -o B_C_F_process1.bed F',
    'bash B.txt C X1 -o B_C_G_process1.bed G',
    'bash B.txt C X1 -o B_C_H_process1.bed H',
    'bash B.txt 1 X1 -o B_1_E_process1.bed E',
    'bash B.txt 1 X1 -o B_1_F_process1.bed F',
    'bash B.txt 1 X1 -o B_1_G_process1.bed G',
    'bash B.txt 1 X1 -o B_1_H_process1.bed H'
]



