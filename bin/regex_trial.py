import re


def checkMyPattern(sequence, Gposition):
    startPostion = Gposition - 3
    pattern = '.{'+ str(startPostion) + '}(?!G).(?!G).G(?!G).(?!G).*'
    prog = re.compile(pattern)
    result = prog.match(sequence)
    if result:
        return True
    else:
        return False

print(checkMyPattern('AAAAATTTTTAAAAATTTGAAA', 19))
print(checkMyPattern('AAAAATTTTTAAAAATTTTGAAA', 19))
print(checkMyPattern('AAAAATTTTTAAAAATTTTGAAA', 20))