import os    
from pipe import pipe
import argument
import bed
import pipeTools
sys.path.append('..')
from commonPipelineMethods import commonPipeline

class pipeline(commonPipeline):
    def __init__(self, input, args = argument.args()):
        pipe.__init__(self, input, args)
    
