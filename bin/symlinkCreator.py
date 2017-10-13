#!/usr/bin/env python
import os
import sys
from glob import glob

pathString = sys.argv[1]
fileList = glob(pathString)

for file in fileList:
	os.system('ln -s ' + file + './' + file.split('/')[-1])
