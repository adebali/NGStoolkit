#!/usr/bin/env python
import os
import sys

fastq = sys.argv[1]
os.system("perl -ne 'y/@/>/;print($_.<>)&&<>&&<>' " + fastq)