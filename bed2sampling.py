#!/usr/bin/env python
import os
import sys
import argparse
import generalUtils

parser = argparse.ArgumentParser(description='random spling of hits')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')
parser.add_argument('-c', required=  True, help='count')
