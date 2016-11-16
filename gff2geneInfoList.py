#!/usr/bin/env python
import generalUtils
import argparse
import gff

parser = argparse.ArgumentParser(description='prints a meta field from gff by order')
parser.add_argument('-i', required=True, help='<Required> input')
parser.add_argument('-o', required=True, help='<Required> output')
parser.add_argument('-f', required=True, help='<Required> field of interest')
args = parser.parse_args()

generalUtils.lineBasedFileOperation(args.i, args.o, gff.getGeneInformationFromGFFline, [args.f])