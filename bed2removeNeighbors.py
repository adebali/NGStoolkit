#!/usr/bin/env python

import argparse
import bed

PARSER = argparse.ArgumentParser(description='removes neighboring intervals')
PARSER.add_argument('-i', required=True, help='input')
PARSER.add_argument('-d', required=True, help='neighbor distance')

ARGS = PARSER.parse_args()
bedFile = ARGS.i
distance = int(ARGS.d)

bed.bed(bedFile).removeNeighbors(distance)
