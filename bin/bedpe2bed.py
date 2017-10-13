#!/usr/bin/env python
import bed
import sys

bed.bedpe(sys.argv[1]).mergeFragments()
