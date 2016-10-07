#!/bin/sh

loopThroughFiles.py -code "bsub python pipeline.py #IN" -files "dataDir/160831_UNC23_0048_000000000-ARR90/*001.fastq"
