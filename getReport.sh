#!/usr/bin/bash

N=1

python pipeline.py -n $N --mock --outputCheck | cut -d' ' -f 3 | 