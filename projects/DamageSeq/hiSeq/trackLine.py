#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Write a trackline')
parser.add_argument('-n', required= True, help='name of the track')
parser.add_argument('-url', required= True, help='url')
parser.add_argument('-type', required= True, help='url')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()
url = args.url
name = args.n
out = open(args.o, 'w')

colors = ["153,50,204", "104,34,139", "0,139,139", "0,255,255"]
if 'ell' in name and 'Plus' in name:
    i = 0
elif 'ell' in name and 'Minus' in name:
    i = 1
elif 'naked' in name and 'Plus' in name:
    i = 2
elif 'naked' in name and 'Minus' in name:
    i = 3
else:
    raise ValueError('unexpected treatment_title ' + name)



out.write('track type=' + args.type + ' name="' + name + '" visibility=2 viewLimits=0:5 color=' + colors[i] + ' windowingFunction=mean smoothingWindow=16 maxHeightPixels=40:32:11 autoScale=off description="' + name + '" bigDataUrl=' + url + '\n')
out.close()