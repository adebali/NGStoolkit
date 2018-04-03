import os
import sys
import urllib2  # the lib that handles the url stuff

class attrdict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.__dict__ = self

# filein = open('encoderoadmap_lasso.114.csv')
target_url = 'http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_lasso/encoderoadmap_lasso.114.csv'

# for line in filein:
for line in urllib2.urlopen(target_url):
    ll = line.split(',')
    enhancer = attrdict()
    target = attrdict()
    enhancer.location = ll[0]
    target.info = ll[1]
    score = float(ll[2])

    enhancer.chromosome = enhancer.location.split(':')[0]
    enhancer.start = int(enhancer.location.split(':')[1].split('-')[0])
    enhancer.end = int(enhancer.location.split(':')[1].split('-')[1])

    target.list = target.info.split('$')
    target.id = target.list[0]
    target.name = target.list[1]
    target.chromosome = target.list[2]
    target.start = int(target.list[3])
    target.strand = target.list[4]


    bedList = [
        enhancer.chromosome,
        str(enhancer.start),
        str(enhancer.end),
        target.name,
        str(score),
        target.strand
    ]

    print('\t'.join(bedList))