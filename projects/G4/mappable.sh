#!/usr/bin/bash

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -e \
    "select chrom, size from hg19.chromInfo" \
| awk -v OFS="\t" '{print $1, 0, $2}' \
| subtractBed -a - -b wgEncodeDukeMapabilityRegionsExcludable.bed.gz \
| sort -k1,1 -k2,2n > hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed