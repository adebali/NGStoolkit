import os
from collections import OrderedDict

# path = "/proj/sancarlb/users/ogun/seq/mm10/mm10genes.bed"
path = "/Users/ogunadebali/Downloads/mm10genes.bed"
#name	chrom	strand	txStart	txEnd	score	name2
d = OrderedDict()
filein = open(path, 'r')
for line in filein:
    if not line.startswith('#'):
        ll = line.strip().split("\t")
        chromosome = ll[1]
        start = int(ll[3])
        end = int(ll[4])
        name = ll[6]
        score = float(ll[5])
        strand = ll[2]

        d[name] = d.get(name, {})
        d[name]['chromosome'] = chromosome
        d[name]['start'] = min(d[name].get('start', start), start)
        d[name]['end'] = min(d[name].get('end', end), end)
        d[name]['name'] = name
        d[name]['score'] = score
        d[name]['strand'] = strand

for name in d:
    obj = d[name]
    fileds = [obj['chromosome'], str(obj['start']), str(obj['end']), obj['name'], str(obj['score']), obj['strand']]
    print('\t'.join(fileds))