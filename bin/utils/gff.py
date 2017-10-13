
def metaLine2metaDict(metaLine):
	metaDict = {}
	fields = metaLine.split(';')
	for field in fields:
		fl = field.split('=')
		subfieldName = fl[0]
		fieldInfo = fl[1]
		metaDict[subfieldName] = fieldInfo
	return metaDict

def getGeneInformationFromGFFline(line, field):
	result = False
	if not line.startswith('#'):
		ll = line.strip().split('\t')
		if len(ll) > 2 and ll[2] == 'gene':
			start = ll[3]
			end = ll[4]
			strand = ll[6]
			metaLine = ll[8]
			metaDict= metaLine2metaDict(metaLine)
			result = metaDict[field]
	return result