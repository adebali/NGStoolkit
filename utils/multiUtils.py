
def filesFunctionField2list(files, func, field):
	theList = []
	for file in files:
		record = {
			"name": file,
			field: func(file)
		}
		theList.append(record)
	return theList

def filesClassFunctionField2list(files, Class, functionString, field):
	theList = []
	for file in files:
		myClass = Class(file)
		method = getattr(myClass, functionString)
		record = {
			"name": file,
			field: method()
		}
		theList.append(record)
	return theList
