def CalcParam(mean,line):
	line = line.split()
	parName = line[0]
	poly = [eval(a) for a in line[1:]]
	parVal = 0
	for i in range(len(poly)):
		parVal += poly[i]*pow(mean,i)
	return [parName,parVal]

def GetParameters(mean,parFileName):

	parDict = {}

	parFile = open(parFileName)
	for line in parFile.readlines():
		par = CalcParam(mean,line)
		parDict[par[0]]=par[1]
	print parDict
	return parDict

