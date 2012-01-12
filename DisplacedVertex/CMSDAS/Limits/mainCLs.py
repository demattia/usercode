from GetWorkspace import *
from GetLimitCLsMass import *
import ROOT as r
import sys,os.path

if len(sys.argv)<3:
	sys.exit('usage python mainCLs.py massValue LeptonType')

r.RooMsgService.instance().deleteStream(1)

M = float(sys.argv[1])
LeptonType = sys.argv[2]

dirName = LeptonType+"CLsMass"
fileName = dirName+"/"+str(M)

os.system('mkdir -p '+dirName)

# w = GetWorkspace(M,LeptonType)
# GetLimit(w,LeptonType,str(M))

# GetLimit(M, eff, efferr, LeptonType, fileName):
GetLimit(M, LeptonType, fileName)
