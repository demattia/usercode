from GetWorkspaceSigma import *
from GetLimitBayesian import *
import ROOT as r
import sys,os.path

if len(sys.argv)<5:
	sys.exit('usage python mainSigma.py massValue LeptonType eff effErr [fileName]')

r.RooMsgService.instance().deleteStream(1)

M = float(sys.argv[1])
LeptonType = sys.argv[2]
eff = float(sys.argv[3])
efferr = float(sys.argv[4])
if len(sys.argv)>5:
	fileName = LeptonType+"Bayesian/"+sys.argv[5]
else:
	fileName = str(M)


# w = GetWorkspaceSigma(M,LeptonType,eff,efferr)
# GetLimit(w,LeptonType,fileName)

GetLimit(M, eff, efferr, LeptonType, fileName)
