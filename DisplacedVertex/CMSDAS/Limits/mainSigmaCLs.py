from GetWorkspaceSigma import *
from GetLimitCLs import *
import ROOT as r
import sys,os.path

# if len(sys.argv)<4:
# 	sys.exit('usage python mainSigmaCLs.py massValue LeptonType eff [fileName]')
# 
# r.RooMsgService.instance().deleteStream(1)
# 
# M = float(sys.argv[1])
# LeptonType = sys.argv[2]
# eff = float(sys.argv[3])
# if len(sys.argv)>4:
# 	fileName = LeptonType+"CLs/"+sys.argv[4]
# else:
# 	fileName = str(M)

if len(sys.argv)<5:
	sys.exit('usage python mainSigma.py massValue LeptonType eff effErr [fileName]')

r.RooMsgService.instance().deleteStream(1)

M = float(sys.argv[1])
LeptonType = sys.argv[2]
eff = float(sys.argv[3])
efferr = float(sys.argv[4])
if len(sys.argv)>5:
	fileName = LeptonType+"CLs/"+sys.argv[5]
else:
	fileName = str(M)

# w = GetWorkspaceSigma(M,LeptonType,eff)
# GetLimit(w,LeptonType,fileName)

GetLimit(M,eff,efferr,LeptonType,fileName)
