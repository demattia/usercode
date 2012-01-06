import sys
import ROOT as r
from GetWorkspace import *

if len(sys.argv)<4:
	sys.exit('usage: python GenerateFakeData.py LeptonType mean nsig nbkg')

M = float(sys.argv[2])
nsig = int(sys.argv[3])
nbkg = int(sys.argv[4])

w = GetWorkspace(M,sys.argv[1])

s = r.RooRealVar("s","s",nsig)
b = r.RooRealVar("b","b",nbkg)

sigPdf = w.pdf("signalPdf")
bkgPdf = w.pdf("bkgPdf")

mass = w.var("mass")
model = r.RooAddPdf("model","model",r.RooArgList(sigPdf,bkgPdf),r.RooArgList(s,b))

data = model.generate(r.RooArgSet(mass),nsig+nbkg)
fdata = r.TFile("data.root",'recreate')
data.Write()
fdata.Close()

frame = mass.frame()
data.plotOn(frame)
frame.Draw()
