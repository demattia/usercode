__author__ = 'demattia'

import ROOT
from utils import *

# Construct combined dataset in (x,sample) and perform simultaneous fit
# Skip the last, overflow, bin from fitting to save time. It does not appear in the final plots.
# Note that the following code assumes an extra bin to build the name of the output file.
ptBins = [26, 30, 50, 80]

from array import array
hEff = ROOT.TH2D("hEff", "hEff", len(ptBins)-1, array('d',ptBins), len(ptBins)-1, array('d',ptBins))

for ptBin1 in range(0, len(ptBins)-1):
    for ptBin2 in range(0, len(ptBins)-1):

        for line in open(buildNamePars("parameters_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBins)+".txt"):
            if line.find("efficiency") != -1:
                eff = float(line.split()[2])
                effErr = float(line.split()[4])
                print "eff["+str(ptBins[ptBin1])+"_"+str(ptBins[ptBin1+1])+", "+str(ptBins[ptBin1])+"_"+str(ptBins[ptBin1+1])+"] = ", eff, "+/-", effErr

        hEff.SetBinContent(ptBin1+1, ptBin2+1, eff)
        hEff.SetBinError(ptBin1+1, ptBin2+1, effErr)

canvas = ROOT.TCanvas("efficiency", "efficiency", 800, 800)
hEff.Draw("COLZTEXTE")
canvas.Print("EfficiencyPlot.pdf")
