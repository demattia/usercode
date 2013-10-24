import ROOT
from utils import *

# Define binning
ptBinsX = [26, 30, 35, 40, 45, 50, 60, 80]
ptBinsY = [26, 30, 35, 40, 45, 50, 60, 80]
#ptBins = [26, 30, 50]


# Plot the fit results

def plotResults(ptBin1, ptBin2, canvas2, canvas3):
    canvas2.cd(find_position(ptBin1, ptBin2, ptBinsX)+1)
    inputAll = ROOT.TFile.Open(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    frame1 = inputAll.FindObjectAny(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame1.Draw()

    canvas3.cd(find_position(ptBin1, ptBin2, ptBinsX)+1)
    inputPass= ROOT.TFile.Open(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    frame2 = inputPass.FindObjectAny(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame2.Draw()


from array import array
hEff = ROOT.TH2D("hEff", "hEff", len(ptBinsX)-1, array('d',ptBinsX), len(ptBinsY)-1, array('d',ptBinsY))

# Canvases for fit results
canvas2 = ROOT.TCanvas("RooFitCanvas", "RooFitCanvas", 800, 800)
canvas2.Divide(len(ptBinsX),len(ptBinsY))
canvas3 = ROOT.TCanvas("RooFitCanvasPass", "RooFitCanvasPass", 800, 800)
canvas3.Divide(len(ptBinsX),len(ptBinsY))

# Construct combined dataset in (x,sample) and perform simultaneous fit
# Skip the last, overflow, bin from fitting to save time. It does not appear in the final plots.
# Note that the following code assumes an extra bin to build the name of the output file.
for ptBin1 in range(0, len(ptBinsX)-1):
    for ptBin2 in range(0, len(ptBinsY)-1):
        plotResults(ptBin1, ptBin2, canvas2, canvas3)
        for line in open(buildNamePars("parameters_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".txt"):
            if line.find("efficiency") != -1: 
                eff = float(line.split()[2])
                effErr = float(line.split()[4])
                print "eff["+str(ptBinsX[ptBin1])+"_"+str(ptBinsX[ptBin1+1])+", "+str(ptBinsY[ptBin1])+"_"+str(ptBinsY[ptBin1+1])+"] = ", eff, "+/-", effErr
        hEff.SetBinContent(ptBin1+1, ptBin2+1, eff)
        hEff.SetBinError(ptBin1+1, ptBin2+1, effErr)

canvas2.Print("RooFitCanvas_updated.pdf")
canvas3.Print("RooFitCanvasPass_updated.pdf")

canvas4 = ROOT.TCanvas("efficiency", "efficiency", 800, 800)
hEff.Draw("COLZTEXTE")
canvas4.Print("Efficiency_updated.pdf")
