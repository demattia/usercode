import ROOT
from utils import *

# Define binning
ptBinsX = [15, 20, 23, 26, 30, 35, 40, 45, 50, 60, 70]
#ptBinsX = [15, 70]
ptBinsY = [20, 10000]


# Plot the fit results

def plotResults(ptBin1, ptBin2, canvas2, canvas3):
    canvas2.cd(find_position_NoOverflow(ptBin1, ptBin2, ptBinsX)+1)
    inputAll = ROOT.TFile.Open(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    inputAll_nomatch = ROOT.TFile.Open(buildNamePars("fitAll_nomatch", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    frame1 = inputAll.FindObjectAny(buildNamePars("fitAll_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame1_nomatch = inputAll_nomatch.FindObjectAny(buildNamePars("fitAll_nomatch", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame1.Draw()
    frame1_nomatch.Draw("same")
    frame1.SetMarkerColor(1)
    frame1_nomatch.SetMarkerColor(2)
    frame1.SetLineColor(1)
    frame1_nomatch.SetLineColor(2)

    canvas3.cd(find_position_NoOverflow(ptBin1, ptBin2, ptBinsX)+1)
    inputPass= ROOT.TFile.Open(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    inputPass_nomatch= ROOT.TFile.Open(buildNamePars("fitPass_nomatch", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY)+".root")
    frame2 = inputPass.FindObjectAny(buildNamePars("fitPass_", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame2_nomatch = inputPass_nomatch.FindObjectAny(buildNamePars("fitPass_nomatch", ptBin1, ptBin1+1, ptBin2, ptBin2+1, ptBinsX, ptBinsY))
    frame2.Draw()
    frame2_nomatch.Draw("same")
    frame2.SetMarkerColor(1)
    frame2_nomatch.SetMarkerColor(2)
    frame2.SetLineColor(1)
    frame2_nomatch.SetLineColor(2)
 #   frame1.Draw("same")


from array import array
hEff = ROOT.TH2D("hEff", "hEff", len(ptBinsX)-1, array('d',ptBinsX), len(ptBinsY)-1, array('d',ptBinsY))
hEff1D = ROOT.TH1D("hEff1D", "hEff1D", len(ptBinsX)-1, array('d',ptBinsX))

# Canvases for fit results
canvas2 = ROOT.TCanvas("RooFitCanvas", "RooFitCanvas", 2000, 300)
canvas2.Divide(len(ptBinsX)-1,len(ptBinsY)-1)
canvas3 = ROOT.TCanvas("RooFitCanvasPass", "RooFitCanvasPass", 2000, 300)
canvas3.Divide(len(ptBinsX)-1,len(ptBinsY)-1)

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
                print "eff["+str(ptBinsX[ptBin1])+"_"+str(ptBinsX[ptBin1+1])+", "+str(ptBinsY[ptBin2])+"_"+str(ptBinsY[ptBin2+1])+"] = ", eff, "+/-", effErr
        hEff.SetBinContent(ptBin1+1, ptBin2+1, eff)
        hEff.SetBinError(ptBin1+1, ptBin2+1, effErr)
        hEff1D.SetBinContent(ptBin1+1, eff)
        hEff1D.SetBinError(ptBin1+1, effErr)


canvas2.Print("RooFitCanvas_updated.pdf")
canvas2.Print("RooFitCanvas_updated.png")
canvas3.Print("RooFitCanvasPass_updated.pdf")
canvas3.Print("RooFitCanvasPass_updated.png")

canvas4 = ROOT.TCanvas("efficiency", "efficiency", 800, 800)
hEff.Draw("COLZTEXTE")
canvas4.Print("Efficiency_updated.pdf")
canvas4.Print("Efficiency_updated.png")

canvas5 = ROOT.TCanvas("efficiency1D", "efficiency1D", 600, 600)
hEff1D.GetXaxis().SetTitle(" probe muon p_{T} [GeV/c]")
hEff1D.GetYaxis().SetTitle("Trigger Efficiency")
hEff1D.Draw()
canvas5.Print("Efficiency1D_updated.pdf")
canvas5.Print("Efficiency1D_updated.png")
