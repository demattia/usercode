# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Prepare plots and tables for the BsMuMu AN

# <rawcell>

# Collection of scripts to produce the plots and tables for the BsMuMu AN.
# If the root trees for data and MC have not been merged yet run "python all.py" to produce them.
# Use the runTMVA.sh script to produce all the necessary files. To run the script do: "source runTMVA.sh" (It will submit 8 BDT trainings in parallel).
# It is better to first compile the macro TMVAClassification.C with "root -l"; ".L TMVAClassification.C+"; ".q".
# The first input is the cd into the directory where all the files are located.
# 
# Summary of all the steps:
# - Produce the PATtuples using the Onia2MuMu package. The PATtuples are currently applying preselection cuts, but no muon-id and trigger.
# - Produce the trees from the PATtuples. The trees apply no additional cut is applied.
# - Run "python all.py" to merge the data trees and then split them in the eventNum%3 subsamples. It also applies the trigger and muon-id requirements from Common/Selection.h.
# - Run "runTMVA.sh" to do the BDT training and produce the TMVA*.root files needed by the macros in this page.

# <headingcell level=2>

# Produce TMVA plots and copy them to latex figures dir

# <codecell>

# Definition of common variables
figuresDir = "/Users/demattia/TMVA-v4.1.2/test/NewTrees/BsMuMuLatex/Figures/"
tablesDir = "/Users/demattia/TMVA-v4.1.2/test/NewTrees/BsMuMuLatex/Tables/"

# <codecell>

# cd /Users/demattia/TMVA-v4.1.2/test/NewTrees/

# <codecell>

from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import THStack

def saveCorrelationMatrixPlots(inputFile, region):
    correlationMatrixS = inputFile.Get("CorrelationMatrixS")
    correlationMatrixB = inputFile.Get("CorrelationMatrixB")
    matrixCanvasS = TCanvas("correlationMatrixS", "correlation matrix signal", 600, 500)
    matrixCanvasS.Draw()
    correlationMatrixS.Draw("BOXCOLTEXT");
    matrixCanvasS.SaveAs(figuresDir+"correlationMatrixS_"+region+".pdf")
    matrixCanvasB = TCanvas("correlationMatrixB", "correlation matrix background", 600, 500)
    matrixCanvasB.Draw()
    correlationMatrixB.Draw("BOXCOLTEXT");
    matrixCanvasB.SaveAs(figuresDir+"correlationMatrixB_"+region+".pdf")

def TMVAPlots(region):
    inputFile = TFile("TMVA_"+region+".root")
    saveCorrelationMatrixPlots(inputFile, region)
    histoDir = inputFile.Get("InputVariables_Id")
    histoList = histoDir.GetListOfKeys()
    for key in histoList:
        histo = key.ReadObj()
        if isinstance(histo, TH1F):
            fullName = histo.GetName()
            name = fullName.split("__")[0]
            hType = fullName.split("__")[1].split("_")[0]
            if hType == "Signal":
                canvas = TCanvas(name+"_canvas")
                stack = THStack(name+"_stack", name)
                canvas.Draw()
                histo.Scale(1/histo.Integral())
                stack.Add(histo)
                histo.SetLineColor(4)
                histo.SetFillColor(4);
                histo.SetFillStyle(3004);
                leg = TLegend(0.5536913,0.770979,0.9345638,0.9702797,"","brNDC")
                leg.AddEntry(histo, hType, "f")
            else:
                histo.Scale(1/histo.Integral())
                stack.Add(histo)
                stack.Draw("nostack")
                stack.SetMaximum(stack.GetMaximum("nostack")*1.35)
                stack.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
                stack.GetYaxis().SetTitle(histo.GetYaxis().GetTitle())
                histo.SetLineColor(2)
                histo.SetFillColor(2);
                histo.SetFillStyle(3005);
                leg.AddEntry(histo, hType, "f")
                leg.SetFillColor(0)
                leg.Draw()
                canvas.SaveAs(figuresDir+name+"_"+region+".pdf")

TMVAPlots("barrel")
TMVAPlots("endcaps")

# <headingcell level=3>

# Build ranking tables

# <codecell>

def prepareTable(region, rankingType):
    outputFileBase = open(tablesDir+"ranking_"+region+".txt", "w")

    outputFileBase.write("\\begin{table} \n")
    outputFileBase.write("  \\begin{tabular}{|l|l|l|} \n")
    outputFileBase.write("    \\hline \n")
    outputFileBase.write("    rank & variable & separation \\\\ \n")
    outputFileBase.write("    \\hline \n")

    rankingFound = False
    for line in open("TMVALog_"+region+".txt"):
        if rankingFound and line.find("Factory") != -1:
            rankingFound = False
        if rankingFound:
            if line.find("----") == -1:
                splittedLine = line.split()
                rank = splittedLine[3]
                variable = splittedLine[5]
                separation = splittedLine[7]
                if line.find(rankingType) != -1:
                    outputFileBase.write("    "+rank+" & "+variable+" & "+separation+" \\\\ \n")
        if line.find("Rank ") != -1:
            rankingFound = True

    outputFileBase.write("    \\hline \n")
    outputFileBase.write("  \\end{tabular} \n")
    outputFileBase.write("  \\label{tab:datasets} \n")
    outputFileBase.write("\\end{table} \n")
    outputFileBase.close()
    
prepareTable("barrel", "IdTransformation")
prepareTable("endcaps", "IdTransformation")


import itertools

def buildTableLine(line):
    splittedLine = line.split()
    variable = splittedLine[5]
    separation = splittedLine[7]
    return variable+" & "+separation
    
def prepareMultiTable(region, rankingType):
    outputFileBase = open(tablesDir+"ranking_"+region+"_"+rankingType+".txt", "w")

    outputFileBase.write("\\begin{table} \n")
    outputFileBase.write("  \\begin{tabular}{|l|l|l|l|l|l|l|} \n")
    outputFileBase.write("    \\hline \n")
    outputFileBase.write("    & \multicolumn{2}{|c|}{0} & \multicolumn{2}{|c|}{1} & \multicolumn{2}{|c|}{2} \\\\ \n")
    outputFileBase.write("    \\hline \n")
    outputFileBase.write("    rank & variable & separation & variable & separation & variable & separation \\\\ \n")
    outputFileBase.write("    \\hline \n")

    file0 = open("TMVALog_"+region+"_0.txt")
    file1 = open("TMVALog_"+region+"_1.txt")
    file2 = open("TMVALog_"+region+"_2.txt")
    rankingFound = False
    for line1, line2, line3 in itertools.izip(file0, file1, file2):
        if rankingFound and line1.find("Factory") != -1:
            rankingFound = False
        if rankingFound:
            if line1.find("----") == -1:
                if line1.find(rankingType) != -1:
                    rank = line1.split()[3]
                    outputFileBase.write("    "+rank+" & ")
                    outputFileBase.write(buildTableLine(line1)+" & ")
                    outputFileBase.write(buildTableLine(line2)+" & ")
                    outputFileBase.write(buildTableLine(line3)+" \\\\ \n")
        if line1.find("Rank ") != -1:
            rankingFound = True

    outputFileBase.write("    \\hline \n")
    outputFileBase.write("  \\end{tabular} \n")
    outputFileBase.write("  \\label{tab:datasets"+region+rankingType+"} \n")
    if rankingType == "IdTransformation":
        outputFileBase.write("  \\caption{Variable ranking for events of the three different event samples in the "+region+" before BDT training.} \n")
    else:
        outputFileBase.write("  \\caption{Variable ranking for events of the three different event samples in the "+region+" after BDT training.} \n")
    outputFileBase.write("\\end{table} \n")
    outputFileBase.close()

prepareMultiTable("barrel", "IdTransformation")
prepareMultiTable("endcaps", "IdTransformation")

prepareMultiTable("barrel", "BDT")
prepareMultiTable("endcaps", "BDT")

# <headingcell level=2>

# Extract number of entries in each tree (= number of candidates after preselection)

# <codecell>

def numEvents(fileName, dirName = ""):
    inputFile = TFile(fileName, "READ")
    h = 0
    if dirName != "":
        h = inputFile.Get(dirName).Get("probe_tree").GetEntries()
    else:
        h = inputFile.Get("probe_tree").GetEntries()
    # print fileName+" entries = "+str(h)
    return h

def writeSamples(region):
    numEvents(region+"_preselection.root")
    for index in range(0,3):
        num = numEvents(region+"_preselection_"+str(index)+".root")
        tableFile.write(str(num))
        if index < 2:
            tableFile.write(" & ")
        else:
            tableFile.write(" \\\\ \n")

# Prepare table
tableFile = open(tablesDir+"numCandidates.txt", "w")
tableFile.write("\\begin{table} \n")
tableFile.write("  \\begin{tabular}{|l|l|l|l|} \n")
tableFile.write("    \\hline \n")
tableFile.write("    Sample & Type 0 & Type 1 & Type 2 \\\\ \n")
tableFile.write("    \\hline \n")

# print "Signal MC"
tableFile.write("    Signal barrel & ")
writeSamples("BsMC12_barrel")
tableFile.write("    Signal endcaps & ")
writeSamples("BsMC12_endcaps")

# print "Sidebands data"
tableFile.write("    Background barrel & ")
writeSamples("Barrel")
tableFile.write("    Background endcaps & ")
writeSamples("Endcaps")

# Closing the table
tableFile.write("    \\hline \n")
tableFile.write("  \\end{tabular} \n")
tableFile.write("  \\label{tab:numCandidates} \n")
tableFile.write("  \\caption{Number of events per type for signal and background events in the barrel and endcap.} \n")
tableFile.write("\\end{table} \n")
tableFile.close()

# This can be used to get the numbers for the full sample
# num = numEvents("BsMC12_barrel_preselection.root")
# num = numEvents("BsMC12_endcaps_preselection.root")
# num = numEvents("Barrel_preselection.root")
# num = numEvents("Endcaps_preselection.root")

# <headingcell level=2>

# Extract the BDT control plots

# <codecell>

from ROOT import TFile
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TROOT

# Remove the statbox
TROOT.gStyle.SetOptStat(0);

def saveBDTControlPlots(region):
    inputFile = TFile("TMVA_"+region+".root")
    dir = inputFile.Get("Method_BDT").Get("BDT")
    canvasBDT = TCanvas("BDTControlPlots"+region, "BDT control plots "+region, 600, 500)
    canvasBDT.Divide(1,3)
    canvasBDT.cd(1)
    histo1 = dir.Get("BoostWeightVsTree")
    histo1.SetLineColor(4)
    histo1.Draw()
    canvasBDT.cd(2)
    histo2 = dir.Get("ErrFractHist")
    histo2.SetLineColor(4)
    histo2.Draw()
    canvasBDT.cd(3)
    histo3 = dir.Get("NodesBeforePruning")
    histo3.SetLineColor(4)
    histo3.Draw()
    canvasBDT.SaveAs(figuresDir+"BDTControlPlots_"+region+".pdf")

saveBDTControlPlots("barrel")
saveBDTControlPlots("endcaps")

# <rawcell>

# Not clear where to get the information for the plots in figure 12 of the AN. Skipping them for now.

# <headingcell level=2>

# Preselection efficiency and expected number of signal events after preselection

# <rawcell>

# This is required in order to estimate the expected number of signal events. This number is used together with the estimate of the background in the TMVA to extract the maximum significance (using a S/sqrt(S+B) figure of merit) and obtain the cut value.
# An analyzer can run on the PATtuples and extract the number of events processed and those passing the filters. It saves these numbers in two histograms with a single bin in a root file.
# The following macro extracts these numbers.

# <codecell>

countersDir = "Trees/BsMC/"

from ROOT import TFile
from ROOT import TH1F

# Bs->mumu branching ratio from SM (see e.g.: http://www.zora.uzh.ch/49239/1/1103.2465v2-1.pdf): (0.32 +/- 0.02)10^-8
# Generated cross section 2082.8 with a filter efficiency of 0.002514 (see here: http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=BsToMuMu_BsFilter_8TeV-pythia6-evtgen).
# sigma = 2082.8
# This sigma is taken from http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=BsToMuMu_8TeV-pythia6-evtgen. This sample does not have BsFilter in the name, but the filter
# efficiency is the same and the cross section is much smaller. Using this the numbers obtained are more in line with expectations.
# NEED TO UNDERSTAND WICH ONE IS THE CORRECT CROSS SECTION
# From the paper: http://cds.cern.ch/record/1360175/files/PhysRevD.84.052008.pdf?version=1 -> ~ Bs*B(Bs->J/psi phi) = 6.9 +/- 0.6 +/- 0.6 nb
# From pdg B(Bs->J/psi phi) = (1.09 + 0.28 − 0.23)*10^−3
# And the expected B(Bs->mumu) = 0.32*10^-8
# Thus the cross section x branching ratio of Bs->mumu at 7TeV is 6.9*0.32/1.09 * 10^-5 nb = 2.03 * 10^-5 nb = 2.03 * 10^-2 pb
# THIS DOES NOT WORK


# Another indication is coming from the MC sample. In table 3 of the reference AN the MC sample is stated to be equivalent to 7405.3/fb.
# The number of entries in the MC file after all preselections (including trigger and muon-id) is of ~ 30000 for the barrel. When scaled
# to the integrated luminosity of the data it agrees better with the values obtained with the 173.57/pb sigma value.


# in /pb
sigma = 173.57
genFilterEff = 0.002514
# in /pb
runA = 805.116
runARecover = 82.136
runB = 4429.
runC1 = 495.003
runCEcalRecover = 134.242
runC2 = 6397.
runD = 7483.

# Full 2012 data
# integratedLuminosity = runA+runARecover+runB+runC1+runCEcalRecover+runC2+runD
# Same data as reference AN
integratedLuminosity = runA+runARecover+runB+runC1+runCEcalRecover+runC2

expectedEvents = sigma*genFilterEff*integratedLuminosity

print "Generated cross section =", sigma
print "Generator filter efficiency =", genFilterEff
print "Integrated luminosity =", integratedLuminosity
print "Expected events in "+str(integratedLuminosity/1000.)+"/fb =", expectedEvents

# Compute prefilter efficiency
inputFile = TFile(countersDir+"counters.root")
histoDir = inputFile.Get("counterReader")
totEvents = histoDir.Get("numProcessedEvents").GetBinContent(1)
prefilterEvents = histoDir.Get("numEventsPassingFilter").GetBinContent(1)
prefilterEfficiency = prefilterEvents/totEvents
print "processed events =", int(totEvents)
print "events passing filters =", int(prefilterEvents)
print "prefilter efficiency =", prefilterEfficiency

# Take number of signal events from tree before trigger and muon-id
preselectedEvents = numEvents("Trees/BsMC/selection_test.root", "detailedDimuonTree")
selectedSignalBarrel = numEvents("BsMC12_barrel_preselection.root")
selectedSignalEndcaps = numEvents("BsMC12_endcaps_preselection.root")
effBarrel = selectedSignalBarrel/float(preselectedEvents)
effEndcaps = selectedSignalEndcaps/float(preselectedEvents)
print "efficiency in the barrel after trigger and muon-id =", effBarrel*prefilterEfficiency
print "efficiency in the endcaps after trigger and muon-id =", effEndcaps*prefilterEfficiency

expectedEventsAfterPrefilterOnly = expectedEvents*prefilterEfficiency
expectedEventsBarrel = expectedEventsAfterPrefilterOnly*effBarrel
expectedEventsEndcaps = expectedEventsAfterPrefilterOnly*effEndcaps

print "-----------------------------------------------------"
print "Expected signal events in the barrel =", expectedEventsBarrel
print "Expected signal events in the endcaps =", expectedEventsEndcaps

# <headingcell level=2>

# Background estimation

# <rawcell>

# To estimate the background take the number of events from the sidebands, average them, and normalize them to the area of the signal region.

# <codecell>

from ROOT import TROOT

def estimateBackground(region):
    inputFile = TFile(region+"_preselection.root")
    inputDir = inputFile.Get("probe_tree")
    inputDir.Draw("mass>>mass")
    massPlot = TROOT.gROOT.FindObject("mass")
    # canvas = TCanvas()
    # canvas.Draw()
    lowSideband = massPlot.Integral(massPlot.FindBin(4.9), massPlot.FindBin(5.2))
    lowSidebandWidth = 5.2-4.9
    # highSideband = massPlot.Integral(massPlot.FindBin(5.45), massPlot.FindBin(5.9))
    # highSidebandWidth = 5.9-5.45
    highSideband = massPlot.Integral(massPlot.FindBin(5.45), massPlot.FindBin(5.75))
    highSidebandWidth = 5.75-5.45

    signalRegionWidth = 5.45-5.2

    # print lowSideband, lowSidebandWidth
    # print highSideband, highSidebandWidth

    # print signalRegionWidth

    # Estimate background events in the signal region
    A = lowSideband
    a = lowSidebandWidth
    b = signalRegionWidth
    C = highSideband
    c = highSidebandWidth

    # background = (A/a + C/c - (A+C)/(a+b+c))*b*(a+b+c)/(a+2*b+c)
    background = b*(A+C)/(a+c)
    # print A/a
    # print C/c
    # print background/b
    print region+" estimated background =", background
    return background

estimatedBackgroundBarrel = estimateBackground("Barrel")
estimatedBackgroundEndcaps = estimateBackground("Endcaps")

# <rawcell>

# With these signal and background estimated the best cuts are:
# - for the barrel 0.1361
# - for the endcaps 0.2163
# Extract the numbers for the optimal cut value. This requires the previous steps for signal and background estimates to have already run.

# <codecell>

import subprocess
def optimalCutBDT(TMVAFileName, expectedSignalEvents, estimatedBackground):
    p = subprocess.Popen(["/Users/demattia/Downloads/root-v5-34-00-patches/bin/root", "-q", "-l", "mvaeffs.C(\""+TMVAFileName+"\","+str(expectedSignalEvents)+","+str(estimatedBackground)+")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.splitlines():
        if line.find("Classifier  ") != -1:
            print line
        if line.find("    BDT") != -1:
            print line
            return line.split()[5]

region = "barrel"
print region.title()+":"
optimalCutBarrel = optimalCutBDT(region, expectedEventsBarrel, estimatedBackgroundBarrel)
os.system("mv mvaeffs_BDT_"+region+".pdf "+figuresDir)
region = "endcaps"
print region.title()+":"
optimalCutBarrel = optimalCutBDT(region, expectedEventsEndcaps, estimatedBackgroundEndcaps)
os.system("mv mvaeffs_BDT_"+region+".pdf "+figuresDir)

# <headingcell level=3>

# Apply the BDT to the data and obtain the mass plot.

# <codecell>

import subprocess
from ROOT import TLine
from ROOT import TGraphErrors
from array import array

def applyBDT(inputFileName, outputFileName, weightDir, cutValue):
    p = subprocess.Popen(["/Users/demattia/Downloads/root-v5-34-00-patches/bin/root", "-q", "-l", "TMVAClassificationApplication.C+(\""+inputFileName+"\",\""+outputFileName+"\",\""+weightDir+"\","+str(cutValue)+")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    # print out
    # print err

def drawMassPlot(region, plotType = ""):
    appFile = TFile(region+"TMVApp"+plotType+".root")
    histoMass = appFile.Get("mass")
    canvasMass = TCanvas("canvasMass"+region, "Mass "+region)
    canvasMass.Draw()
    histoMass.Draw()
    histoMass.GetYaxis().SetRangeUser(0, 10)
    # line1 = TLine(5.2,1,5.3,1)
    # line1.SetLineColor(2)
    # line1.SetLineWidth(2)
    # line1.Draw("same")
    # line2 = TLine(5.3,1.5,5.45,1.5)
    # line2.SetLineColor(4)
    # line2.SetLineWidth(2)
    # line2.Draw("same")
    x = array("d", [5.25])
    y = array("d", [0.6])
    errX = array("d", [0.05])
    errY = array("d", [0.])
    line1 = TGraphErrors(1, x, y, errX, errY)
    line1.Draw("E1")
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.SetMarkerColor(2)
    line1.SetMarkerSize(0)
    x = array("d", [5.375])
    y = array("d", [0.8])
    errX = array("d", [0.075])
    errY = array("d", [0.])
    line2 = TGraphErrors(1, x, y, errX, errY)
    line2.Draw("E1")
    line2.SetLineColor(4)
    line2.SetLineWidth(2)
    line2.SetMarkerColor(4)
    line2.SetMarkerSize(0)
    histoMass.GetXaxis().SetTitle("Mass [GeV/c^{2}]")
    histoMass.GetYaxis().SetTitle("Entries")
    massLegend = TLegend(0.5,0.7,0.8,0.9,"","brNDC")
    massLegend.SetHeader(region)
    massLegend.AddEntry(line1, "B^{0} Signal Region", "l")
    massLegend.AddEntry(line2, "B^{0}_{s} Signal Region", "l")
    massLegend.Draw("same")
    massLegend.SetFillColor(0)
    massLegend.SetLineColor(0)
    canvasMass.SaveAs(figuresDir+region+"MassPlot"+plotType+".pdf")

def drawVariablePlots(region, plotType = ""):
    appFile = TFile(region+"TMVApp"+plotType+".root")
    for var in ["pt", "eta", "fls3d", "alpha", "maxdoca", "pvip", "pvips", "iso", "docatrk", "closetrk", "chi2dof"]:
        canvas = TCanvas("canvas"+var+region+plotType)
        canvas.Draw()
        appFile.Get(var).Draw()
        variablePdfFileName = "AfterBDTCut_"+var+"_"+region+plotType+".pdf"
        canvas.SaveAs(figuresDir+variablePdfFileName)

        
# applyBDT("Barrel_preselection.root", "BarrelTMVApp.root", "barrelWeights/", optimalCutBarrel)
# drawMassPlot("Barrel")
# applyBDT("Endcaps_preselection.root", "EndcapsTMVApp.root", "endcapsWeights/", optimalCutBarrel)
# drawMassPlot("Endcaps")
drawVariablePlots("Barrel")
drawVariablePlots("Endcaps")

# applyBDT("Barrel_preselection_unblinded.root", "BarrelTMVAppUnblinded.root", "barrelWeights/", optimalCutBarrel)
# drawMassPlot("Barrel", "Unblinded")
# applyBDT("Endcaps_preselection_unblinded.root", "EndcapsTMVAppUnblinded.root", "endcapsWeights/", optimalCutBarrel)
# drawMassPlot("Endcaps", "Unblinded")
drawVariablePlots("Barrel", "Unblinded")
drawVariablePlots("Endcaps", "Unblinded")

# <rawcell>

# Figure 19

# <codecell>

class TypeSettings:
    def __init__(self, inputIndex, inputColor, inputOption):
        self.index = inputIndex
        self.color = inputColor
        self.option = inputOption

def drawAppBDTOutputPlots(region):
    canvas = TCanvas("canvasMVA_BDT"+region, "canvasMVA_BDT"+region)
    stackBDT = THStack("ApplicationBDTOutput_"+region, "ApplicationBDTOutput_"+region)
    # histo = []
    canvas.Draw()
    # The loop plots only one histogram. Probably the others are removed or overwritten
    # for plotType in [TypeSettings("", 1, ""), TypeSettings("0", 2, "same"), TypeSettings("1", 3, "same"), TypeSettings("2", 4, "same")]:
    #     appFile = TFile(region+"TMVApp"+plotType.index+".root")
    #     print plotType.option
    #     histo.append(appFile.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+plotType.index))
    #     stackBDT.Add(histo[len(histo)-1])
    #     # histo.Draw()
    appFile = TFile(region+"TMVApp.root")
    histo = appFile.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region)
    histo.Scale(1/histo.GetEntries())
    stackBDT.Add(histo)

    appFile0 = TFile(region+"TMVApp0.root")
    histo0 = appFile0.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+"0")
    histo0.SetLineColor(2)
    histo0.Scale(1/histo0.GetEntries())
    stackBDT.Add(histo0)

    appFile1 = TFile(region+"TMVApp1.root")
    histo1 = appFile1.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+"1")
    histo1.SetLineColor(8)
    histo1.Scale(1/histo1.GetEntries())
    stackBDT.Add(histo1)

    appFile2 = TFile(region+"TMVApp2.root")
    histo2 = appFile2.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+"2")
    histo2.SetLineColor(4)
    histo2.Scale(1/histo2.GetEntries())
    stackBDT.Add(histo2)

    stackBDT.Draw("nostack")

    if region.find("BsMC") != -1:
        applicationBDTLegend = TLegend(0.2,0.7,0.5,0.9,"","brNDC")
        applicationBDTLegend.SetHeader("Bs MC "+region.split("BsMC")[1])
    else:
        applicationBDTLegend = TLegend(0.55,0.7,0.85,0.9,"","brNDC")
        applicationBDTLegend.SetHeader(region)
    applicationBDTLegend.AddEntry(histo, "Full sample", "l")
    applicationBDTLegend.AddEntry(histo0, "Trained on 0, tested on 1, applied on 2", "l")
    applicationBDTLegend.AddEntry(histo1, "Trained on 1, tested on 2, applied on 0", "l")
    applicationBDTLegend.AddEntry(histo2, "Trained on 2, tested on 0, applied on 1", "l")
    applicationBDTLegend.Draw("same")
    applicationBDTLegend.SetFillColor(0)
    applicationBDTLegend.SetLineColor(0)

    # appFile = TFile(region+"TMVApp1.root")
    # appFile.Get("MVA_BDT").Draw()
    canvas.SaveAs(figuresDir+"ApplicationBDTOutput_"+region+".pdf")

# Barrel
# Applied on 2, trained on 0 (tested on 1)
# applyBDT("Barrel_preselection_2.root", "BarrelTMVApp0.root", "barrel0Weights/", optimalCutBarrel)
# Applied on 0, trained on 1 (tested on 2)
# applyBDT("Barrel_preselection_0.root", "BarrelTMVApp1.root", "barrel1Weights/", optimalCutBarrel)
# Applied on 1, trained on 2 (tested on 0)
# applyBDT("Barrel_preselection_1.root", "BarrelTMVApp2.root", "barrel2Weights/", optimalCutBarrel)

drawAppBDTOutputPlots("Barrel")

# applyBDT("BsMC12_barrel_preselection.root", "BsMCBarrelTMVApp.root", "barrelWeights/", optimalCutBarrel)
# applyBDT("BsMC12_barrel_preselection_2.root", "BsMCBarrelTMVApp0.root", "barrel0Weights/", optimalCutBarrel)
# applyBDT("BsMC12_barrel_preselection_0.root", "BsMCBarrelTMVApp1.root", "barrel1Weights/", optimalCutBarrel)
# applyBDT("BsMC12_barrel_preselection_1.root", "BsMCBarrelTMVApp2.root", "barrel2Weights/", optimalCutBarrel)

drawAppBDTOutputPlots("BsMCBarrel")

# Endcaps
# Applied on 2, trained on 0 (tested on 1)
# applyBDT("Endcaps_preselection_2.root", "EndcapsTMVApp0.root", "endcaps0Weights/", optimalCutEndcaps)
# Applied on 0, trained on 1 (tested on 2)
# applyBDT("Endcaps_preselection_0.root", "EndcapsTMVApp1.root", "endcaps1Weights/", optimalCutEndcaps)
# Applied on 1, trained on 2 (tested on 0)
# applyBDT("Endcaps_preselection_1.root", "EndcapsTMVApp2.root", "endcaps2Weights/", optimalCutEndcaps)

drawAppBDTOutputPlots("Endcaps")

# applyBDT("BsMC12_endcaps_preselection.root", "BsMCEndcapsTMVApp.root", "endcapsWeights/", optimalCutEndcaps)
# applyBDT("BsMC12_endcaps_preselection_2.root", "BsMCEndcapsTMVApp0.root", "endcaps0Weights/", optimalCutEndcaps)
# applyBDT("BsMC12_endcaps_preselection_0.root", "BsMCEndcapsTMVApp1.root", "endcaps1Weights/", optimalCutEndcaps)
# applyBDT("BsMC12_endcaps_preselection_1.root", "BsMCEndcapsTMVApp2.root", "endcaps2Weights/", optimalCutEndcaps)

drawAppBDTOutputPlots("BsMCEndcaps")

# <codecell>


