#!/bin/env python
from setdirs import *

"""
Collection of python methods implementing the BsMuMu MVA selection analysis
Produce materials (plots, tables) for the BsMuMu AN

methods implemented execute the following actions:

(part A)
- merges the input dataset trees
- split in barrel/endcaps regions
- split in event number (%3)
- applies pre-selection cuts (Common/Selection.h)
- adds TMVA MuonID information to trees, and applies MuonID selection (AddMuonID.C)
* input: data and MC trees
* output: trees ready to be passed on to TMVA training (runTMVA.sh)
* note: all rootfiles are saved under specified directory <rootDir>

(part B)
- run the TMVA training/optimization
  do the BDT training and produce the TMVA*.root files needed by the ensuing methods

(part C)
using previous inputs:
+ the "merged" data and MC root trees (former "all.py")
+ the TMVA results, obtained by executing "source runTMVA.sh" (which will submit 8 training in parallel)
does
- production of selection results, plots and tables

The analysis steps that need to be executed before are:
 - produce the PATtuples using the Onia2MuMu package; the PATtuples are currently applying preselection cuts, but no muon-id and trigger
 - produce the trees from the PATtuples; no additional cuts are applied

"""

### part A
################################## from previous TMVAPlots.py #####################
#Prepare plots and tables for the BsMuMu AN
#Collection of scripts to produce plots and tables for the BsMuMu analysis documentation.
#Needed input include:
#+ the "merged" data and MC root trees -- this is produced with "python all.py"
#+ the TMVA results, obtained by executing "source runTMVA.sh" (which will submit 8 training in parallel)
#   (note: it is better to first compile the macro TMVAClassification.C with "root -l"; ".L TMVAClassification.C+"; ".q")
# Summary of all previous analysis step:
# - Produce the PATtuples using the Onia2MuMu package. The PATtuples are currently applying preselection cuts, but no muon-id and trigger.
# - Produce the trees from the PATtuples. The trees apply no additional cut is applied.
# - Run "python all.py" to merge the data trees and then split them in the eventNum%3 subsamples. It also applies the trigger and muon-id requirements from Common/Selection.h.
# - Run "runTMVA.sh" to do the BDT training and produce the TMVA*.root files needed by the macros in this page.

# Remove the statbox
TROOT.gStyle.SetOptStat(0);

def saveCorrelationMatrixPlots(inputFile, region):
    print "saveCorrelationMatrixPlots(",inputFile,",",region,")"
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
    inputFile = TFile(rootDir+"TMVA_"+region+".root")
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

# Build ranking tables
def prepareTable(region, rankingType):
    outputFileBase = open(tablesDir+"ranking_"+region+".txt", "w")

    outputFileBase.write("\\begin{table} \n")
    outputFileBase.write("  \\begin{tabular}{|l|l|l|} \n")
    outputFileBase.write("    \\hline \n")
    outputFileBase.write("    rank & variable & separation \\\\ \n")
    outputFileBase.write("    \\hline \n")

    rankingFound = False
    for line in open(logsDir+"TMVALog_"+region+".txt"):
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
    print "created table:", outputFileBase.name

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

    file0 = open(logsDir+"TMVALog_"+region+"_0.txt")
    file1 = open(logsDir+"TMVALog_"+region+"_1.txt")
    file2 = open(logsDir+"TMVALog_"+region+"_2.txt")
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
        outputFileBase.write("  \\caption{Variable ranking for events of the three different event samples in the "+region+" before training.} \n")
    else:
        outputFileBase.write("  \\caption{Variable ranking for events of the three different event samples in the "+region+" after "+rankingType+" training.} \n")
    outputFileBase.write("\\end{table} \n")
    outputFileBase.close()
    print "created table:", outputFileBase.name


# Extract number of entries in each tree (= number of candidates after preselection)
def numEvents(fileName, dirName = ""):
    inputFile = TFile(fileName, "READ")
    h = 0
    if dirName != "":
        h = inputFile.Get(dirName).Get("probe_tree").GetEntries()
    else:
        h = inputFile.Get("probe_tree").GetEntries()
    # print fileName+" entries = "+str(h)
    return h

def writeSamples(file, region):
    numEvents(rootDir+region+"_preselection.root")  ##this call does nothing !!! please fix
    for index in range(0,3):
        num = numEvents(rootDir+region+"_preselection_"+str(index)+".root")
        file.write(str(num))
        if index < 2:
            file.write(" & ")
        else:
            file.write(" \\\\ \n")

# Prepare table
def getNumEvtTable():
    tableFile = open(tablesDir+"numCandidates.txt", "w")
    tableFile.write("\\begin{table} \n")
    tableFile.write("  \\begin{tabular}{|l|l|l|l|} \n")
    tableFile.write("    \\hline \n")
    tableFile.write("    Sample & Type 0 & Type 1 & Type 2 \\\\ \n")
    tableFile.write("    \\hline \n")
    
    # print "Signal MC"
    tableFile.write("    Signal barrel & ")
    writeSamples(tableFile, "BsMC12_barrel")
    tableFile.write("    Signal endcaps & ")
    writeSamples(tableFile, "BsMC12_endcaps")
    
    # print "Sidebands data"
    tableFile.write("    Background barrel & ")
    writeSamples(tableFile, "Barrel")
    tableFile.write("    Background endcaps & ")
    writeSamples(tableFile, "Endcaps")

    # Closing the table
    tableFile.write("    \\hline \n")
    tableFile.write("  \\end{tabular} \n")
    tableFile.write("  \\label{tab:numCandidates} \n")
    tableFile.write("  \\caption{Number of events per type for signal and background events in the barrel and endcap.} \n")
    tableFile.write("\\end{table} \n")
    tableFile.close()
    print "created table:", tableFile.name

# <headingcell level=2>

# Extract the BDT control plots


def saveBDTControlPlots(region):
    print rootDir+"TMVA_"+region+".root"
    inputFile = TFile(rootDir+"TMVA_"+region+".root")
    print inputFile
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


# Not clear where to get the information for the plots in figure 12 of the AN. Skipping them for now.
# Preselection efficiency and expected number of signal events after preselection

# This is required in order to estimate the expected number of signal events. This number is used together with the estimate of the background in the TMVA to extract the maximum significance (using a S/sqrt(S+B) figure of merit) and obtain the cut value.
# An analyzer can run on the PATtuples and extract the number of events processed and those passing the filters. It saves these numbers in two histograms with a single bin in a root file.
# The following macro extracts these numbers.


# Bs->mumu branching ratio from SM (see e.g.: http://www.zora.uzh.ch/49239/1/1103.2465v2-1.pdf): (0.32 +/- 0.02)10^-8
# Generated cross section 2082.8 with a filter efficiency of 0.002514 (see here: http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=BsToMuMu_BsFilter_8TeV-pythia6-evtgen).
# sigma = 2082.8
# This sigma is taken from http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=BsToMuMu_8TeV-pythia6-evtgen. This sample does not have BsFilter in the name, but the filter
# efficiency is the same and the cross section is much smaller. Using this the numbers obtained are more in line with expectations.
# NEED TO UNDERSTAND WICH ONE IS THE CORRECT CROSS SECTION
# From the paper: http://cds.cern.ch/record/1360175/files/PhysRevD.84.052008.pdf?version=1 -> ~ Bs*B(Bs->J/psi phi) = 6.9 +/- 0.6 +/- 0.6 nb
# From pdg B(Bs->J/psi phi) = (1.09 + 0.28 - 0.23)*10^-3
# And the expected B(Bs->mumu) = 0.32*10^-8
# Thus the cross section x branching ratio of Bs->mumu at 7TeV is 6.9*0.32/1.09 * 10^-5 nb = 2.03 * 10^-5 nb = 2.03 * 10^-2 pb
# THIS DOES NOT WORK


# Another indication is coming from the MC sample. In table 3 of the reference AN the MC sample is stated to be equivalent to 7405.3/fb.
# The number of entries in the MC file after all preselections (including trigger and muon-id) is of ~ 30000 for the barrel. When scaled
# to the integrated luminosity of the data it agrees better with the values obtained with the 173.57/pb sigma value.

def estimateSignal(region):

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
    # integratedLuminosity = runA+runARecover+runB+runC1+runCEcalRecover+runC2
    # In /pb
    integratedLuminosity = 18620

    expectedEvents = sigma*genFilterEff*integratedLuminosity

    printit = False
    if printit:
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
    if printit:
        print "processed events =", int(totEvents)
        print "events passing filters =", int(prefilterEvents)
        print "prefilter efficiency =", prefilterEfficiency

    # Take number of signal events from tree before trigger and muon-id
    preselectedEvents = numEvents(countersDir+"selection_test.root", "detailedDimuonTree")
    selectedSignalBarrel = numEvents(rootDir+"BsMC12_barrel_preselection.root")
    selectedSignalEndcaps = numEvents(rootDir+"BsMC12_endcaps_preselection.root")
    effBarrel = selectedSignalBarrel/float(preselectedEvents)
    effEndcaps = selectedSignalEndcaps/float(preselectedEvents)
    if printit:
        print "efficiency in the barrel after trigger and muon-id =", effBarrel*prefilterEfficiency
        print "efficiency in the endcaps after trigger and muon-id =", effEndcaps*prefilterEfficiency
    
    expectedEventsAfterPrefilterOnly = expectedEvents*prefilterEfficiency
    expectedEventsBarrel = expectedEventsAfterPrefilterOnly*effBarrel
    expectedEventsEndcaps = expectedEventsAfterPrefilterOnly*effEndcaps
    
    if printit:
        print "-----------------------------------------------------"
        print "Expected signal events in the barrel =", expectedEventsBarrel
        print "Expected signal events in the endcaps =", expectedEventsEndcaps
    
    #expected = []
    #expected.append(expectedEventsBarrel)
    #expected.append(expectedEventsEndcaps)
    #return expected

    expected = {"barrel":expectedEventsBarrel, "endcaps":expectedEventsEndcaps}    
    print region+" estimated signal =", expected[region]
    return expected[region]


# Background estimation
# To estimate the background take the number of events from the sidebands, average them, and normalize them to the area of the signal region.
def estimateBackground(region):
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    inputFile = TFile(rootDir + Region+"_preselection.root")
    inputDir = inputFile.Get("probe_tree")
    canvas = TCanvas()
    inputDir.Draw("m>>m")
    massPlot = TROOT.gROOT.FindObject("m")
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


# With these signal and background estimated the best cuts are:
# - for the barrel 0.1361
# - for the endcaps 0.2163
# Extract the numbers for the optimal cut value. This requires the previous steps for signal and background estimates to have already run.


#note: deprecated !!!
def optimalCutBDT(TMVAFileName, expectedSignalEvents, estimatedBackground):
    p = subprocess.Popen([rootExecutable, "-b", "-q", "-l", "mvaeffs.C+(\""+TMVAFileName+"\","+str(expectedSignalEvents)+","+str(estimatedBackground)+")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.splitlines():
        if line.find("Classifier  ") != -1:
            print line
        if line.find("    BDT") != -1:
            print line
            return line.split()[5]

#region = "barrel"
#print region.title()+":"
#optimalCutBarrel = optimalCutBDT(region, expectedEventsBarrel, estimatedBackgroundBarrel)
#os.system("mv mvaeffs_BDT_"+region+".pdf "+figuresDir)
#region = "endcaps"
#print region.title()+":"
#optimalCutEndcaps = optimalCutBDT(region, expectedEventsEndcaps, estimatedBackgroundEndcaps)
#os.system("mv mvaeffs_BDT_"+region+".pdf "+figuresDir)


# Extract the optimal BDT cuts for each BDT in the subsamples:
# - merge the TMVA output files for each BDT training
# - run the significance macro to extract the maximum significance
# Since the significance macro only uses the number of entries to recompute internally the efficiencies it is possible to merge the input files.

#get list of significance figures-of-merit
def getSignFomName(file):
    signifn = [] 
    maxSigFile = open(file)
    line = maxSigFile.readline()
    for i in range(3):
        print i, 2*i
        signifn.append(line.split()[2*i].replace(":",""))
        maxSigFile.close()
    return signifn


def printMvaCut(arr,regions,methods):
    signi_fom = getSignFomName(logsDir+"maxsignificance_BDT_barrel.txt")
    print "optimal cuts:",
    for rr in regions:
        print "\n\t",rr,
        for mm in methods:
            print "\n\t\t",mm,
            for ii in range(3):
                print "\t",signi_fom[ii],":% 5.3f" % float(arr[rr][mm][ii]), 



# BDT Application Overlay Plots
# Change this to redo also the BDT application #

class TypeSettings:
    def __init__(self, inputIndex, inputColor, inputOption):
        self.index = inputIndex
        self.color = inputColor
        self.option = inputOption

def applyMVA(inputFileName, outputFileName, weightDir, method, cutValue): 
    inputFileName  = rootDir    + inputFileName
    outputFileName = rootDir    + outputFileName
    weightDir      = weightsDir + weightDir
    print "applying MVA selection ", method, ">",cutValue," input:", inputFileName, " output:", outputFileName
    os.system(rootExecutable + " -q -l TMVAClassificationApplication.C+\(\\\""+inputFileName+"\\\",\\\""+outputFileName+"\\\",\\\""+weightDir+"\\\","+str(cutValue)+",\\\""+method+"\\\"\)")
    #p = subprocess.Popen([rootExecutable, "-b", "-q", "-l", "TMVAClassificationApplication.C+(\""+inputFileName+"\",\""+outputFileName+"\",\""+weightDir+"\","+str(cutValue)+",\""+method+"\")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #out, err = p.communicate()
    # print out
    # print err

## note: method obsolete: replaced by similar (BDR->MVA) method
def applyBDT(inputFileName, outputFileName, weightDir, cutValue):
    p = subprocess.Popen([rootExecutable, "-b", "-q", "-l", "TMVAClassificationApplication.C+(\""+inputFileName+"\",\""+outputFileName+"\",\""+weightDir+"\","+str(cutValue)+")"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    # print out
    print err

def drawAppMVAOutputPlots(region,method,isMC):
    print "executing drawAppMVAOutputPlots ", region, " ", method, " ", isMC
    MCstr = ""
    if isMC:
       MCstr = "BsMC12"
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    HistName = "ApplicationOutput"+MCstr+Region+method
    CvsName  = "canvasMVA_"+MCstr+method+region
    canvas = TCanvas(CvsName,CvsName)
    stackBDT = THStack(HistName,HistName)
    histos = []
    appFiles = []
    canvas.Draw()

    # note: TMVApp file name structure: "TMVApp"+Region+method+sample+".root"
    #                             "BsMC12TMVApp"+Region+method+sample+".root"
    appFile = TFile(rootDir+MCstr+"TMVApp"+Region+method+".root")
    histo = appFile.Get("MVA_"+method).Clone(HistName)
    histo.Scale(1/histo.GetEntries())
    stackBDT.Add(histo)

    for sample in range(3):
        appFiles.append(TFile(rootDir+MCstr+"TMVApp"+Region+method+str(sample)+".root"))
        histos.append(appFiles[sample].Get("MVA_"+method).Clone("HistName"+str(sample)))
        histos[sample].SetLineColor(2*sample)
        histos[sample].Scale(1/histos[sample].GetEntries())
        stackBDT.Add(histos[sample])


    stackBDT.Draw("nostack")
    if isMC:
        applicationBDTLegend = TLegend(0.2,0.7,0.5,0.9,"","brNDC")
        #applicationBDTLegend.SetHeader("Bs MC "+region.split("BsMC")[1])
        applicationBDTLegend.SetHeader("Bs MC "+Region)
    else:
        applicationBDTLegend = TLegend(0.55,0.7,0.85,0.9,"","brNDC")
        applicationBDTLegend.SetHeader(Region)

    applicationBDTLegend.AddEntry(histo, "Full sample", "l")
    applicationBDTLegend.AddEntry(histos[0], "Trained on 0, tested on 1, applied on 2", "l")
    applicationBDTLegend.AddEntry(histos[1], "Trained on 1, tested on 2, applied on 0", "l")
    applicationBDTLegend.AddEntry(histos[2], "Trained on 2, tested on 0, applied on 1", "l")
    applicationBDTLegend.Draw("same")
    applicationBDTLegend.SetFillColor(0)
    applicationBDTLegend.SetLineColor(0)
    canvas.SaveAs(figuresDir+"Application"+method+"Output_"+MCstr+region+".pdf")


## note: method obsolete: replaced by similar (BDT->MVA) method
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
    appFile = TFile(rootDir+region+"TMVApp.root")
    histo = appFile.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region)
    histo.Scale(1/histo.GetEntries())
    stackBDT.Add(histo)

    appFile0 = TFile(rootDir+region+"TMVApp0.root")
    histo0 = appFile0.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+"0")
    histo0.SetLineColor(2)
    histo0.Scale(1/histo0.GetEntries())
    stackBDT.Add(histo0)

    appFile1 = TFile(rootDir+region+"TMVApp1.root")
    histo1 = appFile1.Get("MVA_BDT").Clone("ApplicationBDTOutput_"+region+"1")
    histo1.SetLineColor(8)
    histo1.Scale(1/histo1.GetEntries())
    stackBDT.Add(histo1)

    appFile2 = TFile(rootDir+region+"TMVApp2.root")
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
    canvas.SaveAs(figuresDir+"ApplicationBDTOutput_"+region+".pdf")



#drawAppBDTOutputPlots("Barrel")
#drawAppBDTOutputPlots("BsMCBarrel")
#drawAppBDTOutputPlots("Endcaps")
#drawAppBDTOutputPlots("BsMCEndcaps")

# Draw the mass plot combining the results of the three BDTs
# note: TMVApp file name structure: "TMVApp"+Region+method+sample+".root"
#                             "BsMC12TMVApp"+Region+method+sample+".root"

def drawMassPlot(region, method, plotType = ""):
    print "drawMassPlot ", region, "", method
    # Take the histograms from each application
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    appFile0 = TFile(rootDir+"TMVApp"+Region+method+"0"+plotType+".root")
    appFile1 = TFile(rootDir+"TMVApp"+Region+method+"1"+plotType+".root")
    appFile2 = TFile(rootDir+"TMVApp"+Region+method+"2"+plotType+".root")
    histoMass0 = appFile0.Get("mass")
    histoMass1 = appFile1.Get("mass")
    histoMass2 = appFile2.Get("mass")
    histoMass = histoMass0.Clone(histoMass0.GetName()+"_merged")
    histoMass.Add(histoMass1)
    histoMass.Add(histoMass2)
    #print "Mass entries "+region+" "+str(histoMass.GetEntries())
    #print "Mass integral "+region+" "+str(histoMass.Integral())
    # Draw the merged histogram
    canvasMass = TCanvas("canvasMass"+region, "Mass "+region)
    canvasMass.Draw()
    histoMass.Draw()
    #histoMass.GetYaxis().SetRangeUser(0, 10)
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
    canvasMass.SaveAs(figuresDir+region+method+"MassPlot"+plotType+".pdf")

def drawVariablePlots(region, method, plotType = ""):
    Region = region[:1].upper()+region[1:] #capitalize first letter for retrieving file names
    appFile0 = TFile(rootDir+"TMVApp"+Region+method+"0"+plotType+".root")
    appFile1 = TFile(rootDir+"TMVApp"+Region+method+"1"+plotType+".root")
    appFile2 = TFile(rootDir+"TMVApp"+Region+method+"2"+plotType+".root")
    for var in ["pt", "eta", "fls3d", "alpha", "maxdoca", "pvip", "pvips", "iso", "docatrk", "closetrk", "chi2dof"]:
        canvas = TCanvas("canvas"+var+region+plotType)
        canvas.Draw()
        varPlot = appFile0.Get(var).Clone(var+"_merged")
        varPlot.Add(appFile1.Get(var))
        varPlot.Add(appFile2.Get(var))
        varPlot.Draw()
        variablePdfFileName = "AfterCut_"+var+"_"+region+method+plotType+".pdf"
        canvas.SaveAs(figuresDir+variablePdfFileName)



### part B

def runClassification(methodlist, eventsToTrain):

    print "running tmva classification, for methods ", methodlist
    os.system(rootExecutable+" -l -b -q compileTMVA.C")
    
    for region in regions:
        for isplit in range(-1,3):
            cmd = '\\\"'+eventsToTrain+'\\\",\\\"'+region+'\\\",\\\"'
            log = logsDir+"TMVALog_"+region
            if isplit != -1: 
                cmd += str(isplit)
                log += "_"+str(isplit)
            cmd+= '\\\",\\\"' + methodlist +  '\\\"'
            log+=".txt"
            print cmd,log
            if os.path.exists(log):os.system("rm "+log)
            #os.system(rootExecutable + " -l -b -q TMVAClassification.C\("+cmd+"\) >& "+log+" &")
            os.system(rootExecutable+" -l -b -q TMVAClassification.C+\("+cmd+"\) >& "+log+" &")

    #root -b -q -l compileTMVA.C
    #root -b -q -l TMVAClassification.C\(\"barrel\",\"\",\"BDT,MLP,CutsSA\"\)   >& logs/TMVALog_barrel.txt    &
    #root -b -q -l TMVAClassification.C\(\"barrel\",\"0\",\"BDT,MLP,CutsSA\"\)  >& logs/TMVALog_barrel_0.txt  &
    #root -b -q -l TMVAClassification.C\(\"barrel\",\"1\",\"BDT,MLP,CutsSA\"\)  >& logs/TMVALog_barrel_1.txt  &
    #root -b -q -l TMVAClassification.C\(\"barrel\",\"2\",\"BDT,MLP,CutsSA\"\)  >& logs/TMVALog_barrel_2.txt  &
    #root -b -q -l TMVAClassification.C\(\"endcaps\",\"\",\"BDT,MLP,CutsSA\"\)  >& logs/TMVALog_endcaps.txt   &
    #root -b -q -l TMVAClassification.C\(\"endcaps\",\"0\",\"BDT,MLP,CutsSA\"\) >& logs/TMVALog_endcaps_0.txt &
    #root -b -q -l TMVAClassification.C\(\"endcaps\",\"1\",\"BDT,MLP,CutsSA\"\) >& logs/TMVALog_endcaps_1.txt &
    #root -b -q -l TMVAClassification.C\(\"endcaps\",\"2\",\"BDT,MLP,CutsSA\"\) >& logs/TMVALog_endcaps_2.txt &


### part C
################################## from previous all.py #####################
#- merges the input dataset trees
#- split in barrel/endcaps regions
#- split in event number (%3)
#- applies pre-selection cuts (Common/Selection.h)
#- adds TMVA MuonID information to trees, and applies MuonID selection (AddMuonID.C)
#* input: data and MC trees
#* output: trees ready to be passed on to TMVA training (runTMVA.sh)
#* note: all rootfiles are saved under specified directory <rootDir>


#splitting = 0,1,2  is the value of the rest of event number % 3
#            -1     denotes the full (unsplitted) sample


#def main():
#
#    # Run the selection and splitting
#    for isplit in range(-1,3):
#
#        appendName = "_preselection.root"
#        if isplit != -1:            
#            appendName = appendName.split(".")[0] +"_"+str(isplit)+".root"
#
#        print "processing blinded sample", isplit 
#        applySelectionAndSplit(inputTrees, isplit, maxRun, True)   
#        combineSamples(appendName)
#        #addMuonID(appendName)
#
#        print "processing unblinded sample", isplit 
#        applySelectionAndSplit(inputTrees, isplit, maxRun, False)
#        appendName = appendName.split(".")[0] +"_unblinded"+".root"
#        combineSamples(appendName)
#        #addMuonID(appendName)

#def combineSamplesAndAddMuonID(appendName):
    #combineSamples(appendName)
    ### TBD: DISABLED TEMPORARILY!!!
    #addMuonID(appendName)

def addMuonID(appendName):
    # Add -- and Apply --  MVA muon-id

    for region in regions:

        Region = region[:1].upper()+region[1:]  # capitalize
        fname = rootDir + Region + appendName

        os.system(rootExecutable+" -l -b -q AddMuonID.C+\(\\\""+fname+"\\\"\)")
            #os.system(rootExecutable+" -l -b -q AddMuonID.C+\(\\\"Endcaps"+appendName+"\\\"\)")
        os.system("mv " + fname + "_muonID.root " + fname)
            #os.system("mv Endcaps"+appendName+"_muonID.root Endcaps"+appendName)
    
        # Process MC
        if appendName.find("unblinded") == -1:
            fname = rootDir + "BsMC12_" + region + appendName
                #os.system("mv BsMC_barrel"+appendName+" BsMC12_barrel"+appendName)
                #os.system("mv BsMC_endcaps"+appendName+" BsMC12_endcaps"+appendName)
            os.system(rootExecutable+" -l -b -q AddMuonID.C+\(\\\""+fname+"\\\"\)")
                #os.system(rootExecutable+" -l -b -q AddMuonID.C+\(\\\"BsMC12_endcaps"+appendName+"\\\"\)")
            os.system("mv "+fname+"_muonID.root " + fname)
                #os.system("mv BsMC12_endcaps"+appendName+"_muonID.root BsMC12_endcaps"+appendName)
        

def combineSamples(appendName):

    print "combineSamples(",appendName,")"

    # Combine barrel, endcaps samples
    for region in regions:

        Region = region[:1].upper()+region[1:]  # capitalize
        fname = rootDir + Region + appendName
        if os.path.exists(fname):
            os.system("rm -f "+fname)

        addlist = ""
        for run in runlist:
            addlist += rootDir + run +"_" + region + appendName + " "

        print "hadd-ing files:", addlist,  " created:", fname
        print "hadd "+ fname + " " + addlist
        os.system("hadd "+ fname + " " + addlist)

        #rename MC samples
        if not "unblind" in appendName:
            fname = rootDir + "BsMC_" + region + appendName
            os.system("mv " + fname + " " + str(fname.replace("BsMC","BsMC12")))

    os.system("rm -f "+rootDir+"Run*.root")



def applySelectionAndSplit(inputTrees, splitting, maxRun, blindData = True, cut_based=False):

    print "applySelectionAndSplit(",inputTrees,",", splitting,",", maxRun,",", blindData,",",cut_based,")"

    splitString =""
    if splitting != -1:
        splitString = "_"+str(splitting)
    
    appendName = "_preselection"+splitString+".root"
    # appendNameMC = "_preselection.root"
    
    for tree in inputTrees:
        #cut_based = False
        data = True
        if tree.find("MC") != -1: data = False
        # if data == "0": append = appendNameMC
        # NO SPACES IN THE ROOT COMMAND
        for regidx, region in enumerate(regions):
            outputTree = tree.split("/")[-2]+'_'+region
            if not cut_based:
                # print "Applying preselection cuts"
                outputTree += appendName
            else:
                # print "Applying analysis cuts"
                outputTree += '.root'
            blinding = False
            if data:
                if blindData:
                    blinding = True
                else:
                    outputTree = outputTree.split(".")[0]+"_unblinded.root"
            if cut_based:
                    outputTree = outputTree.split(".")[0]+"_MainCNCSelected.root"                
            outputTree = rootDir + outputTree
            # print 'root -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regionIndex)+','+data+','+cut_based+','+blinding+','+str(splitting)+','+maxRun+'\)'
            # print "applying selection and creating:", outputTree
            os.system(rootExecutable+' -l -b -q cutTree_BsMuMu.C\(\\"'+tree+'\\",\\"'+outputTree+'\\",'+str(regidx)+','+str(data)+','+str(cut_based)+','+str(blinding)+','+str(splitting)+','+maxRun+'\)')

