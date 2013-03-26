# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

Prapare plots and tables for the BsMuMu AN

# <rawcell>

# Collection of scripts to produce the plots and tables for the BsMuMu AN.
# Use the runTMVA.sh script first to produce all the necessary files. To run the script do: "source runTMVA.sh" (It will submit 6 BDT trainings in parallel).
# It is better to first compile the macro TMVAClassification.C with "root -l"; ".L TMVAClassification.C+"; ".q".
# The first input is the cd into the directory where all the files are located.

# <headingcell level=2>

# Produce TMVA plots and copy them to latex figures dir

# <codecell>

cd /Users/demattia/TMVA-v4.1.2/test/NewTrees/

# <codecell>

figuresDir = "/Users/demattia/TMVA-v4.1.2/test/NewTrees/BsMuMuLatex/Figures/"

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
            # print fullName
            name = fullName.split("__")[0]
            hType = fullName.split("__")[1].split("_")[0]
            # print name, hType
            if hType == "Signal":
                canvas = TCanvas(name+"_canvas")
                stack = THStack(name+"_stack", name)
                canvas.Draw()
                histo.Scale(1/histo.Integral())
                stack.Add(histo)
                # histo.Draw()
                histo.SetLineColor(4)
                histo.SetFillColor(4);
                histo.SetFillStyle(3004);
                leg = TLegend(0.5536913,0.770979,0.9345638,0.9702797,"","brNDC")
                leg.AddEntry(histo, hType, "f")
            else:
                histo.Scale(1/histo.Integral())
                stack.Add(histo)
                # histo.Draw("same")
                stack.Draw("nostack")
                stack.SetMaximum(stack.GetMaximum("nostack")*1.35)
                stack.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
                stack.GetYaxis().SetTitle(histo.GetYaxis().GetTitle())
                # canvas.Update()
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

tablesDir = "/Users/demattia/TMVA-v4.1.2/test/NewTrees/BsMuMuLatex/Tables/"

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
                # else:
                #    print line
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
    # rank = splittedLine[3]
    variable = splittedLine[5]
    separation = splittedLine[7]
    # return rank+" & "+variable+" & "+separation
    return variable+" & "+separation
    
def prepareMultiTable(region, rankingType):
    outputFileBase = open(tablesDir+"ranking_barrel.txt", "w")

    outputFileBase.write("\\begin{table} \n")
    outputFileBase.write("  \\begin{tabular}{|l|l|l|l|l|l|l|} \n")
    outputFileBase.write("    \\hline \n")
    outputFileBase.write("    rank & variable & separation \\\\ \n")
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
    outputFileBase.write("  \\label{tab:datasets} \n")
    outputFileBase.write("\\end{table} \n")
    outputFileBase.close()

prepareMultiTable("barrel", "IdTransformation")
prepareMultiTable("endcaps", "IdTransformation")

# <codecell>


