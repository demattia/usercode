from ROOT import TROOT, TLine, TGraphErrors, TLatex
from ROOT import TFile, TH1F, TH2F, TCanvas, THStack, TLegend
from ROOT import gBenchmark, gStyle, gROOT
from array import array


import subprocess
import sys
import os

## select mva methods
regions = ["barrel", "endcaps"]
methods = ["BDT", "MLP", "CutsSA"]
#regions = ["barrel"]
#methods = ["BDT"]
unblind = True

figuresDir  = "BsMuMuLatex/Figures/"
figuresDira = "BsMuMuLatex/Figures/bdt/"
figuresDirb = "BsMuMuLatex/Figures/mlp/"
figuresDirc = "BsMuMuLatex/Figures/cnt/"
figuresDird = "BsMuMuLatex/Figures/mainSel/"
tablesDir = "BsMuMuLatex/Tables/"
plotsDir = "plots/"
logsDir = "logs/"  # note: this is set in TMVAClassification.C
countersDir = "Trees/BsMC/"
weightsDir = "weights/"
rootDir = "rootfiles/"
#rootDir = "rootfiles_ams/"
rootExecutable = "root" #check if this is needed (currently not enforced) for running in specific environments
dirList = ["figuresDir","tablesDir","plotsDir","logsDir","countersDir","weightsDir","rootDir","figuresDira","figuresDirb","figuresDirc","figuresDird"]

# check otherwise create root directories
def checkdirs():
    if not os.path.exists("Trees"):
        print "no input trees found, exiting..."
        exit(3)
    for dir in dirList: #needed are: rootDir,logsDir,plotsDir,weightsDir,figuresDir
        if not os.path.exists(dir):
            print "creating directory:", dir, ":", globals()[dir] 
            os.system("mkdir -p "+globals()[dir])
checkdirs()

# generate c++ header
def genheader():
    hf = "setdirs.h"
    os.system("cp _includes.h "+hf)
    with open(hf,'r+') as f:
        content = f.read()
        #print content
        for dir in dirList: #needed are: rootDir,logsDir,plotsDir,weightsDir,figuresDir
            if dir not in content: 
                f.write("TString "+dir+"(\""+globals()[dir]+"\");\n")
genheader()

# list of input trees
inputTrees = [
    'Trees/Run2012D/selection_test.root',
    'Trees/Run2012A/selection_test.root',
    'Trees/Run2012ARecover/selection_test.root',
    'Trees/Run2012B/selection_test.root',
    'Trees/Run2012C1/selection_test.root',
    'Trees/Run2012DRereco/selection_test.root',
    'Trees/Run2012C2/selection_test.root',
#    'trees-run2012dcleaned/selection_test.root',
    'Trees/BsMC/selection_test.root',
              ]
#runlistA = ["Run2012A","Run2012ARecover","Run2012B","Run2012C1","Run2012DRereco","Run2012C2","Run2012D"]
runlist = []
for idx, tree in enumerate(inputTrees):
    if not "MC" in tree:
        runlist.append(tree.split("/")[-2])
#print "runlist:", runlist

# maxRun = "203002"  # this is the maximum run number (included) considered for data
maxRun = "999999"  # this is the maximum run number (included) considered for data
