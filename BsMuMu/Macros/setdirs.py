from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import THStack
from ROOT import TROOT
from ROOT import TLine
from ROOT import TGraphErrors
from array import array

import subprocess
import sys
import os

## select mva methods
regions = ["barrel", "endcaps"]
methods = ["BDT", "MLP", "CutsSA"]
#regions = ["barrel"]
#methods = ["BDT"]


figuresDir = "BsMuMuLatex/Figures/"
tablesDir = "BsMuMuLatex/Tables/"
plotsDir = "plots/"
logsDir = "logs/"  # note: this is set in TMVAClassification.C
countersDir = "Trees/BsMC/"
weightsDir = "weights/"
rootDir = "rootfiles/"
rootExecutable = "root" #check if this is needed (currently not enforced) for running in specific environments

# check otherwise create root directories
def checkdir(dirname):
    if not os.path.exists(dirname):
        print "creating directory:", dirname
        os.system("mkdir "+dirname)
checkdir(rootDir)
checkdir(logsDir)
checkdir(weightsDir)

# generate c++ header
def genheader():
    hf = "setdirs.h"
    os.system("cp _includes.h "+hf)
    with open(hf,'r+') as f:
        content = f.read()
        print content
        for dir in ["figuresDir","tablesDir","logsDir","countersDir","weightsDir","rootDir"]: #needed are: rootDir,logsDir,plotsDir,weightsDir,figuresDir
            if dir not in content: 
                f.write("TString "+dir+"(\""+globals()[dir]+"\");\n")
genheader()

# list of input trees
inputTrees = ['Trees/Run2012A/selection_test.root',
              'Trees/Run2012ARecover/selection_test.root',
              'Trees/Run2012B/selection_test.root',
              'Trees/Run2012C1/selection_test.root',
              'Trees/Run2012DRereco/selection_test.root',
              'Trees/Run2012C2/selection_test.root',
              'Trees/Run2012D/selection_test.root',
              'Trees/BsMC/selection_test.root',
              ]
#runlistA = ["Run2012A","Run2012ARecover","Run2012B","Run2012C1","Run2012DRereco","Run2012C2","Run2012D"]
runlist = []
for idx, tree in enumerate(inputTrees):
    if not "MC" in tree:
        runlist.append(tree.split("/")[-2])
#print runlist

maxRun = "203002"  # this is the maximum run number (included) considered for data
