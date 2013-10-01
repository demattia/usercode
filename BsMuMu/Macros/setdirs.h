#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStopwatch.h"

#include "Common/setTDRStyle_modified.C"

TString figuresDir("BsMuMuLatex/Figures/");
TString tablesDir("BsMuMuLatex/Tables/");
TString plotsDir("plots/");
TString logsDir("logs/");
TString countersDir("Trees/BsMC/");
TString weightsDir("weights/");
TString rootDir("rootfiles/");
TString figuresDira("BsMuMuLatex/Figures/bdt/");
TString figuresDirb("BsMuMuLatex/Figures/mlp/");
TString figuresDirc("BsMuMuLatex/Figures/cnt/");
TString figuresDird("BsMuMuLatex/Figures/mainSel/");
typedef float mytype;