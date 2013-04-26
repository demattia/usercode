#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include <iostream>
#include <sstream>

TString fileName(const TString & method, const TString & region, const TString & index)
{
  TString fnameA = "TMVA_" + region;
  if(index!="") 
    fnameA += "_"+index;
  fnameA += ".root";
  return fnameA;
}

// Merge the signal and background entries histograms
void mergeTMVAs(TString method="BDT",TString region="barrel")
{
  TH1F * tmva_s[3];
  TH1F * tmva_b[3];
  for( int i=0; i<3; ++i ) {
    std::stringstream ss;
    ss << i;
    TString fname = fileName(method, region, ss.str());
    TFile* input = TFile::Open(fname);

    gDirectory->Cd(fname+":/Method_"+method+"/"+method);
    gDirectory->pwd();
    tmva_s[i] = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_S_high"));
    tmva_b[i] = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_B_high"));
  }
  TFile * outputFile = new TFile("TMVA_"+region+"_merged.root", "RECREATE");
  TH1F * merged_s = tmva_s[0]->Clone();
  TH1F * merged_b = tmva_b[0]->Clone();
  merged_s->Add(tmva_s[1]);
  merged_s->Add(tmva_s[2]);
  merged_b->Add(tmva_b[1]);
  merged_b->Add(tmva_b[2]);
  outputFile->cd();
  merged_s->Write();
  merged_b->Write();
  outputFile->Close();
}
