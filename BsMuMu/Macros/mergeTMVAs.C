#include "setdirs.h"

TString fileName(const TString & method, const TString & region, const TString & index);

// Merge the signal and background entries histograms
void mergeTMVAs(TString method="BDT",TString region="barrel") {

  TH1F * tmva_s[3];
  TH1F * tmva_b[3];
  TString dir("Method_"+method+"/"+method);

  for( int i=0; i<3; ++i ) {
    //std::stringstream ss;ss << i; TString tss = ss.str(); //this code does not compile
    TString fname = fileName(method, region, TString::Format("%d",i));
    TFile* input = TFile::Open(fname);
    //gDirectory->Cd(fname+":/"+dir);
    input->cd(dir);
    //gDirectory->pwd();gDirectory->ls();
    tmva_s[i] = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_S_high"));
    tmva_b[i] = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_B_high"));
  }
  TFile* outputFile = new TFile(rootDir+"TMVA_"+region+"_"+method+"_merged.root", "RECREATE");
  TH1F* merged_s = (TH1F*) tmva_s[0]->Clone();
  TH1F* merged_b = (TH1F*) tmva_b[0]->Clone();
  merged_s->Add(tmva_s[1]);
  merged_s->Add(tmva_s[2]);
  merged_b->Add(tmva_b[1]);
  merged_b->Add(tmva_b[2]);
  outputFile->mkdir(dir);
  outputFile->cd(dir);
  merged_s->Write();
  merged_b->Write();
  outputFile->Close();
  std::cout << "created " << outputFile->GetName() << std::endl;

}

TString fileName(const TString & method, const TString & region, const TString & index) {
  TString fnameA = "TMVA_" + region;
  if(index=="merged") 
    fnameA += "_"+method;
  if(index!="") 
    fnameA += "_"+index;
  fnameA += ".root";
  return rootDir+fnameA;
}
