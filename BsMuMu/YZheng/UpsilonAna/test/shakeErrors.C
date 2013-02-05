#include "TH1.h"
#include "TFile.h"
#include "TROOT.h"

void shakeErrors(TString inFileName, TString histName){
  TFile inFile(inFileName);
  TH2F* h = (TH2F*)gROOT->FindObject(histName);
  if(!h){
    std::cout<<"Could not access histogram!"<<std::endl;
    return;
  }

  TFile lowFile(TString("low_")+inFileName,"recreate");
  TH2F* hLow = (TH2F*)h->Clone(histName);
  for(int i=0; i<h->GetSize(); i++){
    hLow->SetBinContent( i, h->GetBinContent(i)-h->GetBinError(i) );
  }
  hLow->Write();
  lowFile.Close();

  TFile hiFile(TString("hi_")+inFileName,"recreate");
  TH2F* hHi = (TH2F*)h->Clone(histName);
  for(int i=0; i<h->GetSize(); i++){
    hHi->SetBinContent( i, h->GetBinContent(i)+h->GetBinError(i) );
  }
  hHi->Write();
  hiFile.Close();
}

