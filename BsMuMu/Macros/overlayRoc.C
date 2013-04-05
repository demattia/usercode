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
#include "../Common/setTDRStyle_modified.C"

TH2F* getplot(TString filename);
TGraph* getpoint(TString filename);

void overlayRoc(){

  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillColor(0);


  TCanvas c;

  TH2F* roc =getplot("plots/signif_BDT_barrel.root");
  TH2F* roc0=getplot("plots/signif_BDT_barrel_0.root");
  TH2F* roc1=getplot("plots/signif_BDT_barrel_1.root");
  TH2F* roc2=getplot("plots/signif_BDT_barrel_2.root");
  roc->GetXaxis()->SetRangeUser(0.,1);
  roc->GetYaxis()->SetRangeUser(0.99,1.);

  roc ->SetMarkerColor(1);
  roc0->SetMarkerColor(2);
  roc1->SetMarkerColor(3);
  roc2->SetMarkerColor(4);  

  roc ->SetMarkerStyle(24);
  roc0->SetMarkerStyle(25);
  roc1->SetMarkerStyle(26);
  roc2->SetMarkerStyle(27);

  roc ->SetMarkerSize(0.8);
  roc0->SetMarkerSize(0.4);
  roc1->SetMarkerSize(0.4);
  roc2->SetMarkerSize(0.4);

  roc ->Draw();
  roc0->Draw("same");
  roc1->Draw("same");
  roc2->Draw("same");

  /*
  double x, y; 
  x = getpoint("plots/signif_BDT_barrel.root")->GetX()[0];
  y = getpoint("plots/signif_BDT_barrel.root")->GetY()[0];
  printf("cccccc %f   %f\n ",x,y);
  TMarker* tm = new TMarker(x,y, 28);
  tm->SetMarkerSize(1.5);
  tm->SetMarkerColor(2);
  tm->Draw("same");
  */

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  //leg->SetHeader("The Legend Title");
  leg->AddEntry(roc,"combined","p");
  leg->AddEntry(roc0,"0","p");
  leg->AddEntry(roc1,"1","p");
  leg->AddEntry(roc2,"2","p");
  leg->Draw();
  
  c.SaveAs("plots/roc_overlay_barrel.pdf");
}


TH2F* getplot(TString filename){
  TFile::Open(filename);
  TH2F* h = (TH2F*)gROOT->FindObject(TString("roc"));
  return h;
}

TGraph* getpoint(TString filename){
  TFile::Open(filename);
  return (TGraph*)gROOT->FindObject(TString("mark"));
}
