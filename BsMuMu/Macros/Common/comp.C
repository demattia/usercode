#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TROOT.h>

#include "setTDRStyle_modified.C"

//using namespace std; 

void comp() {

  TString fnameA("data_preselection_endcaps.root");
  //TString fnameA("data_barrel_xcheck.root");
  TFile* inputA = TFile::Open( fnameA );
  
  TString fnameB ("../compare/data_endcaps_main.root");
  //TString fnameB ("../compare/data_barrel_main.root");
  TFile* inputB = TFile::Open( fnameB );

  gDirectory->Cd(fnameA+":/");
  TTree* treeA = (TTree*)gROOT->FindObject("probe_tree");
  //treeA->Show();

  gDirectory->Cd(fnameB+":/");
  TTree* treeB = (TTree*)gROOT->FindObject("events");
  //treeB->Show(); 

  TFile outfile("comp.root","recreate");
  //gDirectory->Cd(*outfile+":/");

  gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();


  const int nvar = 11;
  TString varsA[nvar]= {"mass",  "pt", "dca",      "mu1_pt","mu2_pt", "delta3d","delta3d/delta3dErr","cosAlpha",       "ntrk",  "minDca", "isolation"};  //"l3d",   "l3dSig",   "alpha"
  TString varsB[nvar]= {   "m",  "pt", "maxdoca",    "m1pt",  "m2pt",    "pvip",   "pvips",              "cosa",   "closetrk", "docatrk", "iso"      };  //"fl3d",   "fls3d",   "alpha"

  Int_t  nbins[nvar] = {    40,   40,   40,            40,     40,          40,        40,                   40,           20,        40,       40 };   // 40,   40,     40,    
  Double_t min[nvar] = {   4.9,    0,    0,             0,      0,           0,         0,                 0.95,            0,         0,        0 };   //  0,    0,      0,    
  Double_t max[nvar] = {   5.9,   50,  0.1,            50,     50,        0.05,         5,                  1.0,           20,       0.2,        1 };   //100,   10,   0.03,   

  TString dir = "figs/";
  TString ext = ".gif";
  char tmp[100];
  TString varA(""), varB(""), cvsn("");
    
  for(int i=0; i<nvar; i++) {
    std::cout << "variable: " << varsA[i] << " -> " << varsB[i] << std::endl; 
    varA=""; varB=""; cvsn="";
    cvsn.Append(varsA[i].Data());
    varA.Append(varsA[i].Data());
    varB.Append(varsB[i].Data());
    sprintf(tmp,"(%d,%4.2f,%4.2f)",nbins[i],min[i],max[i]);
    TString lim (""); lim.Append(tmp);
    cvsn = dir + varB + ext;
    TCanvas cvs(cvsn);
    varA+=" >> hA";    
    varB+=" >> hB";
    varA+=lim; 
    varB+=lim;
    std::cout << "plot: " <<varA << " " << varB << std::endl << std::flush;
    treeA->Draw(varA.Data());
    treeB->Draw(varB.Data());
    TH1F *hA = (TH1F*)gDirectory->Get("hA");
    TH1F *hB = (TH1F*)gDirectory->Get("hB");
    if(hA)    
      hA->Scale(1./hA->Integral());
    if(hB)
      hB->Scale(1./hB->Integral());
    hA->SetLineColor(kRed);
    hB->SetLineColor(kBlue);
    hA->SetLineWidth(2);
    hB->SetLineWidth(2);

    TLegend* leg = new TLegend(0.5,0.7,0.7,0.8);
    leg->AddEntry(hA,"xcheck","pl");
    leg->AddEntry(hB,"main","pl");

    cvs.cd(); 
    THStack hs(varsA[i].Data(),varsA[i].Data());
    //THStack hs("hs",hA->GetTitle());
    hs.Add(hA);
    hs.Add(hB);
    hs.Draw("nostack");

    leg->Draw("same");
    cvs.SaveAs(cvsn);
  
    hs.Write();
    cvs.Write();
  }
  

  inputA->Close();
  inputB->Close();
  outfile.Close();

  return;
}
