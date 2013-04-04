#include "TString.h"
#include "TH1.h"
#include "TH2.h"
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
#include "setTDRStyle_modified.C"

const double soverb = 63./40361.;
//Barrel: B=40361 S=63
//Endcap: B=70362 S=40


void significance(TString method="BDT") {
  
  //  setTDRStyle(false);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillColor(0);

  TString fnameA("TMVA_barrel.root");
  TFile* inputA = TFile::Open(fnameA);

  gDirectory->Cd(fnameA+":/Method_"+method+"/"+method);
  gDirectory->pwd();
  TH1F* tmva_s = (TH1F*)gROOT->FindObject("MVA_BDT_S_high");
  TH1F* tmva_b = (TH1F*)gROOT->FindObject("MVA_BDT_B_high");

  const int rb =10;
  tmva_s->Rebin(rb);
  tmva_b->Rebin(rb);


  int nbins   = tmva_s -> GetNbinsX();
  double xmin = tmva_s -> GetBinLowEdge(1);
  double xmax = tmva_s -> GetBinLowEdge(nbins+1);

  double normS = tmva_s -> Integral(1,nbins);
  double sumS;//  = tmva_s -> Integral(tmva_s->FindBin(0),nbins);
  double normB = tmva_b -> Integral(1,nbins);
  double sumB;//  = tmva_b -> Integral(tmva_b->FindBin(0),nbins);
  double signi; 

  double sign_scale = normS/normB/soverb;

  printf("check bins: %d  %f %f %f %f \n", nbins, xmin, xmax, sumS, normS);

  TH1F* effS = new TH1F("effS","effS",nbins,xmin,xmax);
  TH1F* effB = new TH1F("effB","effB",nbins,xmin,xmax);
  TH1F* rejB = new TH1F("rejB","rejB",nbins,xmin,xmax);
  TH1F* sign = new TH1F("sign","sign",nbins,xmin,xmax);

  TH2F* roc  = new TH2F("roc","roc",nbins,0,1.02,nbins,0,1.02);

  for(int i=0; i<nbins; i++) {

    sumS  = tmva_s -> Integral(i+1,nbins);
    sumB  = tmva_b -> Integral(i+1,nbins);

    rejB->SetBinContent(i, 1.-sumB);
    effS->SetBinContent(i, normS?sumS/normS:0);
    effB->SetBinContent(i, normB?sumB/normB:0);
    rejB->SetBinContent(i, normB?1.-sumB/normB:1);

    roc->Fill(normS?sumS/normS:0,normB?1.-sumB/normB:1);

    signi = sumS/sqrt(sumS + sign_scale*sumB);
    sign->SetBinContent(i, signi);

  }


  TCanvas c;
  effS->SetTitle("");
  effS->GetYaxis()->SetTitle("#epsilon_{S} ,  1-#epsilon_{B}  ");
  effS->GetXaxis()->SetTitle(method+ " >     ");
  effS->Draw();
  rejB->Draw("same");
  c.Update();  

  Float_t rightmax = 1.02*sign->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  sign->SetLineColor(8);
  sign->Scale(scale);
  sign->SetLineWidth(3);

  sign->Draw("same");

  int sig_max_bin = sign->GetMaximumBin(); 
  double sig_max_mva = effS->GetBinCenter(sig_max_bin);
  double sig_max = sign->GetBinContent(sig_max_bin);
  double sig_max_effS = effS->GetBinContent(sig_max_bin);
  double sig_max_effB = effB->GetBinContent(sig_max_bin);
  // printf("maximum significance: %f bin:%d mva:%f %f %f \n", sig_max, sig_max_bin, sig_max_mva, sig_max_effS, sig_max_effB);
 char ss[100];
 sprintf(ss, "at significance maximum: %s=%4.3f, #epsilon_{S}=%4.3f, #epsilon_{B}=%5.4f",method.Data(), sig_max_mva, sig_max_effS, sig_max_effB);

  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
  	  gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor (8);
   axis->SetLabelColor(8);
   axis->SetTextColor (8);
   axis->SetTitle("Significance");
   axis->Draw();

   TPaveText *tp = new TPaveText(0.3,0.9,0.9,0.95,"brNDC");
   tp->SetBorderSize(0);
   tp->AddText(ss);
   tp->Draw();
   
   TLine* tl = new TLine(sig_max_mva,0,sig_max_mva, gPad->GetUymax());
   tl->SetLineStyle(3);
   tl->Draw("same");

   c.SaveAs("eff.gif");

   TCanvas c2;
   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("#epsilon_{S}");
   roc->GetYaxis()->SetTitle("1-#epsilon_{B}");
   roc->Draw();

   TMarker* tm = new TMarker(sig_max_effS,1.-sig_max_effB, 28);
   tm->Draw("same");

   c2.SaveAs("roc.gif");

   inputA->Close();
}
