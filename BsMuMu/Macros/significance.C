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
#include "Common/setTDRStyle_modified.C"
#include <fstream>

const double nsig_barrel  = 60.;
const double nbkg_barrel  = 28884.;
const double nsig_endcap  = 35.;
const double nbkg_endcap  = 35392.;

void significance(TString method="BDT",TString region="barrel", TString index="", const int subdir = 1) {
  
  double nsig, nbkg;
  if(region=="barrel") {
    nsig = nsig_barrel;
    nbkg = nbkg_barrel;
  } else if (region=="endcaps") {
    nsig = nsig_endcap;
    nbkg = nbkg_endcap;
  } else  {cout<<"wrong input"<<endl; exit(-1);}

  cout << "processing " << method << " for " << region << endl;

  TString name = method + "_" + region;
  if(index!="") 
    name += "_"+index;

  TString fnameA = "TMVA_" + region;
  if(index!="") 
    fnameA += "_"+index;
  fnameA += ".root";
  TFile* inputA = TFile::Open(fnameA);

  if( subdir ) gDirectory->Cd(fnameA+":/Method_"+method+"/"+method);
  gDirectory->pwd();
  TH1F* tmva_s = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_S_high"));
  TH1F* tmva_b = (TH1F*)gROOT->FindObject(TString("MVA_"+method+"_B_high"));

  const int rb = 250;
  // const int rb = 1;
  tmva_s->Rebin(rb);
  tmva_b->Rebin(rb);

  int nbins   = tmva_s -> GetNbinsX();
  double xmin = tmva_s -> GetBinLowEdge(1);
  double xmax = tmva_s -> GetBinLowEdge(nbins+1);

  double normS = tmva_s -> Integral(1,nbins);
  double normB = tmva_b -> Integral(1,nbins);
  tmva_s->Scale(nsig/normS);
  tmva_b->Scale(nbkg/normB);
  normS = tmva_s -> Integral(1,nbins);
  normB = tmva_b -> Integral(1,nbins);

  double sumS(0), sumB(0), signi(0), signi1(0), signi2(0); 

  printf("check bins: %d  %f %f %f %f \n", nbins, xmin, xmax, sumS, normS);

  TH1F* effS = new TH1F("effS","effS",nbins,xmin,xmax);
  TH1F* effB = new TH1F("effB","effB",nbins,xmin,xmax);
  TH1F* rejB = new TH1F("rejB","rejB",nbins,xmin,xmax);
  TH1F* sign = new TH1F("sign","sign",nbins,xmin,xmax);
  TH1F* sign1 = new TH1F("sign1","sign1",nbins,xmin,xmax);
  TH1F* sign2 = new TH1F("sign2","sign2",nbins,xmin,xmax);

  TH2F* roc  = new TH2F("roc","roc",nbins,0,1.02,nbins,0,1.02);

  for(int i=0; i<nbins+1; i++) {

    sumS  = tmva_s -> Integral(i,nbins);
    sumB  = tmva_b -> Integral(i,nbins);

    rejB->SetBinContent(i, 1.-sumB);
    effS->SetBinContent(i, normS?sumS/normS:0);
    effB->SetBinContent(i, normB?sumB/normB:0);
    rejB->SetBinContent(i, normB?1.-sumB/normB:1);

    roc->Fill(normS?sumS/normS:0,normB?1.-sumB/normB:1);

    signi = sumS/sqrt(sumS + sumB);
    sign->SetBinContent(i, signi);

    signi1 = sumB?sumS/sqrt(sumB):0;
    sign1->SetBinContent(i, signi1);

    signi2 = sumS/(sqrt(sumB)+0.5);
    sign2->SetBinContent(i, signi2);

    printf("--> signif : %f  %f %f\n", 
	   signi, signi1, signi2 );

  }

  int sig_max_bin = sign->GetMaximumBin(); 
  double sig_max_mva = effS->GetBinCenter(sig_max_bin);
  double sig_max = sign->GetBinContent(sig_max_bin);
  double sig_max_effS = effS->GetBinContent(sig_max_bin);
  double sig_max_effB = effB->GetBinContent(sig_max_bin);
  //printf("maximum significance: %f bin:%d mva:%f %f %f \n", sig_max, sig_max_bin, sig_max_mva, sig_max_effS, sig_max_effB);
  char ss[100];
  sprintf(ss, "%s,  max. significance: %3.1f, %s>%5.4f, #epsilon_{S}=%4.3f, #epsilon_{B}=%5.4f",region.Data(),sig_max,method.Data(), sig_max_mva, sig_max_effS, sig_max_effB);

  ofstream outputTxt("maxsignificance_"+region+".txt");
  outputTxt << ss;
  outputTxt.close();

  //setTDRStyle(false);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(0); 

  TCanvas c;
  c.SetGrid();
  effS->SetTitle("");
  effS->GetYaxis()->SetTitle("#epsilon_{S}   |   1-#epsilon_{B}  ");
  effS->GetXaxis()->SetTitle(method+ " >     ");
  effS->SetLineColor(kBlue);
  effS->Draw("c");
  rejB->SetLineColor(kRed);
  rejB->Draw("same c");
  c.Update();  

  Float_t rightmax = 1.02*sign1->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;

  sign->Scale(scale);
  sign1->Scale(scale);
  sign2->Scale(scale);

  sign ->SetLineColor(8);
  sign1->SetLineColor(38);
  sign2->SetLineColor(42);

  sign ->SetLineWidth(2);
  sign1->SetLineWidth(2);
  sign2->SetLineWidth(2);

  sign1->SetLineStyle(3);
  sign2->SetLineStyle(2);

  sign ->Draw("same c");
  sign1->Draw("same c");
  sign2->Draw("same c");

  TLegend* leg = new TLegend(0.12,0.65,0.45,0.75);
  //leg->SetHeader("The Legend Title");
  char ss0[50], ss1[50],ss2[50];
  sprintf(ss0, "S/#sqrt{S+B},     max.@ %s>%5.4f", method.Data(), effS->GetBinCenter(sign ->GetMaximumBin()));
  sprintf(ss1, "S/#sqrt{B},          max.@ %s>%5.4f", method.Data(), effS->GetBinCenter(sign1->GetMaximumBin()));
  sprintf(ss2, "S/(#sqrt{B}+0.5), max.@ %s>%5.4f", method.Data(), effS->GetBinCenter(sign2->GetMaximumBin()));
  leg->AddEntry(sign ,ss0,"l");
  leg->AddEntry(sign1,ss1,"l");
  leg->AddEntry(sign2,ss2,"l");
  leg->Draw();


  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
  	  gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor (8);
   axis->SetLabelColor(8);
   axis->SetTextColor (8);
   axis->SetTitle("Significance");
   axis->Draw();

   TPaveText *tp = new TPaveText(0.2,0.9,0.9,0.95,"brNDC");
   tp->SetBorderSize(0);
   tp->AddText(ss);
   tp->Draw("same");
   
   TLine* tl = new TLine(sig_max_mva,0,sig_max_mva, gPad->GetUymax());
   tl->SetLineStyle(4);
   tl->Draw("same");

   TString ext(".pdf");

   c.SaveAs(TString("plots/"+name+"_eff"+ext));

   TCanvas c2;
   roc->SetTitle("");
   roc->GetXaxis()->SetTitle("#epsilon_{S}");
   roc->GetYaxis()->SetTitle("1-#epsilon_{B}");
   //roc->SetMarkerStyle(7);
   roc->Draw();
   //tp->Draw();
   TMarker* tm = new TMarker(sig_max_effS,1.-sig_max_effB, 28);
   tm->SetMarkerSize(1.5);
   tm->SetMarkerColor(2);
   tm->Draw("same");
   c2.SaveAs(TString("plots/"+name+"_roc"+ext));

   TCanvas c3;
   roc->GetXaxis()->SetRangeUser(0.,1);
   roc->GetYaxis()->SetRangeUser(0.99,1.);
   roc->Draw();
   //tp->Draw("same");
   tm->Draw("same");
   c3.SaveAs(TString("plots/"+name+"_roc_zoom"+ext));

   TGraph* mark = new TGraph(1);
   mark->SetPoint(1,sig_max_effS,1.-sig_max_effB);
   mark->SetName("marker");
   TString outfileName("plots/signif_"+name+".root");
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   roc->Write();   
   mark->Write();   
   outputFile->Close();

   
   inputA->Close();
}
