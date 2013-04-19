#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TText.h"
#include <iostream>
#include "setTDRStyle_modified.C"

struct Vars
{
  Vars(const TString & inVarA, const TString & inVarB, const Int_t inNbins, const Double_t & inMin, const Double_t & inMax)
    : varA(inVarA), varB(inVarB), nbins(inNbins), min(inMin), max(inMax) {}
  // varA is the name in the parallel analysis tree, while varB in the reference tree from the main analysis
  TString varA;
  TString varB;
  Int_t nbins;
  Double_t min;
  Double_t max;
};


void comp(const TString & parallelAnalysisFile = "Barrel_preselection.root",
	  const TString & mainAnalysisFile = "data_main_barrel.root")
{
  TString fnameA(parallelAnalysisFile);
  TFile* inputA = TFile::Open( fnameA );
  
  TString fnameB (mainAnalysisFile);
  TFile* inputB = TFile::Open( fnameB );

  gDirectory->Cd(fnameA+":/");
  TTree* treeA = (TTree*)gROOT->FindObject("probe_tree");
  //treeB->Show();

  gDirectory->Cd(fnameB+":/");
  TTree* treeB = (TTree*)gROOT->FindObject("events");
  //treeA->Show();

  TFile outfile("comp.root","recreate");
  //gDirectory->Cd(*outfile+":/");

  // gROOT->LoadMacro("setTDRStyle_modified.C");
  setTDRStyle();


  std::vector<Vars> vars;
  vars.push_back(Vars("mass", "m", 40, 4.9, 5.9));
  vars.push_back(Vars("pt", "pt", 40, 0, 50));
  vars.push_back(Vars("dca", "maxdoca", 40, 0, 0.1));
  vars.push_back(Vars("mu1_pt", "m1pt", 40, 0, 50));
  vars.push_back(Vars("mu2_pt", "m2pt", 40, 0, 50));
  vars.push_back(Vars("delta3d", "pvip", 40, 0, 0.05));
  vars.push_back(Vars("delta3d/delta3dErr", "pvips", 40, 0, 5));
  vars.push_back(Vars("cosAlpha3D", "cosa", 40, 0.95, 1.));
  vars.push_back(Vars("ntrk20", "closetrk", 40, 0, 20));
  vars.push_back(Vars("minDca", "docatrk", 40, 0, 0.2));
  vars.push_back(Vars("isolation", "iso", 40, 0, 1));
  vars.push_back(Vars("l3d", "fl3d", 40, 0, 5));
  vars.push_back(Vars("l3dsig", "fls3d", 40, 0, 10));
  vars.push_back(Vars("pvw8", "pvw8", 40, 0, 1.1));
  // vars.push_back(Vars("alpha3D", "alpha"));

  TString dir = "figs/";
  TString ext = ".pdf";
  char tmp[100];

  TString varA(""), varB(""), cvsn("");
  // TCanvas cvs[nvar];

  Double_t ks = 0;
  Int_t i=0;
  TH1F *KTest = new TH1F("KTest","KTest",50,0,1);

  // for(int i=0; i<nvar; i++) {
  for( std::vector<Vars>::const_iterator it = vars.begin(); it != vars.end(); ++it ) {
    std::cout << "variable: " << it->varA << " -> " << it->varB << std::endl; // << flush;
    varA=""; varB=""; cvsn="";
    cvsn.Append(it->varA.Data());
    varA.Append(it->varA.Data());
    varB.Append(it->varB.Data());
    sprintf(tmp,"(%d,%4.2f,%4.2f)",it->nbins,it->min,it->max);
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

    hA->Sumw2();
    hB->Sumw2();
    if(hA)    
      //hA->Scale(1.);
      hA->Scale(1./hA->Integral());
    if(hB)
      //hB->Scale(1. * hA->Integral()/hB->Integral());
      hB->Scale(1./hB->Integral());
    hA->SetLineColor(kRed);
    hB->SetLineColor(kBlue);

    hA->SetLineWidth(2);
    hB->SetLineWidth(2);
    TLegend leg(0.5,0.7,0.7,0.8);
    leg.AddEntry(hA,"xcheck","pl");
    leg.AddEntry(hB,"main","pl");

    cvs.cd(); 
    THStack hs("hs",hA->GetTitle());
    hs.Add(hA);
    hs.Add(hB);
    hs.Draw("HISTnostack");

	cout<<hA->GetBinContent(5)<<endl;
	cout<<hB->GetBinContent(5)<<endl;
    ks=hA->KolmogorovTest(hB,"D");
    cout<<"KS = "<<ks<<endl;
    TString probatext = Form( "P(KS) = %5.3g", ks );
    TText* tt = new TText( 0.55, 0.65, probatext );
    tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad();
    KTest->Fill(ks);

    leg.Draw("same");
    cvs.SaveAs(cvsn);
  }

  TCanvas Kcvs(dir+"KS");
  Kcvs.cd();
  KTest->SetTitle("KolmogorovTest");
  KTest->SetFillColor(kYellow);
  KTest->Draw("");
  Kcvs.SaveAs(dir+"KS"+ext);

  return;
}


///// Urs variables
/*

run             = 0
 json            = 0
 evt             = 0
 ls              = 0
 tm              = 0
 pr              = 0
 procid          = 0
 hlt             = 0
 pvn             = 0
 cb              = 0
 rr              = 0
 bdt             = 0
 npv             = 0
 pvw8            = 0
 gmuid           = 0
 gmupt           = 0
 gmueta          = 0
 gtqual          = 0
 gtpt            = 0
 gteta           = 0
 w8mu            = 0
 w8tr            = 0
 pvlip           = 0
 pvlips          = 0
 pvlip2          = 0
 pvlips2         = 0
 pvip            = 0
 pvips           = 0
 q               = 0
 type            = 0
 pt              = 0
 eta             = 0
 phi             = 0
 tau             = 0
 m               = 0
 cm              = 0
 cosa            = 0
 alpha           = 0
 iso             = 0
 isotrk          = 0
 closetrk        = 0
 chi2            = 0
 dof             = 0
 prob            = 0
 fls3d           = 0
 fl3d            = 0
 flxy            = 0
 fl3dE           = 0
 flsxy           = 0
 docatrk         = 0
 docatrkbdt      = 0
 maxdoca         = 0
 lip             = 0
 lipE            = 0
 tip             = 0
 tipE            = 0
 osiso           = 0
 osreliso        = 0
 osmpt           = 0
 osmptrel        = 0
 osmdr           = 0
 m1q             = 0
 m1id            = 0
 m1pt            = 0
 m1eta           = 0
 m1phi           = 0
 m1ip            = 0
 m1gt            = 0
 m1pix           = 0
 m1bpix          = 0
 m1bpixl1        = 0
 m1chi2          = 0
 m1pv            = 0
 m2q             = 0
 m2id            = 0
 m2pt            = 0
 m2eta           = 0
 m2phi           = 0
 m2ip            = 0
 m2gt            = 0
 m2pix           = 0
 m2bpix          = 0
 m2bpixl1        = 0
 m2chi2          = 0
 m2pv            = 0
 mudist          = 0
 mudeltar        = 0
 g1pt            = 0
 g2pt            = 0
 g1eta           = 0
 g2eta           = 0
 g1phi           = 0
 g2phi           = 0
 gmass           = 0
 gtau            = 0
 t1pt            = 0
 t1eta           = 0
 t2pt            = 0
 t2eta           = 0
 hm1pt           = 0
 hm1eta          = 0
 hm1phi          = 0
 hm2pt           = 0
 hm2eta          = 0
 hm2phi          = 0


OURS  -------------------------

 NChi2           = 0
 PDGid           = 0
 charge          = 0
 cosAlpha        = 0
 countTksOfPV    = 0
 ctauErrPV       = 0
 ctauPV          = 0
 dca             = 0
 dcaxy           = 0
 delta3d         = 0
 delta3dErr      = 0
 isolation       = 0
 mass            = 0
 minDca          = 0
 mu1_R03emEt     = 0
 mu1_R03hadEt    = 0
 mu1_R03sumPt    = 0
 mu1_R05emEt     = 0
 mu1_R05hadEt    = 0
 mu1_R05sumPt    = 0
 mu1_charge      = 0
 mu1_dB          = 0
 mu1_dxy         = 0
 mu1_dz          = 0
 mu1_eta         = 0
 mu1_nChi2       = 0
 mu1_nMuSegs     = 0
 mu1_nMuSegsCln  = 0
 mu1_nPixHits    = 0
 mu1_nTrHits     = 0
 mu1_phi         = 0
 mu1_pt          = 0
 mu2_R03emEt     = 0
 mu2_R03hadEt    = 0
 mu2_R03sumPt    = 0
 mu2_R05emEt     = 0
 mu2_R05hadEt    = 0
 mu2_R05sumPt    = 0
 mu2_charge      = 0
 mu2_dB          = 0
 mu2_dxy         = 0
 mu2_dz          = 0
 mu2_eta         = 0
 mu2_nChi2       = 0
 mu2_nMuSegs     = 0
 mu2_nMuSegsCln  = 0
 mu2_nPixHits    = 0
 mu2_nTrHits     = 0
 mu2_phi         = 0
 mu2_pt          = 0
 ntrk20            = 0
 pt              = 0
 sumPTPV         = 0
 vProb           = 0
 vertexWeight    = 0
 y               = 0
 mu1_Calo        = 0
 mu1_GM          = 0
 mu1_GMPT        = 0
 mu1_HLT_DoubleMu2BarrelBsL3 = 0
 mu1_HLT_DoubleMu2BsL3 = 0
 mu1_HLT_DoubleMu2Dimuon6BsL3 = 0
 mu1_HLT_DoubleMu3BsL3 = 0
 mu1_HLT_VertexmumuFilterBs345 = 0
 mu1_HLT_VertexmumuFilterBs3p545 = 0
 mu1_HLT_VertexmumuFilterBs4 = 0
 mu1_HLT_VertexmumuFilterBs47 = 0
 mu1_HLT_VertexmumuFilterBs6 = 0
 mu1_TM          = 0
 mu1_TMLSAT      = 0
 mu1_TMOSAT      = 0
 mu1_TMOST       = 0
 mu2_Calo        = 0
 mu2_GM          = 0
 mu2_GMPT        = 0
 mu2_HLT_DoubleMu2BarrelBsL3 = 0
 mu2_HLT_DoubleMu2BsL3 = 0
 mu2_HLT_DoubleMu2Dimuon6BsL3 = 0
 mu2_HLT_DoubleMu3BsL3 = 0
 mu2_HLT_VertexmumuFilterBs345 = 0
 mu2_HLT_VertexmumuFilterBs3p545 = 0
 mu2_HLT_VertexmumuFilterBs4 = 0
 mu2_HLT_VertexmumuFilterBs47 = 0
 mu2_HLT_VertexmumuFilterBs6 = 0
 mu2_TM          = 0
 mu2_TMLSAT      = 0
 mu2_TMOSAT      = 0
 mu2_TMOST       = 0
 run             = 0
 lumi            = 0
 event           = 0


*/
