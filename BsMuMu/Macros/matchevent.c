#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <iostream>

struct Comp
{
  Comp(const TString & mainName, const TString & xCheckName, const unsigned int bins, const double & min, const double & max) :
    name(mainName), diff(0.), absDiff(0.), relDiff(0.), relAbsDiff(0.)
  {
    main = new TH1F(mainName, mainName, bins, min, max);
    xCheck = new TH1F(xCheckName, xCheckName, bins, min, max);
    uxCheck = new TH1F("unmatched_"+xCheckName, "unmatched "+xCheckName, bins, min, max);
  }
  void fill(const double & xCheckValue, const double & mainValue)
  {
    xCheck->Fill(xCheckValue);
    main->Fill(mainValue);
    // Compute differences
    diff += mainValue - xCheckValue;
    absDiff += fabs(mainValue - xCheckValue);
    if( mainValue != 0 ) {
      relDiff += (mainValue - xCheckValue)/mainValue;
      relAbsDiff += fabs(mainValue - xCheckValue)/mainValue;
    }
  }
  void fill(const double & xCheckValue)
  {
    uxCheck->Fill(xCheckValue);
  }
  void makeplot(TCanvas *c, TCanvas *unmatched_c,
                TH1F * histoDiff, TH1F * histoAbsDiff,
                TH1F * histoRelDiff, TH1F * histoRelAbsDiff,
                int n)
  {
    // Draw matched
    c->cd(n);
    main->SetLineColor(kRed);
    main->GetXaxis()->SetTitle(name);
    main->Draw();
    xCheck->SetLineColor(kBlue);
    xCheck->SetLineStyle(kDashed);
    xCheck->Draw("same");
    // Draw unmatched
    unmatched_c->cd(n);
    uxCheck->SetLineColor(kRed);
    uxCheck->GetXaxis()->SetTitle(name);
    uxCheck->Draw();
    // Differences
    histoDiff->SetBinContent(n, diff);
    histoAbsDiff->SetBinContent(n, absDiff);
    histoRelDiff->SetBinContent(n, relDiff);
    histoRelAbsDiff->SetBinContent(n, relAbsDiff);
  }

  TString name;
  TH1F * main;
  TH1F * xCheck;
  TH1F * uxCheck;
  double diff;
  double absDiff;
  double relDiff;
  double relAbsDiff;
};

void matchevent()
{
  unsigned int maxEvents = 5;
  TFile *fxchk = new TFile("BsMC12_preselection.root");
  // TFile *fxchk = new TFile("check_barrel_preselection.root");
  TFile *fmain = new TFile("MainAnalysis/MC_afterCuts.root");

  TTree *Txchk = (TTree*)fxchk->Get("probe_tree");
  TTree *Tmain = (TTree*)fmain->Get("events");

  Float_t xpt, xeta, mu1_pt, mu2_pt, mu1_iso, mu2_iso, l3d, l3dsig, mass, dca, delta3d, cosAlpha3D, minDca, isolation, NChi2;
  Float_t ntrk20;
  UInt_t event, xrun;
  Double_t pt, eta, m1pt, m2pt, m1iso, m2iso, fl3d, fls3d, m, maxdoca, pvip, cosa, docatrk, iso, chi2, dof;
  Int_t closetrk;
  Long64_t evt, run;

  Txchk->SetBranchAddress("pt",&xpt);
  Txchk->SetBranchAddress("eta",&xeta);
  Txchk->SetBranchAddress("mu1_pt",&mu1_pt);
  Txchk->SetBranchAddress("mu2_pt",&mu2_pt);
  Txchk->SetBranchAddress("mu1_iso",&mu1_iso);
  Txchk->SetBranchAddress("mu2_iso",&mu2_iso);
  Txchk->SetBranchAddress("l3d",&l3d);
  Txchk->SetBranchAddress("l3dsig",&l3dsig);
  Txchk->SetBranchAddress("mass",&mass);
  Txchk->SetBranchAddress("ntrk20",&ntrk20);
  Txchk->SetBranchAddress("dca",&dca);
  Txchk->SetBranchAddress("delta3d",&delta3d);
  Txchk->SetBranchAddress("cosAlpha3D",&cosAlpha3D);
  Txchk->SetBranchAddress("minDca",&minDca);
  Txchk->SetBranchAddress("isolation",&isolation);
  Txchk->SetBranchAddress("NChi2", &NChi2);
  Txchk->SetBranchAddress("event",&event);
  Txchk->SetBranchAddress("run",&xrun);

  Tmain->SetBranchAddress("pt",&pt);
  Tmain->SetBranchAddress("eta",&eta);
  Tmain->SetBranchAddress("m1pt",&m1pt);
  Tmain->SetBranchAddress("m2pt",&m2pt);
  Tmain->SetBranchAddress("m1iso",&m1iso);
  Tmain->SetBranchAddress("m2iso",&m2iso);
  Tmain->SetBranchAddress("fl3d",&fl3d);
  Tmain->SetBranchAddress("fls3d",&fls3d);
  Tmain->SetBranchAddress("m",&m);
  Tmain->SetBranchAddress("closetrk",&closetrk);
  Tmain->SetBranchAddress("maxdoca",&maxdoca);
  Tmain->SetBranchAddress("pvip",&pvip);
  Tmain->SetBranchAddress("cosa",&cosa);
  Tmain->SetBranchAddress("docatrk",&docatrk);
  Tmain->SetBranchAddress("iso",&iso);
  Tmain->SetBranchAddress("chi2", &chi2);
  Tmain->SetBranchAddress("dof", &dof);
  Tmain->SetBranchAddress("evt",&evt);
  Tmain->SetBranchAddress("run",&run);

  Comp cpt("xpt", "pt", 100, 0, 50);
  Comp ceta("xeta", "eta", 100, -2.4, 2.4);
  Comp cmu1pt("mu1_pt", "mu1pt", 100, 0, 50);
  Comp cmu2pt("mu2_pt", "mu2pt", 100, 0, 50);
  Comp cmu1iso("mu1_iso", "mu1iso", 100, 0, 10); 
  Comp cmu2iso("mu2_iso", "mu2iso", 100, 0, 10); 
  Comp cl3d("l3d","fl3d",40,0,0.5);
  Comp cl3dsig("l3dsig","fls3d",40,0,10);
  Comp cmass("mass","m",40,4.9,5.9);
  Comp cntrk20("ntrk20","closetrk",40,0,20);
  Comp cdca("dca","maxdoca",40,0,0.03);
  Comp cdelta3d("delta3d","pivp",40,0,0.2);
  Comp ccosAlpha3D("cosAlpha3D","cosa",40,0.9,1);
  Comp cminDca("minDca","docatrk",40,0,0.06);
  Comp cisolation("isolation","iso",40,0,1);
  Comp cchi2dof("chi2dof","NChi2",100,0,10);

  // Int_t nentries_xchk = (Int_t)Txchk->GetEntries();
  Int_t nentries_main = (Int_t)Tmain->GetEntries();

  for (int k=0; k<maxEvents; k++) {
    Txchk->GetEntry(k);
    cout<<event<<" "<<xrun<<" "<<xpt<<" "<<xeta<<" "<<mu1_pt<<" "<<mu2_pt<<" "<<mu1_iso<<" "<<mu2_iso<<" "<<l3d<<" "<<l3dsig<<" "<<mass<<" "<<dca<<" "<<delta3d<<" "<<cosAlpha3D<<" "<<minDca<<" "<<isolation<<" "<<ntrk20<<" "<<NChi2<<endl;
    // cout<<event<<" "<<xpt<<" "<<xeta<<" "<<mu1_pt<<" "<<mu2_pt<<" "<<l3d<<" "<<l3dsig<<" "<<mass<<" "<<dca<<" "<<delta3d<<" "<<cosAlpha3D<<" "<<minDca<<" "<<isolation<<" "<<ntrk20<<" "<<NChi2<<endl;

    bool match = false;
    for (int j=0; j<nentries_main; j++) {
      Tmain->GetEntry(j);
      // cout<<"run: " << run <<endl;
      if (run == xrun && evt==event) {
        //cout<<"found the matched event: "<<endl;
        cout<<evt<<" "<<run<<" "<<pt<<" "<<eta<<" "<<m1pt<<" "<<m2pt<<" "<<m1iso<<" "<<m2iso<<" "<<fl3d<<" "<<fls3d<<" "<<m<<" "<<maxdoca<<" "<<pvip<<" "<<cosa<<" "<<docatrk<<" "<<iso<<" "<<closetrk<<" "<<chi2/dof<<endl<<endl;
        // cout<<evt<<" "<<pt<<" "<<eta<<" "<<m1pt<<" "<<m2pt<<" "<<fl3d<<" "<<fls3d<<" "<<m<<" "<<maxdoca<<" "<<pvip<<" "<<cosa<<" "<<docatrk<<" "<<iso<<" "<<closetrk<<" "<<chi2/dof<<endl<<endl;

        cpt.fill(xpt, pt);
        ceta.fill(xeta, eta);
        if( m1pt > m2pt ) {
          cmu1pt.fill(mu1_pt, m1pt);
          cmu1iso.fill(mu1_iso, m1iso);  
          cmu2pt.fill(mu2_pt, m2pt);
          cmu2iso.fill(mu2_iso, m2iso);
        }
        else {
          cmu1pt.fill(mu1_pt, m2pt);
          cmu1iso.fill(mu1_iso, m2iso);
          cmu2pt.fill(mu2_pt, m1pt);
          cmu2iso.fill(mu2_iso, m1iso);
        }
        cl3d.fill(l3d, fl3d);
        cl3dsig.fill(l3dsig, fls3d);
        cmass.fill(mass, m);
        cntrk20.fill(ntrk20, closetrk);
        cdca.fill(dca, maxdoca);
        cdelta3d.fill(delta3d, pvip);
        ccosAlpha3D.fill(cosAlpha3D, cosa);
        cminDca.fill(minDca, docatrk);
        cisolation.fill(isolation, iso);
        cchi2dof.fill(chi2/dof, NChi2);

        match = true;
        break;
      }
    }

    if (!match) {
      cout<<"!!!!!***** No Matched Event in the Main Analysis Tree *****!!!!!"<<endl<<endl;
      cpt.fill(xpt);
      ceta.fill(xeta);
      cmu1pt.fill(mu1_pt);
      cmu2pt.fill(mu2_pt);
      cmu1iso.fill(mu1_iso);
      cmu2iso.fill(mu2_iso);
      cl3d.fill(l3d);
      cl3dsig.fill(l3dsig);
      cmass.fill(mass);
      cntrk20.fill(ntrk20);
      cdca.fill(dca);
      cdelta3d.fill(delta3d);
      ccosAlpha3D.fill(cosAlpha3D);
      cminDca.fill(minDca);
      cisolation.fill(isolation);
      cisolation.fill(NChi2);
    }
  }

  // Matched
  TCanvas *c = new TCanvas("c","c",1600,1200);
  c->Divide(4,4);
  TCanvas *unmatched_c = new TCanvas("unmatched_c","unmatched c",1600,1200);
  unmatched_c->Divide(4,4);

  TCanvas *cDiff = new TCanvas("cDiff","cDiff",1600,1200);
  cDiff->Divide(2,2);
  const unsigned int nVar = 16;
  TH1F * histoDiff = new TH1F("diff", "diff", nVar, 1, nVar);
  TH1F * histoAbsDiff = new TH1F("absDiff", "absDiff", nVar, 1, nVar);
  TH1F * histoRelDiff = new TH1F("relDiff", "relDiff", nVar, 1, nVar);
  TH1F * histoRelAbsDiff = new TH1F("relAbsDiff", "relAbsDiff", nVar, 1, nVar);
  char *ABC[nVar] = {"pt", "eta", "mu1_pt", "mu2_pt", "fl3d", "fls3d", "m", "closetrk", "maxdoca", "pvip", "cosa", "docatrk", "iso", "chi2dof", "mu1_iso", "mu2_iso"};

  cpt.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 1);
  ceta.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 2);
  cmu1pt.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 3);
  cmu2pt.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 4);
  cl3d.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 5);
  cl3dsig.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 6);
  cmass.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 7);
  cntrk20.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 8);
  cdca.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 9);
  cdelta3d.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 10);
  ccosAlpha3D.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 11);
  cminDca.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 12);
  cisolation.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 13);
  cchi2dof.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 14);
  cmu1iso.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 15);
  cmu2iso.makeplot(c, unmatched_c, histoDiff, histoAbsDiff, histoRelDiff, histoRelAbsDiff, 16);

  c->SaveAs("plots/data_barrel.gif");
  unmatched_c->SaveAs("plots/unmatched_data_barrel.gif");

  cDiff->cd(1);
  histoDiff->Draw();
  for (unsigned int i=1;i<=nVar;i++) histoDiff->GetXaxis()->SetBinLabel(i,ABC[i-1]);
  cDiff->cd(2);
  histoAbsDiff->Draw();
  for (unsigned int i=1;i<=nVar;i++) histoAbsDiff->GetXaxis()->SetBinLabel(i,ABC[i-1]);
  cDiff->cd(3);
  histoRelDiff->Draw();
  for (unsigned int i=1;i<=nVar;i++) histoRelDiff->GetXaxis()->SetBinLabel(i,ABC[i-1]);
  cDiff->cd(4);
  histoRelAbsDiff->Draw();
  for (unsigned int i=1;i<=nVar;i++) histoRelAbsDiff->GetXaxis()->SetBinLabel(i,ABC[i-1]);
}
