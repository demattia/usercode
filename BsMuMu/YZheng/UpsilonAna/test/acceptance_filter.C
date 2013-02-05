#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#endif

void acceptance_filter(const char* file = "gen_mix_wofilter.root", const char* outputName = "acceptance_filter_mixed.root") {
  gROOT->Reset();
  Double_t ptBins[18] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,30.,1e08};
  Double_t ETABIN[11] = {-2.0,-1.8,-1.6,-1.4,-1.0,0.0,1.0,1.4,1.6,1.8,2.0};
  Double_t PTBins[18] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,30.,50.};
  genPtRap = new TH2F("genUps","",100,0,50,60,-3,3);
  genPtRap->Sumw2();
  genPtRap_binned = new TH2F("genUps_binned","",17,ptBins,10,ETABIN);
  genPtRap_binned->Sumw2();
  TH2F* recoPtRap = (TH2F*)genPtRap->Clone("recoUps");
  TH2F* recoPtRap_binned = (TH2F*)genPtRap_binned->Clone("recoUps_binned");

  h_ptUps_binned = new TH1F ("ptUps_binned","ptUps_binned",17,PTBins);
  h_ptUps_binned->Sumw2();
  h_ptUpsFilter_binned = new TH1F ("ptUpsFilter_binned","ptUpsFilter_binned",17,PTBins);
  h_ptUpsFilter_binned->Sumw2();
  h_ptUps = new TH1F ("ptUps","ptUps",100,0.0,50.0);
  h_ptUps->Sumw2();
  h_ptUpsFilter = new TH1F ("ptUpsFilter","ptUpsFilter",100,0.0,50.0);
  h_ptUpsFilter->Sumw2();

   Int_t genUpsSize;
   Float_t genUpsPt[2];
   Float_t genUpsEta[2];
   Float_t genUpsPhi[2];
   Float_t genUpsRap[2];

   Int_t genMuSize;
   Float_t genMuPt[2];
   Float_t genMuEta[2];
   Float_t genMuPhi[2];
   Int_t genMuCharge[2];
   Int_t recoMuSize;
   Float_t recoMuPt[2];
   Float_t recoMuEta[2];
   Float_t recoMuPhi[2];
   Int_t recoMuCharge[2];

  TFile* f1 = new TFile(file);

  f1->cd("");

  TTree* t = (TTree*)gDirectory->Get("UpsTree");
  t->SetBranchAddress("genMuSize",&genMuSize);
  t->SetBranchAddress("genUpsSize",&genUpsSize);
  t->SetBranchAddress("recoMuSize",&recoMuSize);
  t->SetBranchAddress("genMuPt",genMuPt);
  t->SetBranchAddress("genMuEta",genMuEta);
  t->SetBranchAddress("genMuPhi",genMuPhi);
  t->SetBranchAddress("genMuCharge",genMuCharge);
  t->SetBranchAddress("recoMuPt",recoMuPt);
  t->SetBranchAddress("recoMuEta",recoMuEta);
  t->SetBranchAddress("recoMuPhi",recoMuPhi);
  t->SetBranchAddress("recoMuCharge",recoMuCharge);
  t->SetBranchAddress("genUpsPt",genUpsPt);
  t->SetBranchAddress("genUpsEta",genUpsEta);
  t->SetBranchAddress("genUpsPhi",genUpsPhi);
  t->SetBranchAddress("genUpsRapidity",genUpsRap);

  Int_t nentries = 0;

  nentries = Int_t(t->GetEntries());
  cout<<"nentries "<<nentries<<endl;
  for ( int i = 0; i < nentries; i++ ) {
    t->GetEntry(i);
    if(i % 1000 == 0) printf("event %d\n", i);
    TLorentzVector genUps; genUps.SetPtEtaPhiM(genUpsPt[0], genUpsEta[0], genUpsPhi[0], 9.46);
    TLorentzRotation boost(-genUps.BoostVector());
    int mp = genMuCharge[0]>0 ? 0 : 1;
    TLorentzVector genMuPlus; genMuPlus.SetPtEtaPhiM(genMuPt[mp], genMuEta[mp], genMuPhi[mp], 0.106);
    genMuPlus *= boost;
    Float_t cosThetaStar = genMuPlus.Vect().Dot(genUps.Vect())/genMuPlus.Vect().Mag()/genUps.Vect().Mag();
    Float_t alpha = 1;//0,-1,1,-0.34,0.34
    Float_t weight = 1 + alpha * cosThetaStar * cosThetaStar;
    h_ptUps_binned -> Fill(genUpsPt[0],weight);
    h_ptUps -> Fill(genUpsPt[0],weight);
    genPtRap->Fill( genUpsPt[0], genUpsRap[0], weight );
    genPtRap_binned->Fill( genUpsPt[0], genUpsRap[0], weight );
    if( genMuSize ==2 && genMuPt[0] > 2.5 && genMuPt[1] > 2.5 && fabs(genMuEta[0]) < 2.5 && fabs(genMuEta[1]) < 2.5){
        h_ptUpsFilter -> Fill(genUpsPt[0],weight);
        h_ptUpsFilter_binned -> Fill(genUpsPt[0],weight);
        recoPtRap->Fill( genUpsPt[0], genUpsRap[0], weight );
        recoPtRap_binned->Fill( genUpsPt[0], genUpsRap[0], weight );
    }
  }//event

  TH2F* acc = (TH2F*)genPtRap->Clone("acceptance2D");
  acc->Divide(recoPtRap,genPtRap,1,1,"B");
  TH1F* acc1d = (TH1F*)h_ptUps->Clone("acceptance1D");
  acc1d->Divide(h_ptUpsFilter,h_ptUps,1,1,"B");
  TH2F* acc_binned = (TH2F*)genPtRap_binned->Clone("acceptance2D_binned");
  acc_binned->Divide(recoPtRap_binned,genPtRap_binned,1,1,"B");
  TH1F* acc1d_binned  = (TH1F*)h_ptUpsFilter_binned->Clone("acceptance1D_binned ");
  acc1d_binned->Divide(h_ptUpsFilter_binned,h_ptUps_binned, 1,1,"B");

  TFile *output =  new TFile(outputName,"RECREATE");
  output->cd();
  acc->Write();
  acc1d->Write();
  acc1d_binned->Write();
  acc_binned->Write();
  output->Write();
  output->Close();
}

