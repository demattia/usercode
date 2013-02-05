#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <iostream>

#define NMAX 100

void acceptance(const char* file = "upsilonGun.root", const char* outputName = "acceptance.root") {
  Int_t genUpsSize;
  Float_t genUpsPt[NMAX];
  Float_t genUpsEta[NMAX];
  Float_t genUpsPhi[NMAX];
  Float_t genUpsRapidity[NMAX];
  Int_t genMuSize;
  Float_t genMuPt[NMAX];
  Float_t genMuEta[NMAX];
  Float_t genMuPhi[NMAX];
  Int_t genMuCharge[NMAX];
  Int_t recoMuSize;
  Float_t recoMuPt[NMAX];
  Float_t recoMuEta[NMAX];
  Float_t recoMuPhi[NMAX];
  Int_t recoMuCharge[NMAX];

  TFile* f1 = new TFile(file);
  TTree* t = (TTree*)f1->Get("UpsTree");
  t->SetBranchAddress("genUpsSize",&genUpsSize);
  t->SetBranchAddress("genUpsPt",genUpsPt);
  t->SetBranchAddress("genUpsEta",genUpsEta);
  t->SetBranchAddress("genUpsPhi",genUpsPhi);
  t->SetBranchAddress("genUpsRapidity",genUpsRapidity);
  t->SetBranchAddress("genMuSize",&genMuSize);
  t->SetBranchAddress("genMuPt",genMuPt);
  t->SetBranchAddress("genMuEta",genMuEta);
  t->SetBranchAddress("genMuPhi",genMuPhi);
  t->SetBranchAddress("genMuCharge",genMuCharge);
  t->SetBranchAddress("recoMuSize",&recoMuSize);
  t->SetBranchAddress("recoMuPt",recoMuPt);
  t->SetBranchAddress("recoMuEta",recoMuEta);
  t->SetBranchAddress("recoMuPhi",recoMuPhi);
  t->SetBranchAddress("recoMuCharge",recoMuCharge);

  TH2F* genPtRap = new TH2F("genUps","",4,0,20,4,-2,2);
  TH1F* genPt = new TH1F("genUps1D","",4,0,20);
  genPtRap->Sumw2();
  genPt->Sumw2();
  genPtRap->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");
  genPt->SetTitle(";Upsilon pT (GeV/c);Acceptance;");
  TH2F* recoPtRap = (TH2F*)genPtRap->Clone("recoUps");
  TH1F* recoPt = (TH1F*)genPt->Clone("recoUps1D");
  recoPtRap->Sumw2();
  recoPt->Sumw2();

  for(int i=0; i<t->GetEntries(); i++){
    if(i%100000 == 0)
      std::cout<<i<<std::endl;
    t->GetEntry(i);

    // calculate cosTheta
    TLorentzVector genUps; genUps.SetPtEtaPhiM(genUpsPt[0], genUpsEta[0], genUpsPhi[0], 9.46);
    TLorentzRotation boost(-genUps.BoostVector());
    int mp = genMuCharge[0]>0 ? 0 : 1;
    TLorentzVector genMuPlus; genMuPlus.SetPtEtaPhiM(genMuPt[mp], genMuEta[mp], genMuPhi[mp], 0.106);
    genMuPlus *= boost;
    Float_t cosThetaStar = genMuPlus.Vect().Dot(genUps.Vect())/genMuPlus.Vect().Mag()/genUps.Vect().Mag();

    // set the weight
    Float_t weight = 1;

    genPtRap->Fill( genUpsPt[0], genUpsRapidity[0], weight );
    genPt->Fill( genUpsPt[0], weight );

    Float_t recoUpsPt = 0;
    Float_t recoUpsRapidity = 0;
    double minDeltaM = 1000;
    for(int tr1=0; tr1<recoMuSize; tr1++){
      for(int tr2=tr1+1; tr2<recoMuSize; tr2++){
        if ( recoMuCharge[tr1]*recoMuCharge[tr2] == -1 && recoMuPt[tr1] >= 3.5 && recoMuPt[tr2] >= 3.5 && fabs(recoMuEta[tr1]) < 2.1 && fabs(recoMuEta[tr2]) < 2.1 ){
          TLorentzVector mu1; mu1.SetPtEtaPhiM(recoMuPt[tr1], recoMuEta[tr1], recoMuPhi[tr1], 0.106);
          TLorentzVector mu2; mu2.SetPtEtaPhiM(recoMuPt[tr2], recoMuEta[tr2], recoMuPhi[tr2], 0.106);
          TLorentzVector recoUps(mu1 + mu2);
          double deltaM = fabs(recoUps.M()-9.46);
          if( deltaM < minDeltaM ){
            recoUpsPt = recoUps.Pt();
            recoUpsRapidity = recoUps.Rapidity();
            minDeltaM = deltaM;
          }
        }
      }
    }
    if( minDeltaM < 1.0 ){
      recoPtRap->Fill( recoUpsPt, recoUpsRapidity, weight );
      recoPt->Fill( recoUpsPt, weight );
    }
  }

  TFile out(outputName,"recreate");
  TH2F* acc = (TH2F*)genPtRap->Clone("acceptance");
  TH1F* acc1D = (TH1F*)genPt->Clone("acceptance1D");
  acc->Sumw2();
  acc1D->Sumw2();
  acc->Divide(recoPtRap,genPtRap,1,1,"B");
  acc1D->Divide(recoPt,genPt,1,1,"B");
  acc->Write();
  acc1D->Write();
  out.Close();
}

int main(int argc, const char** argv){
  acceptance(argv[1],argv[2]);
}
