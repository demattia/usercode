#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <iostream>

#define NMAX 50

class ScaleCorrection
{
public:

  ScaleCorrection() :
    ms_a0(3.8e-4), ms_a0_e(1.9e-4), ms_a1(0.), ms_a1_e(0.),
    ms_a2(3.0e-4), ms_a2_e(0.7e-4), ms_a3(0.), ms_a3_e(0.)
  {
  }

  inline double scale(const double & mupt, const double & mueta, const int ptscale)
  {
    if( ptscale == 0 ) {
      return (1.0 + ms_a0 + ms_a1*fabs(mueta) + ms_a2*mueta*mueta + ms_a3*mupt)*mupt;
    }
    else if( ptscale == 1 ) {
      return (1.0 + (ms_a0 + ms_a0_e) + (ms_a1 + ms_a1_e)*fabs(mueta) + (ms_a2 + ms_a2_e)*mueta*mueta + (ms_a3 + ms_a3_e)*mupt)*mupt;
    }
    else if( ptscale == 2 ) {
      return (1.0 + (ms_a0 - ms_a0_e) + (ms_a1 - ms_a1_e)*fabs(mueta) + (ms_a2 - ms_a2_e)*mueta*mueta + (ms_a3 - ms_a3_e)*mupt)*mupt;
    }
  }

private:
  //mass scale correction
  double ms_a0;
  double ms_a0_e;
  double ms_a1;
  double ms_a1_e;
  double ms_a2;
  double ms_a2_e;
  double ms_a3;
  double ms_a3_e;
};

void acceptance(const std::string ptscale = "", const char* file = "../MCupsilonTree_pt0-50.root", const char* outputName = "acceptance0_50_newAccCuts_basehistos", int num=1) {

  ScaleCorrection scaleCorrection;

  int ptscaleType = -1;
  if( ptscale == "scale" ) {
    std::cout << "scale corrections will be applied" << std::endl;
    ptscaleType = 0;
  }
  else if( ptscale == "scaleUp" ) {
    std::cout << "scale + errors corrections will be applied" << std::endl;
    ptscaleType = 1;
  }
  else if( ptscale == "scaleDown" ) {
    std::cout << "scale - errors corrections will be applied" << std::endl;
    ptscaleType = 2;
  }
  else if( ptscale != "" ) {
    std::cout << "Error: ptscale mode " << ptscale << " invalid. Available modes are:" << std::endl;
    std::cout << "- \"\": (empty) for no corrections" << std::endl;
    std::cout << "- \"scale\": apply the scale corrections" << std::endl;
    std::cout << "- \"scaleUp\": apply the scale + errors corrections" << std::endl;
    std::cout << "- \"scaleDown\": apply the scale - errors corrections" << std::endl;
    exit(1);
  }
  else {
    std::cout << "NO scale corrections will be applied" << std::endl;
  }

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

  cout << "entries=" << t->GetEntries() << endl;


  TH2F* genPtRap1 = new TH2F("genUps1","",50,0,50,30,0,3.0);
  TH2F* genPtRap2 = new TH2F("genUps2","",50,0,50,30,0,3.0);
  TH2F* genPtRap3 = new TH2F("genUps3","",50,0,50,30,0,3.0);
  TH2F* genPtRap4 = new TH2F("genUps4","",50,0,50,30,0,3.0);
  TH2F* genPtRap5 = new TH2F("genUps5","",50,0,50,30,0,3.0);

  genPtRap1->Sumw2();
  genPtRap2->Sumw2();
  genPtRap3->Sumw2();
  genPtRap4->Sumw2();
  genPtRap5->Sumw2();

  genPtRap1->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");
  genPtRap2->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");
  genPtRap3->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");
  genPtRap4->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");
  genPtRap5->SetTitle(";Upsilon pT (GeV/c);Upsilon Rapidity;");

  TH2F* recoPtRap1 = (TH2F*)genPtRap1->Clone("recoUps1");
  TH2F* recoPtRap2 = (TH2F*)genPtRap2->Clone("recoUps2");
  TH2F* recoPtRap3 = (TH2F*)genPtRap3->Clone("recoUps3");
  TH2F* recoPtRap4 = (TH2F*)genPtRap4->Clone("recoUps4");
  TH2F* recoPtRap5 = (TH2F*)genPtRap5->Clone("recoUps5");

  recoPtRap1->Sumw2();
  recoPtRap2->Sumw2();
  recoPtRap3->Sumw2();
  recoPtRap4->Sumw2();
  recoPtRap5->Sumw2();

  Double_t ptBins[15] = {0,1.5,3,4,5,6,7,8,9,10,12,14,17,20,50};

  TH1F* genPt1  = new TH1F("genPt1","",14,ptBins);
  genPt1->Sumw2();
  TH1F* recoPt1 = (TH1F*)genPt1->Clone("recoUpsPt1");
  recoPt1->Sumw2();
  TH1F* genPt2  = new TH1F("genPt2","",14,ptBins);
  genPt2->Sumw2();
  TH1F* recoPt2 = (TH1F*)genPt2->Clone("recoUpsPt2");
  recoPt2->Sumw2();
  TH1F* genPt3  = new TH1F("genPt3","",14,ptBins);
  genPt3->Sumw2();
  TH1F* recoPt3 = (TH1F*)genPt3->Clone("recoUpsPt3");
  recoPt3->Sumw2();
  TH1F* genPt4  = new TH1F("genPt4","",14,ptBins);
  genPt4->Sumw2();
  TH1F* recoPt4 = (TH1F*)genPt4->Clone("recoUpsPt4");
  recoPt4->Sumw2();
  TH1F* genPt5  = new TH1F("genPt5","",14,ptBins);
  genPt5->Sumw2();
  TH1F* recoPt5 = (TH1F*)genPt5->Clone("recoUpsPt5");
  recoPt5->Sumw2();

  Float_t w1,w2,w3,w4,w5;


  double E=30; double pz = sqrt(E*E - 0.938272*0.938272);
  TLorentzVector h1; h1.SetPxPyPzE(0,0,pz,E);
  TLorentzVector h2; h2.SetPxPyPzE(0,0,-pz,E);
  TLorentzVector genUps;
  TLorentzVector genMuPlus;
  int mp;
  Float_t cosThetaStarHel;
  TVector3 zCS;
  Float_t cosThetaStarCS;



  for(int i=0; i<t->GetEntries(); i++){
    if(i%10000 == 0)
      std::cout<<i<<std::endl;
    t->GetEntry(i);
    double genupspt=genUpsPt[0];
    double genupsrap=fabs(genUpsRapidity[0]);

    // calculate cosTheta helicity
    genUps.SetPtEtaPhiM(genupspt, genUpsEta[0], genUpsPhi[0], 9.46);
    TLorentzRotation boost(-genUps.BoostVector());
    mp = genMuCharge[0]>0 ? 0 : 1;
    genMuPlus.SetPtEtaPhiM(genMuPt[mp], genMuEta[mp], genMuPhi[mp], 0.106);
    genMuPlus *= boost;
    cosThetaStarHel = genMuPlus.Vect().Dot(genUps.Vect())/(genMuPlus.Vect().Mag()*genUps.Vect().Mag());

    // calculate cosTheta CS
    h1.SetPxPyPzE(0,0,pz,E);
    h2.SetPxPyPzE(0,0,-pz,E);
    h1*=boost;
    h2*=boost;
    zCS = ( h1.Vect().Unit() - h2.Vect().Unit() ).Unit();
    cosThetaStarCS = genMuPlus.Vect().Dot(zCS)/genMuPlus.Vect().Mag();

    // setup the weights
    w1 = 1;
    w2 = 1 + cosThetaStarHel*cosThetaStarHel;
    w3 = 1 - cosThetaStarHel*cosThetaStarHel;
    w4 = 1 + cosThetaStarCS*cosThetaStarCS;
    w5 = 1 - cosThetaStarCS*cosThetaStarCS;

    genPtRap1->Fill( genupspt, genupsrap, w1 );
    genPtRap2->Fill( genupspt, genupsrap, w2 );
    genPtRap3->Fill( genupspt, genupsrap, w3 );
    genPtRap4->Fill( genupspt, genupsrap, w4 );
    genPtRap5->Fill( genupspt, genupsrap, w5 );
    genPt1->Fill( genupspt,w1 );
    genPt2->Fill( genupspt,w2 );
    genPt3->Fill( genupspt,w3 );
    genPt4->Fill( genupspt,w4 );
    genPt5->Fill( genupspt,w5 );

    Float_t recoUpsPt = 0;
    Float_t recoUpsRapidity = 0;
    double minDeltaM = 1000;
    double mass =0;
    double tmpmass=0;
    double mupt1,mupt2,mueta1,mueta2;
    double sign;
    int match1,match2;

    
    if(recoMuSize < 2) continue;

    for(int tr1=0; tr1<recoMuSize; tr1++) {
      mupt1=recoMuPt[tr1];
      mueta1=recoMuEta[tr1];
      // apply scale corrections if needed
      if( ptscaleType != -1 ) {
	mupt1 = scaleCorrection.scale(mupt1, mueta1, ptscaleType);
      }
      if(  (mupt1 > 3.75 && fabs(mueta1) < 0.8) || (mupt1 > 3.5 && fabs(mueta1) >= 0.8 && fabs(mueta1) < 1.6) || (mupt1 > 3.0 && fabs(mueta1) >= 1.6 && fabs(mueta1) < 2.4) ) {
	for(int tr2=tr1+1; tr2<recoMuSize; tr2++){
	  mupt2=recoMuPt[tr2];
	  mueta2=recoMuEta[tr2];
	  // apply scale corrections if needed
	  if( ptscaleType != -1 ) {
	    mupt2 = scaleCorrection.scale(mupt2, mueta2, ptscaleType);
	  }
	  if ( recoMuCharge[tr1]*recoMuCharge[tr2] == -1 && ((mupt2 > 3.75 && fabs(mueta2) < 0.8) || (mupt2 > 3.5 && fabs(mueta2) >= 0.8 && fabs(mueta2) < 1.6) || (mupt2 > 3.0 && fabs(mueta2) >= 1.6 && fabs(mueta2) < 2.4))) {
	    TLorentzVector mu1; mu1.SetPtEtaPhiM(mupt1, mueta1, recoMuPhi[tr1], 0.1057);
	    TLorentzVector mu2; mu2.SetPtEtaPhiM(mupt2, mueta2, recoMuPhi[tr2], 0.1057);
	    TLorentzVector recoUps(mu1 + mu2);
	    tmpmass = recoUps.M();
 	    double deltaM = fabs(tmpmass-9.46);
	    //		cout<<"deltaM"<<deltaM<<endl;
 	    if( deltaM < minDeltaM ){
 	      recoUpsPt = recoUps.Pt();
 	      recoUpsRapidity = fabs(recoUps.Rapidity());
 	      minDeltaM = deltaM;
	      mass = tmpmass;
	    }
	  }
	}
      }
    }

    if( mass > 8.0 && mass < 12.0){
      recoPtRap1->Fill( recoUpsPt, recoUpsRapidity, w1 );
      recoPtRap2->Fill( recoUpsPt, recoUpsRapidity, w2 );
      recoPtRap3->Fill( recoUpsPt, recoUpsRapidity, w3 );
      recoPtRap4->Fill( recoUpsPt, recoUpsRapidity, w4 );
      recoPtRap5->Fill( recoUpsPt, recoUpsRapidity, w5 );
      recoPt1->Fill( recoUpsPt,w1);	
      recoPt2->Fill( recoUpsPt,w2);
      recoPt3->Fill( recoUpsPt,w3);
      recoPt4->Fill( recoUpsPt,w4);
      recoPt5->Fill( recoUpsPt,w5);
    }
  
    

  } // loop over tree entries

  

  char outfilename[100];
  sprintf(outfilename,"%s_%i.root",outputName,num);
  TFile out(outfilename,"recreate");

  genPtRap1->Write();
  genPtRap2->Write();
  genPtRap3->Write();
  genPtRap4->Write();
  genPtRap5->Write();
  genPt1->Write();
  genPt2->Write();
  genPt3->Write();
  genPt4->Write();
  genPt5->Write();

  recoPtRap1->Write();
  recoPtRap2->Write();
  recoPtRap3->Write();
  recoPtRap4->Write();
  recoPtRap5->Write();
  recoPt1->Write();
  recoPt2->Write();
  recoPt3->Write();
  recoPt4->Write();
  recoPt5->Write();

  out.Close();
}

