#define YieldTree_cxx
#include "YieldTree.h"
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"

#include "commonVar_ups.inc"

//dimuon related histograms
    TTree* t;
    Float_t invariantMass;
    Float_t upsPt;
    Float_t upsRapidity;
    Float_t genUpsPt;
    Float_t genUpsRap;
    Float_t genMinMuPt;
    Float_t genMinMuEta;
    Float_t genPosMuPt;
    Float_t genPosMuEta;

    Float_t muPlusPt;
    Float_t muPlusEta;
    Float_t muMinusPt;
    Float_t muMinusEta;
    Float_t McInvMass;
    Int_t GenRef;
    Int_t SampleFlag;
//=======================================
void YieldTree::Loop(Bool_t removeQQ, Bool_t matchMC, Bool_t printGoodEvents)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  //    nentries = 10000;

  Long64_t nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
 
    if(jentry % 1000 == 0) printf("event %d/%d\n", (Int_t) jentry, (Int_t) nentries);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    SampleFlag = 0;
//    if(removeQQ && Mc_QQ_size > 0) continue; //Mc_QQ_size contains the nb. of MC J/psi's

    //loop over all quarkonia

    //now fill the histos for the "best" QQ pair only
    Int_t iDBestQQ = -1;
    if(Reco_QQ_size > 0){
      iDBestQQ = theBestQQ();
      if(iDBestQQ >= 0){

	TLorentzVector *Reco_QQ = (TLorentzVector *) Reco_QQ_4mom->At(iDBestQQ);
             Double_t mass_Reco_QQ = Reco_QQ->M();
             Double_t pT_Reco_QQ = Reco_QQ->Pt();
             Double_t rap_Reco_QQ = Reco_QQ->Rapidity();

	Bool_t isMatched[2] = {kFALSE, kFALSE};//only relevant for MC; for data must be set to kTRUE
        Bool_t isMatchedHLTMu3[2] = {kFALSE, kFALSE};
        Bool_t isMatchedHLTMu5[2] = {kFALSE, kFALSE};

	  if(Mc_QQmupl_indx[0] < 2 && Mc_QQmumi_indx[0] < 2){ //sanity check
              TLorentzVector * gen_mu_1,  * gen_mu_2;
              gen_mu_1 = (TLorentzVector *) Mc_mu_4mom->At(Mc_QQmupl_indx[0]);
              gen_mu_2 = (TLorentzVector *) Mc_mu_4mom->At(Mc_QQmumi_indx[0]);
             genMinMuPt = gen_mu_2->Pt();
             genMinMuEta = gen_mu_2->Eta();
             genPosMuPt = gen_mu_1->Pt();
             genPosMuEta = gen_mu_1->Eta();

          if(Mc_QQ_size > 0){
             TLorentzVector *Mc_QQ = (TLorentzVector *) Mc_QQ_4mom->At(0);
             McInvMass = Mc_QQ->M();
             genUpsPt = Mc_QQ->Pt();
             genUpsRap = Mc_QQ->Rapidity();
          }   
	    switch (Reco_QQ_type[iDBestQQ]) {
	    case 0: // global-global
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                 if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
                   isMatchedHLTMu3[0] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                 if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR )
                  isMatchedHLTMu5[0] = kTRUE;
                }             
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR )
                  isMatchedHLTMu3[1] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
                  isMatchedHLTMu5[1] = kTRUE;
                }
	      }
	      break;
	    case 1: // global-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_glb_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR )
                  isMatchedHLTMu3[0] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
                  isMatchedHLTMu5[0] = kTRUE;
                }
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR )
                  isMatchedHLTMu3[1] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
                  isMatchedHLTMu5[1] = kTRUE;
                }
	      }
	      break;
	    case 2: // tracker-tracker
	      //first MC muon
	      if(Reco_QQ_muhpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[0] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR )
                  isMatchedHLTMu3[0] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ])) < MAX_deltaR)
                  isMatchedHLTMu5[0] = kTRUE;
                }
	      }
	      //second MC muon
	      if(Reco_QQ_mulpt[iDBestQQ] < Reco_mu_trk_size){//sanity check
               if(Mc_QQ_size > 0){
		if(deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmupl_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR || 
		   deltaR((TLorentzVector*)Mc_mu_4mom->At(Mc_QQmumi_indx[0]), 
			  (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
		  isMatched[1] = kTRUE;
               }
                for(int nHltMu3 = 0; nHltMu3 < HLT1Mu3_L3_size; nHltMu3++){
                if(deltaR((TLorentzVector*)HLT1Mu3_L3_4mom->At(nHltMu3),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR) 
                  isMatchedHLTMu3[1] = kTRUE;
                }
                for(int nHltMu5 = 0; nHltMu5 < HLT1Mu5_L3_size; nHltMu5++){
                if(deltaR((TLorentzVector*)HLT1Mu5_L3_4mom->At(nHltMu5),
                          (TLorentzVector*)Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ])) < MAX_deltaR)
                  isMatchedHLTMu5[1] = kTRUE;
                }
	      }
	      break;
	    default:
	      break;
	    }	    

	  }
//	}
	//check the HLT matching before continuing...
        if((isMatchedHLTMu3[0] == kTRUE || isMatchedHLTMu3[1] == kTRUE) && fabs(rap_Reco_QQ) < 2.0){
            TLorentzVector *Reco_mu_1, *Reco_mu_2;
            if(Reco_QQ_type[iDBestQQ] == 0){//gl-gl
              if(Reco_mu_glb_charge[Reco_QQ_muhpt[iDBestQQ]] == -1){
              Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_2 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }else{
              Reco_mu_2 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }
            }
            else if(Reco_QQ_type[iDBestQQ] == 1){//gl-tr
              if(Reco_mu_glb_charge[Reco_QQ_muhpt[iDBestQQ]] == -1){
              Reco_mu_1 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_2 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }else{
              Reco_mu_2 = (TLorentzVector *) Reco_mu_glb_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_1 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }
            }
            else if(Reco_QQ_type[iDBestQQ] == 2){//tr-tr
              if(Reco_mu_trk_charge[Reco_QQ_muhpt[iDBestQQ]] == -1){
              Reco_mu_1 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_2 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }else{
              Reco_mu_2 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_muhpt[iDBestQQ]);
              Reco_mu_1 = (TLorentzVector *) Reco_mu_trk_4mom->At(Reco_QQ_mulpt[iDBestQQ]);
              }
            }
            Double_t etaMu[2] = {Reco_mu_1->Eta(), Reco_mu_2->Eta()};
            Double_t pTMu[2]  = {Reco_mu_1->Pt(),        Reco_mu_2->Pt()};
          if(Reco_mu_1->Pt() > 3.5 && Reco_mu_2->Pt() >3.5 && fabs(Reco_mu_1->Eta()) <2.4 && fabs(Reco_mu_2->Eta())<2.4){ 
	  //mass histograms
	    invariantMass = mass_Reco_QQ;

	  //fill pT and rap histos only for the signal:
	    if(mass_Reco_QQ > 7.0 && mass_Reco_QQ < 12.0){
	      //pT histograms
	      invariantMass = mass_Reco_QQ;
              upsPt = pT_Reco_QQ;
	      //rap histograms
	      upsRapidity = rap_Reco_QQ;
	    //===============================================================
	    //single muon histograms for best dimuon within 3.0 < M < 3.2 GeV
	    //and for dimuons with M > 2 GeV
	    //===============================================================
	      muPlusPt = pTMu[1];
              muPlusEta = etaMu[1];
              muMinusPt = pTMu[0];
              muMinusEta = etaMu[0];    
            if(isMatched[0] == kTRUE && isMatched[1] == kTRUE) GenRef = 1;
            else GenRef = 0;
            t->Fill();
          }//mass 
         }//mu pt/eta cut
	}//matched: relevant for MC
      }//best QQ
    }//Reco_QQ_size > 0
//   t->Fill();
  }//events
}

//=============================
int YieldTree::theBestQQ() {
    
  int theBest = -1;
  float thehighestPt = -1.;
 
  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 0 ) { //global-global

      int thehptMu = Reco_QQ_muhpt[iqq];   if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];   if (thelptMu >= Reco_mu_glb_size) continue;
      if (Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi && 
//          abs(int(Reco_QQ_lxy[iqq]/Reco_QQ_lxyErr[iqq]))<= 3 &&//lxy
	  Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb && 
	  (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel))&&
// || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) && 
	  fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	  Reco_mu_glb_nhitstrack[thelptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thelptMu] < MAX_normchi2_glb &&
	  (((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel))&&
 //|| ((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thelptMu] == 1)) &&
	  fabs(Reco_mu_glb_d0[thelptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thelptMu]) < MAX_dz_trk
	  ) {
	return iqq;
      }
    }
  }

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 1 ) { //global-tracker
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;
      
      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
//           abs(int(Reco_QQ_lxy[iqq]/Reco_QQ_lxyErr[iqq]))<= 3 &&//lxy
	   Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel))&&
// || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) &&
	   fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,10))/(int)pow(2,10) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,9))/(int)pow(2,9) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk  ) {

	// if pairs containing tracker muons have both mu emmited at |eta| > 2.0
	// check if they are within deltaPhi < 10 degrees (0.17 mrad) and 
	// in this case place request that BOTH have at least two segment matches
	TLorentzVector *glbMu = (TLorentzVector*)Reco_mu_glb_4mom->At(thehptMu);
	TLorentzVector *trkMu = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
/*	if (((glbMu->Eta() > 2.0 && trkMu->Eta() > 2.0) ||
	     (glbMu->Eta() < -2.0 && trkMu->Eta() < -2.0)) &&
	    ((glbMu->Phi() - trkMu->Phi()) < 10.*TMath::Pi()/180. ||
	     ((trkMu->Phi() - glbMu->Phi()) < 10.*TMath::Pi()/180.)) 
	    )
	  continue;
	else{
*/	  TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
	  if (theTrMumom->Perp() > thehighestPt) {
	    thehighestPt = theTrMumom->Perp();
	    theBest = iqq;
	  }
//	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 2 ) { //tracker-tracker
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_trk_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;
      printf("significance %d/\n", (Int_t) (Reco_QQ_lxy[iqq]/Reco_QQ_lxyErr[iqq]));
      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
//           abs(int(Reco_QQ_lxy[iqq]/Reco_QQ_lxyErr[iqq]))<= 3 && //lxy
   	   Reco_mu_trk_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,10))/(int)pow(2,10) > 0 || (Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,9))/(int)pow(2,9) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thehptMu] + Reco_mu_trk_nhitsPixE[thehptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thehptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thehptMu]) < MAX_dz_trk &&
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,10))/(int)pow(2,10) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,9))/(int)pow(2,9) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk 
	   ) {
	
	// if pairs containing tracker muons have both mu emmited at |eta| > 2.0
	// check if they are within deltaPhi < 10 degrees (0.17 mrad) and 
	// in this case place request that BOTH have at least two segment matches
	TLorentzVector *trkMu1 = (TLorentzVector*)Reco_mu_trk_4mom->At(thehptMu);
	TLorentzVector *trkMu2 = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
/*	if (((trkMu1->Eta() > 2.0 && trkMu2->Eta() > 2.0) ||
	     (trkMu1->Eta() < -2.0 && trkMu2->Eta() < -2.0)) &&
	    ((trkMu1->Phi() - trkMu2->Phi()) < 10.*TMath::Pi()/180. ||
	     (trkMu2->Phi() - trkMu1->Phi()) < 10.*TMath::Pi()/180.)
	      )
	  continue;
	else{
*/	  TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thehptMu);
	  if (theTrMumom->Perp() > thehighestPt) {
	    thehighestPt = theTrMumom->Perp();
	    theBest = iqq;
	  }
//	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

}
//=========================================
Double_t YieldTree::deltaR(TLorentzVector* t, TLorentzVector* u){

  return sqrt(pow(t->Eta()-u->Eta(),2) +pow(PhiInRange(t->Phi()-u->Phi()),2));
}
//=========================================
double YieldTree::PhiInRange(double phi){

      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}

