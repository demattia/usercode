#define AnalyzeTree_cxx
#include "AnalyzeTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMath.h>
#include <map>

inline double deltaPhi(const double & phi1, const double & phi2)
{
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

inline double deltaR(const double & eta1, const double & eta2, const double & phi1, const double & phi2)
{
  return sqrt(pow(eta1-eta2, 2) + pow(deltaPhi(phi1,phi2), 2));
}

/// Take them by copy so that they can be manipulated inside
template <class T1, class T2>
std::map<const T1*, const T2*> associate(const std::vector<T1> * tracks, const std::vector<T2> * genParticles)
{
  typename std::map<const T1*, const T2*> associationMap;
  double deltaRCut = 0.1;

  // Copy to a vector of pointers
  typename std::vector<const T2*> genParticlePtrs;
  for( typename std::vector<T2>::const_iterator gen = genParticles->begin(); gen != genParticles->end(); ++gen ) {
    genParticlePtrs.push_back(&*gen);
  }

  typename std::vector<T1>::const_iterator tk = tracks->begin();
  for( ; tk != tracks->end(); ++tk ) {
    double deltaRMin = 1000.;
    const T2 * genMatch = 0;
    typename std::vector<const T2*>::iterator gen = genParticlePtrs.begin();
    for( ; gen != genParticlePtrs.end(); ++gen ) {
      double deltaRValue = deltaR((*gen)->eta, tk->eta, (*gen)->phi, tk->phi);
      if( deltaRValue < deltaRCut && deltaRValue < deltaRMin ) {
        deltaRMin = deltaRValue;
        genMatch = *gen;
      }
    }
    associationMap.insert(std::pair<const T1*, const T2*>(&*tk, genMatch));
    // Remove the associated genParticle from the vector
    genParticlePtrs.erase(std::remove(genParticlePtrs.begin(), genParticlePtrs.end(), genMatch), genParticlePtrs.end());
  }
  return associationMap;
}


void AnalyzeTree::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  TFile * outputFile = new TFile("histos_cut_old.root", "RECREATE");
  // TFile * outputFile = new TFile("histos_cut_new.root", "RECREATE");
  TH1F * hDeltaR = new TH1F("hDeltaR", "#Delta R(reco, gen)", 100, 0, 0.2);
  TH1F * hGenD0 = new TH1F("hGenD0", "gen vs d0", 100, 0, 400);
  TH1F * hGenPt = new TH1F("hGenPt", "gen vs p_{T}", 110, 0, 110);
  TH1F * hPassingVsD0 = new TH1F("hPassingVsD0", "passing vs d0", 100, 0, 400);
  TH1F * hFailingVsD0 = new TH1F("hFailingVsD0", "failing vs d0", 100, 0, 400);
  TH1F * hTotalVsD0 = new TH1F("hTotalVsD0", "total vs D0", 100, 0, 400);

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;


    for( std::vector<Track>::const_iterator it = tracks->begin(); it != tracks->end(); ++it ) {
      if( it->nValidHits > 1 &&  it->pt > 23 ) {
	hPassingVsD0->Fill(it->dxy);
      }
      else {
	hFailingVsD0->Fill(it->dxy);
      }
    }

    for( std::vector<GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
      hGenD0->Fill(it->dxy);
      hGenPt->Fill(it->pt);
    }


    std::map<const GenParticle*, const Track*> associationMap = associate(genParticles, tracks);
    
    for( std::map<const GenParticle*, const Track*>::const_iterator it = associationMap.begin();
         it != associationMap.end(); ++it ) {
      // Skip unassociated tracks
      if( it->second != 0 ) {
    	// Use this to avoid any cut
    	// if( true ) {
    	if( it->second->nHits >= 6 && it->second->pt > 23 ) {
    	  hDeltaR->Fill(deltaR(it->first->eta, it->second->eta, it->first->phi, it->second->phi));
    	  // hPassingVsD0->Fill(it->first->dxy);
    	}
    	// else {
    	//   hFailingVsD0->Fill(it->first->dxy);
    	// }
      }
    }
  }
  outputFile->cd();
  hDeltaR->Write();
  hPassingVsD0->Write();
  hFailingVsD0->Write();

  hTotalVsD0->Add(hPassingVsD0, hFailingVsD0);
  TH1F * hEffVsD0 = (TH1F*)hPassingVsD0->Clone("hEffVsD0");
  hEffVsD0->Divide(hTotalVsD0);
  hEffVsD0->Write();

  outputFile->Write();
  outputFile->Close();
}
