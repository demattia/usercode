#define CosmicMuonAnalyzer_cxx
#include "CosmicMuonAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

inline double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

bool passTrackCut(const Track * track)
{
  // if( track->phi > 0 ) return false;

  // return true;

  // if( track->trackQuality && (fabs(track->eta) < 2.0) && (track->pt > 25) && (track->nValidHits > 6) ) {
  // if( (fabs(track->eta) < 2.0) && (track->pt > 25) && (track->nValidHits > 6) ) {
  if( track->trackQuality ) {
    return true;
  }
  return false;
}

template <class T>
void fillEff(const double & value, TH1F * totalVsD0, const T & tracks, TH1F * passingVsD0)
{
  // Twice because we expect two tracks
  // double muonDxy = fabs((*muons)[0].dxy);
  // double muonDxy = fabs(returnValue->get((*muons)[0]));
  totalVsD0->Fill(value);
  totalVsD0->Fill(value);

  if( tracks->size() > 0 ) {
    bool firstPass = true;
    double firstPhi = 0.;
    std::vector<Track>::const_iterator it = tracks->begin();
    for( ; it != tracks->end(); ++it ) {
      // Always check that the second track (if any) is back to back. Discard tracks that are close by.
      // We are assuming the back-to-back topology of cosmic tracks.
      // PiOver2 just to check that they are not very close. The deltaPhi should be either ~0 or ~pi.
      if( firstPass || fabs(deltaPhi(firstPhi, it->phi)) > TMath::PiOver2() ) {
	if( passTrackCut( &*it ) ) {
	  passingVsD0->Fill(value);
	  if( firstPass ) firstPhi = it->phi;
	  else break;
	  firstPass = false;
	}
      }
    }
  }
}

double polynomial( const double & pt )
{
  return ( 0.1+1.91364 - 0.0211496*pt + 0.0000906055*pt*pt - 0.000000130650*pt*pt*pt );
}

/// Returns the maximum dxy error value allowed for the cuts
double dxyErrMax( const double & pt )
{
  double dxyErrMax = 1.;
  if(pt < 200 ) dxyErrMax = polynomial(pt);
  else dxyErrMax = polynomial(200);
  return std::min(dxyErrMax, 1.);
}


bool passMuonCut(const Track * muon)
{
  if( (muon->nValidHits >= 0)
      && ( muon->dtStationsWithValidHits + muon->cscStationsWithValidHits > 1 )
      && (fabs(muon->eta) < 2.)
      && (fabs(muon->dxyError) < dxyErrMax(muon->pt))
      && (fabs(muon->dzError) < dxyErrMax(muon->pt)) // Note: the use of the same function is intentional.
      ) {

    if( fabs(muon->dz) < 10. ) {
      return true;
    }
  }
  return false;
}

void CosmicMuonAnalyzer::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  int nBins = 50;

  TFile * outputFile = new TFile("eff.root", "RECREATE");
  outputFile->cd();

  // Vs D0
  TH1F * passingVsD0 = new TH1F("passingVsD0", "tracks found vs |d0|", nBins, 0., 100.);
  TH1F * totalVsD0 = new TH1F("totalVsD0", "cosmics found vs |d0| (filled twice)", nBins, 0., 100.);
  TH1F * effVsD0 = (TH1F*)passingVsD0->Clone("effVsD0");

  TH1F * passingVsGenD0 = 0;
  TH1F * totalVsGenD0 = 0;
  TH1F * effVsGenD0 = 0;
  if( MC ) {
    passingVsGenD0 = new TH1F("passingVsGenD0", "tracks found vs |genD0|", nBins, 0., 100.);
    totalVsGenD0 = new TH1F("totalVsGenD0", "genParticles found vs |genD0| (filled twice)", nBins, 0., 100.);
    effVsGenD0 = (TH1F*)passingVsGenD0->Clone("effVsGenD0");
  }

  // Vs Eta
  TH1F * passingVsEta = new TH1F("passingVsEta", "tracks found vs #eta", nBins, -2.5, 2.5);
  TH1F * totalVsEta = new TH1F("totalVsEta", "cosmics found vs #eta (filled twice)", nBins, -2.5, 2.5);
  TH1F * effVsEta = (TH1F*)passingVsEta->Clone("effVsEta");

  TH1F * passingVsGenEta = 0;
  TH1F * totalVsGenEta = 0;
  TH1F * effVsGenEta = 0;
  if( MC ) {
    passingVsGenEta = new TH1F("passingVsGenEta", "tracks found vs gen #eta", nBins, -2.5, 2.5);
    totalVsGenEta = new TH1F("totalVsGenEta", "genParticles found vs gen #eta (filled twice)", nBins, -2.5, 2.5);
    effVsGenEta = (TH1F*)passingVsGenEta->Clone("effVsGenEta");
  }

  TH1F * etaRelError = new TH1F("etaRelError", "relative error on #eta", nBins, 0., 1.5);
  TH1F * dxyRelError = new TH1F("dxyRelError", "relative error on d_{0}", nBins, 0., 1.5);
  TH1F * dzRelError = new TH1F("dzRelError", "relative error on z_{0}", nBins, 0., 1.5);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( jentry % 1000 == 0 ) {
      std::cout << "Analyzed " << jentry << " events" << std::endl;
    }

    // Compute efficiency by counting how many tracks are reconstructed when two (opposite) are expected
    // Check this efficiency for tracks passing the analysis cuts.
    // Use only events in which exactly one cosmic1LegMuon has been reconstructed.
    if( muons->size() != 1 ) continue;
    if( !passMuonCut(&((*muons)[0])) ) continue;

    // Efficiency vs d0
    fillEff(fabs((*muons)[0].dxy), totalVsD0, tracks, passingVsD0);
    if( MC ) {
      fillEff(fabs((*genParticles)[0].dxy), totalVsGenD0, tracks, passingVsGenD0);
    }

    // Apply also a dxy cut for the efficiency vs eta
    if( (*muons)[0].dxy > 5. ) continue;
    // Efficiency vs eta
    fillEff((*muons)[0].eta, totalVsEta, tracks, passingVsEta);
    if( MC ) {
      fillEff((*genParticles)[0].eta, totalVsGenEta, tracks, passingVsGenEta);
    }

    etaRelError->Fill((*muons)[0].etaError/fabs((*muons)[0].eta));
    dxyRelError->Fill((*muons)[0].dxyError/fabs((*muons)[0].dxy));
    dzRelError->Fill((*muons)[0].dzError/fabs((*muons)[0].dz));

  }

  // Efficiency vs d0
  passingVsD0->Sumw2();
  totalVsD0->Sumw2();
  effVsD0->Divide(passingVsD0, totalVsD0, 1., 1., "B");
  if( MC ) {
    passingVsGenD0->Sumw2();
    totalVsGenD0->Sumw2();
    effVsGenD0->Divide(passingVsGenD0, totalVsGenD0, 1., 1., "B");
  }
  // Efficiency vs eta
  passingVsEta->Sumw2();
  totalVsEta->Sumw2();
  effVsEta->Divide(passingVsEta, totalVsEta, 1., 1., "B");
  if( MC ) {
    passingVsGenEta->Sumw2();
    totalVsGenEta->Sumw2();
    effVsGenEta->Divide(passingVsGenEta, totalVsGenEta, 1., 1., "B");
  }

  outputFile->Write();
  outputFile->Close();
}
