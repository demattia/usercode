#define treeAnalyzer_cxx
#include "treeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TLorentzVector.h>
#include "Plots.h"

TLorentzVector convert(const double & pt, const double & eta, const double & phi, const double & mass)
{
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  double tmp = 2*atan(exp(-eta));
  double pz = pt*cos(tmp)/sin(tmp);
  double E  = sqrt(px*px+py*py+pz*pz+mass*mass);

  return TLorentzVector(px,py,pz,E);
}

void treeAnalyzer::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Plots beforeAllCuts("_beforeAllCuts");
  Plots afterTriggerAndMETCut("_afterTriggerAndMETCut");
  Plots afterAllCuts("_afterAllCuts");

  TH1F * hZLeptonPtH = new TH1F("hZLeptonPtH", "Z lepton pt high", 100, 0., 100.);
  TH1F * hZLeptonPtL = new TH1F("hZLeptonPtL", "Z lepton pt low", 100, 0., 100.);
  TH1F * hThirdLeptonPt = new TH1F("hThirdLeptonPt", "pt of the third lepton", 100, 0., 100.);

  std::map<std::string, int> effCounter;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    bool triggered = false;
    std::vector<std::string>::const_iterator tr = triggers->begin();
    for( ; tr != triggers->end(); ++tr ) {
      if( tr->find("HLT_Mu17_Mu8_v") != std::string::npos ) {
	triggered = true;
      }
    }
    beforeAllCuts.Fill(0., 0., 0., candidates->MET, 0., weight_);

    effCounter["total"] += 1;

    if( !triggered ) continue;

    effCounter["triggered"] += 1;

    if(candidates->MET < 30) continue;

    effCounter["MET"] += 1;

    afterTriggerAndMETCut.Fill(0., 0., 0., candidates->MET, 0., weight_);

    double massDist = 100.;
    std::vector<TreeCandidate>::const_iterator cand = candidates->candidates_.end();
    std::vector<TreeCandidate>::const_iterator it = candidates->candidates_.begin();
    for( ; it != candidates->candidates_.end(); ++it ) {
      if( it->isStandAloneL || it->isStandAloneH ) continue;
      if( !it->isGlobalMuonL || !it->isGlobalMuonH ) continue;
      if( it->leptonChargeL == it->leptonChargeH ) continue;
      if( fabs(it->leptonEtaL) > 2.4 || fabs(it->leptonEtaH) > 2.4 ) continue;
      if( !(it->triggerMatchL != 0 && it->triggerMatchH != 0 && it->triggerMatchL != it->triggerMatchH) ) continue;
      if( it->leptonPtL < 10. || it->leptonPtH < 20. ) continue;
      if( it->leptonD0L > 0.2 || it->leptonD0H > 0.2 ) continue;
      // Attempt at an isolation cut, not the same as the one used in the WZ analysis
      // if( it->leptonIsoL > 4. || it->leptonIsoH > 4. ) continue;
      // Isolation cut used in the WZ analysis
      int indexL = it->leptonIndexL;
      int indexH = it->leptonIndexH;
      if( ( candidates->leptons_[indexL].trackIso + candidates->leptons_[indexL].ecalIso + candidates->leptons_[indexL].hcalIso )/it->leptonPtL > 0.15 ||
	  ( candidates->leptons_[indexH].trackIso + candidates->leptons_[indexH].ecalIso + candidates->leptons_[indexH].hcalIso )/it->leptonPtH > 0.15 ) continue;

      // Additional requirements on track quality from the WZ analysis
      if( candidates->leptons_[indexL].normChi2 > 10 ) continue;
      if( candidates->leptons_[indexL].numberOfValidTrackerHits < 10 ) continue;
      if( candidates->leptons_[indexL].numberOfValidPixelHits < 1 ) continue;
      if( candidates->leptons_[indexL].numberOfValidMuonHits < 1 ) continue;
      if( candidates->leptons_[indexL].numberOfMatchedStations < 2 ) continue;
      if( candidates->leptons_[indexH].normChi2 > 10 ) continue;
      if( candidates->leptons_[indexH].numberOfValidTrackerHits < 10 ) continue;
      if( candidates->leptons_[indexH].numberOfValidPixelHits < 1 ) continue;
      if( candidates->leptons_[indexH].numberOfValidMuonHits < 1 ) continue;
      if( candidates->leptons_[indexH].numberOfMatchedStations < 2 ) continue;

      // Z Mass window cut
      if( it->corrDileptonMass < 60. || it->corrDileptonMass > 120. ) continue;
      if( fabs(it->corrDileptonMass - 91.) < massDist ) {
	massDist = fabs(it->corrDileptonMass - 91.);
	cand = it;
      }
    }

    if( cand == candidates->candidates_.end() ) continue;

    hZLeptonPtH->Fill(cand->leptonPtH);
    hZLeptonPtL->Fill(cand->leptonPtL);

    effCounter["Candidate"] += 1;
    double leptonPt = 0.;
    std::vector<TreeLepton>::const_iterator lepton = candidates->leptons_.begin();
    std::vector<TreeLepton>::const_iterator selectedLepton = candidates->leptons_.end();
    for( ; lepton != candidates->leptons_.end(); ++lepton ) {
      if( (lepton->index != cand->leptonIndexL) && (lepton->index != cand->leptonIndexH) ) {
	if( !(lepton->isGlobalMuon) ) continue;
	if( lepton->isStandAlone ) continue;
	if( lepton->pt < 10. ) continue;
	if( fabs(lepton->eta) > 2.4 ) continue;
	// Simplified isolation
	// if( fabs(lepton->iso) > 4. ) continue;
	// Tight isolation
	if( ( lepton->trackIso + lepton->ecalIso + lepton->hcalIso )/lepton->pt > 0.10 ) continue;
      // Additional requirements on track quality from the WZ analysis
      if( lepton->normChi2 > 10 ) continue;
      if( lepton->numberOfValidTrackerHits < 10 ) continue;
      if( lepton->numberOfValidPixelHits < 1 ) continue;
      if( lepton->numberOfValidMuonHits < 1 ) continue;
      if( lepton->numberOfMatchedStations < 2 ) continue;
	if( lepton->pt > leptonPt ) {
	  selectedLepton = lepton;
	  leptonPt = lepton->pt;
	}
      }
    }

    hThirdLeptonPt->Fill(leptonPt);

    if( leptonPt < 20 ) continue;
    effCounter["LeptonPt"] += 1;
    if( selectedLepton != candidates->leptons_.end() ) {
      effCounter["Lepton"] += 1;
      TLorentzVector dilepton(convert(cand->ptCorr, cand->etaCorr, cand->phiCorr, cand->corrDileptonMass));
      TLorentzVector singleLepton(convert(selectedLepton->pt, selectedLepton->eta, selectedLepton->phi, 0.105658));
      afterAllCuts.Fill(dilepton.M(), (dilepton+singleLepton).M(), cand->leptonPtL+cand->leptonPtH+leptonPt, candidates->MET,
			sqrt(pow(selectedLepton->vx - cand->vx, 2) + pow(selectedLepton->vy - cand->vy, 2) + pow(selectedLepton->vz - cand->vz, 2)),
			weight_);

    }
    // std::cout << "eta = " << cand->eta << std::endl;
  }
  TString type("_muons");
  if( electrons_ ) type = "_electrons";
  TFile outputFile(dirName_+"weighted"+type+".root", "RECREATE");
  outputFile.cd();

  hZLeptonPtL->Write();
  hZLeptonPtH->Write();
  hThirdLeptonPt->Write();

  beforeAllCuts.Write();
  afterTriggerAndMETCut.Write();
  afterAllCuts.Write();

  std::map<std::string, int>::const_iterator it = effCounter.begin();
  for( ; it != effCounter.end(); ++it ) {
    std::cout << it->first << " = " << it->second << std::endl;
  }

  std::cout << "Total number of events = " << effCounter["total"] << std::endl;
  std::cout << "Efficiency of the cuts. Each step includes all the previous cuts." << std::endl;
  std::cout << "Trigger requirement efficiency = " << effCounter["triggered"]/float(effCounter["total"]) << std::endl;
  std::cout << "MET cut efficiency = " << effCounter["MET"]/float(effCounter["total"]) << std::endl;
  std::cout << "Candidate efficiency = " << effCounter["Candidate"]/float(effCounter["total"]) << std::endl;
  std::cout << "Third lepton pt cut efficiency = " << effCounter["LeptonPt"]/float(effCounter["total"]) << std::endl;
}
