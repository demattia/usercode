#ifndef HADTOPANDHIGGSMASSPROBABILITY_H
#define HADTOPANDHIGGSMASSPROBABILITY_H

/**
 * class HadTopAndHiggsMassProbability HadTopAndHiggsMassProbability.cc AnalysisExamples/ttHMEtplusJetsAnalyzer/src/HadTopAndHiggsMassProbability.cc
 * Package:    ttHMEtplusJetsAnalyzer
 * Class:      ttHMEtplusJetsAnalyzer
 *
 * Original Author: Marco De Mattia
 * Creation date: 30/6/2008
 * Mail: demattia@pd.infn.it
 *
 */

// EDM includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHdecaysCounter.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"

#include "TFile.h"
#include "TProfile.h"

using namespace std;
using namespace anaobj;
using namespace edm;

// ------------

class HadTopAndHiggsMassProbability : public edm::EDAnalyzer {
public:
  explicit HadTopAndHiggsMassProbability(const edm::ParameterSet&);
  ~HadTopAndHiggsMassProbability();

private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int eventCounter_;

  edm::ParameterSet conf_;

  ttHdecaysCounter * countTTHdecays_;

  edm::InputTag offlineJetLabel_;
  edm::InputTag MCParticleLabel_;
  double jetEtCut_;
  double jetEtaCut_;

  TFile * outputFile_;

  TProfile * jetVsMCpEt_;
  TH1F * higgsMassTrue_;
  TH1F * trueHiggsPairEt_;
  TH1F * trueHiggsPairEta_;
  TH1F * trueHiggsPairDR_;
  TH1F * falseHiggsPairEt_;
  TH1F * falseHiggsPairEta_;
  TH1F * falseHiggsPairDR_;
  TH1F * hadronicTopMassTrue_;
  TH1F * trueHadronicTopTripletEt_;
  TH1F * trueHadronicTopTripletEta_;
  TH1F * trueHadronicTopTripletDphiHiggsHadronicTop_;
  TH1F * falseHadronicTopTripletEt_;
  TH1F * falseHadronicTopTripletEta_;
  TH1F * falseHadronicTopTripletDphiHiggsHadronicTop_;

  // Bins for the probability matrix
  unsigned int higgsEtBinNum_;
  double higgsEtBinSize_;
  unsigned int higgsEtaBinNum_;
  double higgsEtaBinSize_;
  unsigned int higgsDRbinNum_;
  double higgsDRbinSize_;
  unsigned int hadronicTopEtBinNum_;
  double hadronicTopEtBinSize_;
  unsigned int hadronicTopEtaBinNum_;
  double hadronicTopEtaBinSize_;
  unsigned int dPhiHiggsHadronicTopBinNum_;
  double dPhiHiggsHadronicTopBinSize_;

  // This will be the multidimensional array
  unsigned int *** trueH_;
  unsigned int *** falseH_;
  unsigned int *** trueHadronicTop_;
  unsigned int *** falseHadronicTop_;
};

#endif // HADTOPANDHIGGSMASSPROBABILITY_H
