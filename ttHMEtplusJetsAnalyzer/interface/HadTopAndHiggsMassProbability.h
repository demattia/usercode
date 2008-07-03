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

  // Bins for the probability matrix
  unsigned int etBinNum_;
  double etBinSize_;
  unsigned int etaBinNum_;
  double etaBinSize_;
  unsigned int dRbinNum_;
  double dRbinSize_;

  // This will be the multidimensional array
  unsigned int *** trueH_;
  unsigned int *** falseH_;
};

#endif // HADTOPANDHIGGSMASSPROBABILITY_H
