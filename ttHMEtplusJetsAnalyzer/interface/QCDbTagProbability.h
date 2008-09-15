#ifndef QCDBTAGPROBABILITY_H
#define QCDBTAGPROBABILITY_H

/**
 * class QCDbTagProbability QCDbTagProbability.cc AnalysisExamples/ttHMEtplusJetsAnalyzer/src/QCDbTagProbability.cc
 * Package:    QCDbTagProbability
 * Class:      QCDbTagProbability
 *
 * Original Author: Marco De Mattia
 * Creation date: 4/7/2008
 * Mail: demattia@pd.infn.it
 *
 * This class builds the tag-matrix for QCD.
 * It only considers at most 8 good-jets.
 * ATTENTION:
 * it also draws the histograms that can be
 * used for the QCD 'multiplication', but
 * they should be produced separately on
 * each qcd bin and then merged rescaling them
 * with the appropriate cross sections.
 *
 */

// EDM includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

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

class QCDbTagProbability : public edm::EDAnalyzer {
public:
  explicit QCDbTagProbability(const edm::ParameterSet&);
  ~QCDbTagProbability();

private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int eventCounter_;

  edm::ParameterSet conf_;

  edm::InputTag offlineJetLabel_;
  edm::InputTag MCParticleLabel_;
  double jetEtCut_;
  double jetEtaCut_;
  unsigned int maxConsideredJets_;

  string outputProbabilityFileName_;

  TFile * outputFile_;

  TProfile * jetVsMCpEt_;
  TH1F * taggedJetEt_;
  TH1F * taggedJetEta_;
  TH1F * taggedJetS1_;
  TH1F * taggedJetTagMass_;
  TH1F * taggedJetDiscriminatorHighEff_;
  TH1F * notTaggedJetEt_;
  TH1F * notTaggedJetEta_;
  TH1F * notTaggedJetS1_;
  TH1F * notTaggedJetTagMass_;
  TH1F * notTaggedJetDiscriminatorHighEff_;

  // Bins for the probability matrix
  unsigned int etBinNum_;
  double etBinSize_;
  unsigned int etaBinNum_;
  double etaBinSize_;
  unsigned int s1BinNum_;
  double s1BinSize_;

  unsigned int totEvents_;

  // This will be the multidimensional array
  unsigned int *** taggedJet_;
  unsigned int *** notTaggedJet_;
};

#endif // QCDBTAGPROBABILITY_H
