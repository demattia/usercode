#ifndef TTHMETPLUSJETSANALYZER_H
#define TTHMETPLUSJETSANALYZER_H

/**
 * class ttHMEtplusJetsAnalyzer ttHMEtplusJetsAnalyzer.cc AnalysisExamples/ttHMEtplusJetsAnalyzer/src/ttHMEtplusJetsAnalyzer.cc
 * Package:    ttHMEtplusJetsAnalyzer
 * Class:      ttHMEtplusJetsAnalyzer
 *
 * Original Author: Marco De Mattia
 * Creation date: 29/6/2008
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
// ------------

// For the level 1 trigger
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHdecaysCounter.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/Summary.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"
#include "AnalysisExamples/AnalysisClasses/interface/SimpleJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/Particle.h"

#include <string>

using namespace std;
using namespace anaobj;
using namespace edm;

// Class declaration
// -----------------
class ttHMEtplusJetsAnalyzer : public edm::EDAnalyzer {
public:
  explicit ttHMEtplusJetsAnalyzer(const edm::ParameterSet&);
  ~ttHMEtplusJetsAnalyzer();
  /// To check if a good muon is found in the event
  bool goodMuon( const GlobalMuonCollection & globalMuons, const OfflineJetCollection & caloJets );
  /// To check if a good electron is found in the event
  bool goodElectron( const SimpleElectronCollection & simpleElectrons, const OfflineJetCollection & caloJets );

private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  /** Used to fill the matrices for the probability of Higgs jet pairs, top jet triplets and also the b-tag probability matrix for qcd
   * takes by reference the *** because it creates them with new inside and assign the pointers to new values
   */
  void fillProbabilityMatrices(const string & probabilityFileName, unsigned int * binNum, double * binSize, unsigned int ***& trueArray, unsigned int ***& falseArray);
  /// Used to evaluate the ratios and select the combinations of jets
  double evalHiggsPairProbability(const Particle<const OfflineJet> & higgsCandidate) const;

  edm::ParameterSet conf_;

  edm::InputTag cenJetLabel_;
  edm::InputTag forJetLabel_;
  edm::InputTag tauJetLabel_;
  edm::InputTag l1MEtLabel_;
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag MCParticleLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag summaryLabel_;
  bool withL1ForwardJets_;
  string higgsFileName_;
  string hadronicTopFileName_;
  string qcdFileName_;
  double jetEtCut_;
  double jetEtaCut_;

  int eventCounter_;

  // Declare as static so that only one exists, even if more
  // than one TDAna object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;

  int l1Eff_;

  ttHdecaysCounter * countTTHdecays_;

  // This will be the multidimensional array
  unsigned int *** trueH_;
  unsigned int *** falseH_;
  unsigned int *** trueHadronicTop_;
  unsigned int *** falseHadronicTop_;
  unsigned int *** taggedJet_;
  unsigned int *** notTaggedJet_;

  unsigned int higgsBinNum_[3];
  double higgsBinSize_[3];
  unsigned int hadronicTopBinNum_[3];
  double hadronicTopBinSize_[3];
  unsigned int qcdBinNum_[3];
  double qcdBinSize_[3];
};


#endif // TTHMETPLUSJETSANALYZER_H
