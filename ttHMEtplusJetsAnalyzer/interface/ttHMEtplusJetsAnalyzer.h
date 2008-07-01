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

  int eventCounter_;

  // Declare as static so that only one exists, even if more
  // than one TDAna object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;

  int l1Eff_;

  ttHdecaysCounter * countTTHdecays_;
};


#endif // TTHMETPLUSJETSANALYZER_H
