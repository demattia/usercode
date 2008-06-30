#ifndef TTHMETPLUSJETSANALYZER_CC
#define TTHMETPLUSJETSANALYZER_CC

/**
 * Original Author: Marco De Mattia
 * Creation date: 29/6/2008
 * Mail: demattia@pd.infn.it
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHMEtplusJetsAnalyzer.h"

#include <FWCore/Framework/interface/Frameworkfwd.h>

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

L1Trig ttHMEtplusJetsAnalyzer::L1Trigger;

// Constructors and destructor
// ---------------------------
ttHMEtplusJetsAnalyzer::ttHMEtplusJetsAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  cenJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "CenJets" ) ),
  forJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "ForJets" ) ),
  tauJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "TauJets" ) ),
  l1MEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "L1MEt" ) ),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  globalMuonLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons" ) ),
  simpleElectronLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons" ) ),
  simpleTauLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus" ) ),
  summaryLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "Summary" ) ),
  withL1ForwardJets_( iConfig.getUntrackedParameter<bool>("withL1ForwardJets") ),
  eventcounter_(0),
  l1Eff_(0)
{
  countTTHdecays_ = new ttHdecaysCounter("ttHdecays.txt");
}

ttHMEtplusJetsAnalyzer::~ttHMEtplusJetsAnalyzer()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void ttHMEtplusJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  // L1 Calo
  // -------
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  try {
    iEvent.getByLabel(cenJetLabel_, l1eCenJets);
    iEvent.getByLabel(forJetLabel_, l1eForJets);
    iEvent.getByLabel(tauJetLabel_, l1eTauJets);
    iEvent.getByLabel(l1MEtLabel_, l1eEtMiss);
  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  // Global muons
  // SimpleElectrons
  // SimpleTaus
  // Summary
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  edm::Handle < SimpleTauCollection > simpleTaus;
  edm::Handle < Summary > summary;
  try {
    iEvent.getByLabel( globalMuonLabel_, globalMuons );
    iEvent.getByLabel( simpleElectronLabel_, simpleElectrons );
    iEvent.getByLabel( simpleTauLabel_, simpleTaus );
    iEvent.getByLabel( summaryLabel_, summary );
  }
  catch (...) {
    std::cerr << "One of the remaining collections cannot be found" << endl;
  }
  // --------------------------------------------
  // --- end loading handles with collections ---
  // --------------------------------------------

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }
  // ------------------- //
  // L1 Trigger response //
  // ------------------- //

  // Must be fed level1: cen;for;tau;MEt in this order.
  L1Trigger.Fill(*l1eCenJets, *l1eForJets, *l1eTauJets, l1eEtMiss->et());
  // The bool asks to use or not the level 1 forward jets in the trigger.
  // The first bool is multijet, the second met+jet response.
  pair<bool, bool> responsePair(L1Trigger.Response(withL1ForwardJets_));

  // Evaluate level 1 efficiency for multijet and met+jet
  if (responsePair.first || responsePair.second) ++l1Eff_;

  // ----------------------------------- //
  // Determine the type of event from MC //
  // ----------------------------------- //

  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  countTTHdecays_->countDecays(*MCpartons);
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void ttHMEtplusJetsAnalyzer::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void ttHMEtplusJetsAnalyzer::endJob() {

  countTTHdecays_->writeDecays();
  delete countTTHdecays_;
}

#endif // TTHMETPLUSJETSANALYZER_CC

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ttHMEtplusJetsAnalyzer);
