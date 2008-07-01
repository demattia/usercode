#ifndef TTHMETPLUSJETSANALYZER_CC
#define TTHMETPLUSJETSANALYZER_CC

/**
 * Original Author: Marco De Mattia
 * Creation date: 29/6/2008
 * Mail: demattia@pd.infn.it
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHMEtplusJetsAnalyzer.h"

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
  eventCounter_(0),
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

  eventCounter_++;
  if ( eventCounter_/100 == float(eventCounter_)/100. ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
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
  // The first index is for the Higgs decay, the second for the ttbar.
  // Check ttHdecaysCounter for the decay codes.
  pair<int,int> decayType(countTTHdecays_->countDecays(*MCpartons));
  // ----------------------------------- //

  // ------------------------------------------------ //
  // Check if a good electron or muon is in the event //
  // ------------------------------------------------ //
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );
  //  bool goodElectronFound = goodElectron( *simpleElectrons, *caloJets);
  //  bool goodMuonFound = goodMuon( *globalMuons, *caloJets);
  // ------------------------------------------------ //

  // Associate partons to offlineJets only for H->bb/cc and tt->leptons+4jets decays
  int higgsDecayType = decayType.first;
  int ttDecayType = decayType.second;
  if ( (higgsDecayType == 0 || higgsDecayType == 1) && (ttDecayType == 11 || ttDecayType == 101 || ttDecayType == 1001) ) {
    cout << "ttH->MEt+4Jets decay" << endl;
  }

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

// See if we have a good high-Pt electron or muon
// ----------------------------------------------
bool ttHMEtplusJetsAnalyzer::goodMuon( const GlobalMuonCollection & globalMuons, const OfflineJetCollection & caloJets ) {
  bool muonEvent = false;
  for ( GlobalMuonCollection::const_iterator muon = globalMuons.begin(); 
	muon != globalMuons.end() && !muonEvent; ++muon ) {
    if ( muon->pt()>25 && fabs(muon->eta())<2.5 ) { 
      // See if there are jets closer than 0.5 from this one
      // ---------------------------------------------------
      double dRmin=25.;
      for ( OfflineJetCollection::const_iterator cal = caloJets.begin(); 
	    cal != caloJets.end(); ++cal ) {
	if ( cal->et()>25. ) { 
	  double dR = DeltaR(cal->eta(), cal->phi(), muon->eta(), muon->phi());
	  if ( dR<dRmin ) dRmin=dR;
	}
      }
      if ( dRmin>0.25 ) muonEvent = true;
    }
  }
  return muonEvent;
}
bool ttHMEtplusJetsAnalyzer::goodElectron( const SimpleElectronCollection & simpleElectrons, const OfflineJetCollection & caloJets ) {
  bool elecEvent = false;
  for ( SimpleElectronCollection::const_iterator elec = simpleElectrons.begin(); 
	elec != simpleElectrons.end() && !elecEvent; ++elec ) {
    if ( elec->et()>30 && fabs(elec->eta())<2.5 && elec->hadOverEm()<0.05 ) {
      // See if there are jets closer than 0.5 from this one
      // ---------------------------------------------------
      double dRmin=25.;
      for ( OfflineJetCollection::const_iterator cal = caloJets.begin(); 
	    cal != caloJets.end(); ++cal ) {
	if ( cal->et()>25. && cal->emEnergyFraction()<0.95 ) { 
	  double dR = DeltaR(cal->eta(), cal->phi(), elec->eta(), elec->phi());
	  if ( dR<dRmin ) dRmin=dR;
	}
      }
      if ( dRmin>0.25 ) elecEvent = true;
    }
  }
  return elecEvent;
}


#endif // TTHMETPLUSJETSANALYZER_CC

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ttHMEtplusJetsAnalyzer);
