#ifndef TTHMETPLUSJETSANALYZER_CC
#define TTHMETPLUSJETSANALYZER_CC

/**
 * Original Author: Marco De Mattia
 * Creation date: 29/6/2008
 * Mail: demattia@pd.infn.it
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHMEtplusJetsAnalyzer.h"
#include <fstream>
#include <sstream>

using namespace edm;
using namespace std;
using namespace anaobj;

// Small structs used in the analyze method. Defined in a nameless namespace so that they only exist in this file.
namespace {
  struct pairStruct {
    pairStruct( const OfflineJet * JET1, const OfflineJet * JET2 ) {
      //, const double & PROB ) {
      jet1 = JET1;
      jet2 = JET2;
    }
    double prob;
    const OfflineJet * jet1;
    const OfflineJet * jet2;
    /// To sort the structs
    bool operator< ( const pairStruct& b ) const {
      return prob < b.prob;
    }
  };
  struct tripletStruct : public pairStruct {
    tripletStruct( const OfflineJet * JET1, const OfflineJet * JET2, const OfflineJet * JET3 ) : pairStruct(JET1, JET2) {
      jet3 = JET3;
    }
    const OfflineJet * jet3;
    bool operator< ( const tripletStruct& b ) const {
      return prob < b.prob;
    }
  };
}

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
  vtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "vtxLabel" ) ),
  withL1ForwardJets_( iConfig.getUntrackedParameter<bool>("withL1ForwardJets") ),
  vtxAssoc_( iConfig.getUntrackedParameter<bool>("VtxAssoc") ),
  higgsFileName_(iConfig.getUntrackedParameter<string>("HiggsFileName") ),
  hadronicTopFileName_(iConfig.getUntrackedParameter<string>("HadronicTopFileName") ),
  qcdFileName_(iConfig.getUntrackedParameter<string>("QCDfileName") ),
  qcdHistoFileName_(iConfig.getUntrackedParameter<string>("QCDhistoFileName") ),
  jetEtCut_(iConfig.getUntrackedParameter<double>("JetEtCut") ),
  jetEtaCut_(iConfig.getUntrackedParameter<double>("JetEtaCut") ),

  countTTHdecaysFileName_(iConfig.getUntrackedParameter<string>("CountTTHdecaysFileName") ),
  countTTHdecays2tagsFileName_(iConfig.getUntrackedParameter<string>("CountTTHdecays2tagsFileName") ),
 
  outputFileName_(iConfig.getUntrackedParameter<string>("OutputFileName") ),
  tmvaSuffix_(iConfig.getUntrackedParameter<string>("TMVAsuffix") ),
  useTagMatrixForQCD_(iConfig.getUntrackedParameter<bool>("UseTagMatrixForQCD") ),
  eventCounter_(0),
  l1Eff_(0)
{

  countTTHdecays_ = new ttHdecaysCounter(countTTHdecaysFileName_);
  countTTHdecays2tags_ = new ttHdecaysCounter(countTTHdecays2tagsFileName_);

  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  eventVariablesPresel_ = new EventVariables(higgsFileName_, hadronicTopFileName_, qcdFileName_, "presel", outputFile_, true, tmvaSuffix_);

  // Production of pseudo-events for qcd with 2 b-tags.
  if ( useTagMatrixForQCD_ ) {
    qcdbTagMatrixMultiplier_ = new QCDbTagMatrix(higgsFileName_, hadronicTopFileName_, qcdFileName_, "2tags", outputFile_, true, qcdHistoFileName_, 2, tmvaSuffix_);
  }
  else {
    eventVariables2Tags_ = new EventVariables(higgsFileName_, hadronicTopFileName_, qcdFileName_, "2tags_tagMatrix", outputFile_, true, tmvaSuffix_);
  }

  jetVertexAssociator_ = new JetVertexAssociator(jetEtCut_,jetEtaCut_);

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

  // Take the offlineMEt
  edm::Handle<OfflineMEt> offlineMEt;
  iEvent.getByLabel( offlineMEtLabel_, offlineMEt );

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

  //--------------------------
  // Reco vertex collection
  //--------------------------
  Handle<BaseVertexCollection> recVtxs;
  BaseVertexCollection recoVertexes;
  try {
    iEvent.getByLabel( vtxLabel_ , recVtxs);
  } 
  catch (...) {
    std::cerr << "Could not find the recovertex collection" << std::endl;
    return;
  } 
  recoVertexes = *(recVtxs.product());

  // Select only jets with Et>(JetEtCut_)GeV and Eta<(JetEtaCut_) and write them in the goodJets collection
  vector<const OfflineJet *> goodJets;
  vector<const OfflineJet *> goodbTaggedJets;

  // Consider only the first 8 most energetic (Et) jets
  // --------------------------------------------------
  int numGoodJets = 0;

  // Jet-Vertex Algorithm
  // Need to define this outside the if clause, otherwise the objects will not survive it and the
  // pointers stored in goodJets will be invalid.
  OfflineJetCollection associatedJet;
  if(vtxAssoc_){
    associatedJet = jetVertexAssociator_->associate(*(caloJets.product()),recoVertexes);
    for ( OfflineJetCollection::const_iterator assocJetIt = associatedJet.begin(); assocJetIt != associatedJet.end() && numGoodJets < 8; ++assocJetIt) {
      //non metto i tagli in et ed eta sono  fatti internamente dall'algoritmo di associazione
      goodJets.push_back(&(*assocJetIt));
      ++numGoodJets ;

      // -------------------------------------------------------------- //
      // -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
      // -------------------------------------------------------------- //
      // Consider as tagged those jets with highEff > 5.3.
      // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
      // Set the b-tag cut value
      float medium = 5.3;
      if ( assocJetIt->discriminatorHighEff()>medium ) goodbTaggedJets.push_back(&(*assocJetIt));
    }
  }else{
    for ( OfflineJetCollection::const_iterator allJetIt = caloJets->begin(); allJetIt != caloJets->end() && numGoodJets < 8; ++allJetIt ) {
      if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) {
	goodJets.push_back(&(*allJetIt));
        ++numGoodJets;
      	
	// -------------------------------------------------------------- //
	// -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
	// -------------------------------------------------------------- //
	// Consider as tagged those jets with highEff > 5.3.
	// high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
	// Set the b-tag cut value
	float medium = 5.3;
	if ( allJetIt->discriminatorHighEff()>medium ) goodbTaggedJets.push_back(&(*allJetIt));
      }
    }
  }

  if (responsePair.first) {

    // Preselection
    if ( goodJets.size() >= 5 && offlineMEt->corrL3mEtSig() > 3.0 ) {
      // Fill variables after the preselection
      eventVariablesPresel_->fill( goodJets, goodbTaggedJets, &(*offlineMEt) );

      // If it is QCD, do not cut on the b-tags, but loop on all the combinations of jets
      if ( useTagMatrixForQCD_ ) {
        cout << "calling multiply for event = " << eventCounter_ << endl;
        qcdbTagMatrixMultiplier_->multiply( goodJets, &(*offlineMEt) );
      }
      else {

        // Require at least two b-tags
        if ( goodbTaggedJets.size() >= 2 ) {

          //    vector<double> eventVariablesVector(
          eventVariables2Tags_->fill( goodJets, goodbTaggedJets, &(*offlineMEt) );
          //    );

        } // end if at least two b-tags
      }
    } // end preselection cuts
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
  delete eventVariablesPresel_;
  if ( useTagMatrixForQCD_ ) { delete qcdbTagMatrixMultiplier_; }
  else { delete eventVariables2Tags_; }
  outputFile_->Write();
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
