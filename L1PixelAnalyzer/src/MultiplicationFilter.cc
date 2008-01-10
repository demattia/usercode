//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/MultiplicationFilter.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimplePixelJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/Summary.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"

// ModJet class used in the multiplication
#include "AnalysisExamples/AnalysisObjects/interface/ModOfflineJet.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

#include "TRandom.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MultiplicationFilter::MultiplicationFilter(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  cenJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "CenJets" ) ),
  forJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "ForJets" ) ),
  tauJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "TauJets" ) ),
  l1MEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "L1MEt" ) ),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  simplePixelJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimplePixelJets" ) ),
  globalMuonLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons" ) ),
  simpleElectronLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons" ) ),
  simpleTauLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus" ) ),
  summaryLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "Summary" ) ),
  genJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GenJets" ) ),
  minMultiplicationEt_( iConfig.getUntrackedParameter<double>( "MinMultiplicationEt" ) ),
  mEtAlpha_( iConfig.getUntrackedParameter<double>( "MEtAlpha" ) ),
  seed_( iConfig.getUntrackedParameter<int>( "Seed" ) ),
  outputFileName_( iConfig.getUntrackedParameter<string>( "OutputFileName" ) ),
  eventType_( iConfig.getParameter<unsigned int>( "EventType" ) )
{
  //now do what ever initialization is needed

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventCounter_=0;

  // Multiplier, pass the minimum et of jets to be changed and alpha factor for MEt
  multiplier = new Multiplier( 0.3, minMultiplicationEt_, mEtAlpha_ );

  // White background for the canvases
  gROOT->SetStyle("Plain");

  // Count the total number of generated and written events
  nTotChanged_ = 0;
  nTotWritten_ = 0;
  nTotGoodDiscarded_ = 0;

  // Offline
  produces<anaobj::OfflineJetCollection>( offlineJetLabel_.label() );
  produces<anaobj::OfflineMEt>( offlineMEtLabel_.label() );
  // Summary
  produces<anaobj::Summary>( summaryLabel_.label() );
}


MultiplicationFilter::~MultiplicationFilter()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  cout << "Total number of generated events = " << nTotChanged_ << endl;
  cout << "Total number of written events = " << nTotWritten_ << endl;

}

//
// member functions
//

// ------------ method called to for each event  ------------
bool
MultiplicationFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  // Offline Jet and MEt
  // -------------------
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );
  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  // Summary
  // -------
  edm::Handle < Summary > summary;
  iEvent.getByLabel( summaryLabel_, summary );

  eventCounter_++;
  if ( eventCounter_/100 == float(eventCounter_)/100. ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
  }

  // GenJets
  // -------
  edm::Handle < GenJetCollection > genJets;
  iEvent.getByLabel( genJetLabel_, genJets );

  // -------------- //
  // Multiplication //
  // -------------- //

  vector<ModOfflineJet> * modOfflineJetPtr = multiplier->initialize( *caloJets, *caloMET );
  OfflineMEt * offlineMEt = multiplier->offlineMEt();

  multiplier->fillGen( *genJets );

  int numChanged = 0;

  // Change the et of offlineJets associated to genJets
  gRandom->SetSeed(seed_);
  multiplier->generate();
  numChanged = multiplier->numChanged();

  // Proceed only of at least one jet was changed.
  // This means that if no offlineJet was matched with a GenJet
  // the event will not saved.
  int numGoodJets = 0;
  if ( numChanged != 0 ) {

    ++nTotChanged_;

    // New collection of offlineJets to be saved
    auto_ptr<OfflineJetCollection> vecOfflineJetsPtr( new OfflineJetCollection );

    vector<ModOfflineJet>::const_iterator caloJets_it = modOfflineJetPtr->begin();
    for( ; caloJets_it != modOfflineJetPtr->end(); ++caloJets_it ) {
      float jetEt = caloJets_it->et();

      // Save the new jet
      vecOfflineJetsPtr->push_back( OfflineJet( jetEt, caloJets_it->eta(), caloJets_it->phi(), caloJets_it->et(),
                                                caloJets_it->emEnergyFraction(),
                                                caloJets_it->discriminatorHighEff(), caloJets_it->discriminatorHighPur(),
                                                caloJets_it->jetMass(),
                                                caloJets_it->tkNumS1(), caloJets_it->tkSumPtS1(), caloJets_it->tagTkMassS1(),
                                                caloJets_it->tkNumS2(), caloJets_it->tkSumPtS2(), caloJets_it->tagTkMassS2(),
                                                caloJets_it->tkNumS3(), caloJets_it->tkSumPtS3(), caloJets_it->tagTkMassS3() ) );

      // Evalute the number of goodJets
      if ( jetEt > 25. && fabs(caloJets_it->eta()) < 3. ) ++numGoodJets;
    }

    // Store the event only if it has at least 5 jets with Et > 25 GeV and eta < 3.0
    if ( numGoodJets >= 5 ) {
      auto_ptr<OfflineMEt> offlineMEtPtr( new OfflineMEt( offlineMEt->et(), offlineMEt->phi(), offlineMEt->sumEt(), offlineMEt->mEtSig() ) );

      // Summary
      // -------
      int eventType = eventType_;
      // Multiply by 100 the eventType of multiplied events
      auto_ptr<Summary> summaryPtr( new Summary( eventType, summary->jetNum(), summary->hTobb(), summary->semileptonic() ) );

      iEvent.put( vecOfflineJetsPtr, offlineJetLabel_.label() );
      iEvent.put( offlineMEtPtr, offlineMEtLabel_.label() );
      iEvent.put( summaryPtr, summaryLabel_.label() );

    }
  }// end if numChanged != 0

  // If no jets were changed, do not add the event, but evaluate the cuts
  else {
    vector<ModOfflineJet>::const_iterator caloJets_it = modOfflineJetPtr->begin();
    for( ; caloJets_it != modOfflineJetPtr->end(); ++caloJets_it ) {
      if ( caloJets_it->et() > 25. && fabs(caloJets_it->eta()) < 3. ) ++numGoodJets;
    }
    if ( numGoodJets >= 5 ) ++nTotGoodDiscarded_;
  }

  // Store the event only if it has at least 5 jets with Et > 25 GeV and eta < 3.0
  bool write = false;
// Note: num
  if ( numGoodJets >= 5 && numChanged != 0 ) {
    write = true;
    ++nTotWritten_;
  }

  return write;
//  return true;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void MultiplicationFilter::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void MultiplicationFilter::endJob() {
  ofstream outputFile( outputFileName_.c_str() );
  outputFile << "Total events = " << eventCounter_ << endl;
  outputFile << "Generated events = " << nTotChanged_ << endl;
  outputFile << "Writted events = " << nTotWritten_ << endl;
  outputFile << "Events passing the cuts but not regenerated (discarded) = " << nTotGoodDiscarded_ << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MultiplicationFilter);
