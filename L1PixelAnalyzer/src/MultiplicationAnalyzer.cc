// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/MultiplicationAnalyzer.h"

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
#include "AnalysisExamples/AnalysisObjects/interface/ModJet.h"

// Multiplier
#include "AnalysisExamples/L1PixelAnalyzer/interface/Multiplier.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------

L1Trig MultiplicationAnalyzer::L1Trigger;

// Constructors and destructor
// ---------------------------
MultiplicationAnalyzer::MultiplicationAnalyzer(const edm::ParameterSet& iConfig) :
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
  genJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GenJets" ) )
{

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd();

  // White background for the canvases
  // ---------------------------------
  gROOT->SetStyle("Plain");

}

MultiplicationAnalyzer::~MultiplicationAnalyzer()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  OutputFile->Write();
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void MultiplicationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace anaobj;

  // L1 Calo
  // -------
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  // Should we get rid of this try/catch?
  // ------------------------------------
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
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;

  // SimpleElectrons
  // ---------------
  edm::Handle < SimpleElectronCollection > simpleElectrons;

  // SimpleTaus
  // ----------
  edm::Handle < SimpleTauCollection > simpleTaus;

  // Summary
  // -------
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

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }


  // Offline analysis
  // ----------------

  // MEt
  // ---
  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  // HiVariables
  // -----------
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Take genParticleCandidates
  // --------------------------
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  // Take the mothers
  // ----------------
  MCParticleCollection::const_iterator MCp = MCpartons->begin();

  // Pixel jets
  // ----------
  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );
  const SimplePixelJetCollection pixelJets = *(pixelJetsHandle.product());

  // GenJets
  // -------
  edm::Handle < GenJetCollection > genJets;
  iEvent.getByLabel( genJetLabel_, genJets );

  // -------------------
  // end initializations
  // -------------------

  // Constants
  // ---------
  const float pi = 3.1415926;
  const float drmin = 0.3;                   // matching cone jets/partons
  const float drmin2 = 0.09;                 // matching cone jets/partons
  const float drmin2ext = drmin2*sqrt(2);    // outside ring for background subtraction
  const float seed_ass_ptmin = 10;           // minimum pt of clusters associated to jets
  const float jet_ass_ptmin = 10;            // minimum pt of jets associated to clusters
  const float l1jet_ass_ptmin = 10;          // minimum pt of l1 jets associated to jets

  // Create new collections with jets satisfying the thresholds

  // L1Jets (central and tau)
  BaseJetCollection goodL1Jets;
  BaseJetCollection::const_iterator l1Jet_it = l1eCenJets->begin();
  for( ; l1Jet_it != l1eCenJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  for( l1Jet_it = l1eTauJets->begin(); l1Jet_it != l1eTauJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  for( l1Jet_it = l1eForJets->begin(); l1Jet_it != l1eForJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  sort( goodL1Jets.rbegin(), goodL1Jets.rend(), EtSort<BaseJet>() );

  // offlineJets
  OfflineJetCollection goodOfflineJets;
  OfflineJetCollection::const_iterator offlineJet_it = caloJets->begin();
  for( ; offlineJet_it != caloJets->end(); ++offlineJet_it ) {
    //    if ( offlineJet_it->et() > jet_ass_ptmin && fabs(offlineJet_it->eta()) < max_eta ) {
    if ( offlineJet_it->et() > jet_ass_ptmin ) {
      goodOfflineJets.push_back(*offlineJet_it);
    }
  }
  sort( goodOfflineJets.rbegin(), goodOfflineJets.rend(), EtSort<OfflineJet>() );

  // GenJets
  GenJetCollection goodGenJets;
  GenJetCollection::const_iterator genJet_it = genJets->begin();
  for( ; genJet_it != genJets->end(); ++genJet_it ) {
    //    if ( genJet_it->et() > seed_ass_ptmin && fabs(genJet_it->eta()) < max_eta ) {
    if ( genJet_it->et() > seed_ass_ptmin ) {
      goodGenJets.push_back(*genJet_it);
    }
  }
  sort( goodGenJets.rbegin(), goodGenJets.rend(), EtSort<GenJet>() );


  Multiplier multiplier( *l1eCenJets, *l1eTauJets, *l1eForJets, *caloJets );







  // ATTENTION
  // If the vectors passed for the association are altered, the map is
  // invalidated (contains pointers to the elements of the vectors).
  // ---------

  // Associate L1Jets with offlineJets
  AssociatorEt<BaseJet, OfflineJet> l1OffAssociator( drmin );
  auto_ptr<map<const BaseJet*, const OfflineJet*> > l1JetOfflineJetMap( l1OffAssociator.Associate( goodL1Jets, goodOfflineJets ) );

  // Fill the vector of offlineJets associated to L1Jets
  OfflineJetCollection goodAssocOfflineJets;
  map<const BaseJet*, const OfflineJet*>::const_iterator l1OffMap_it = l1JetOfflineJetMap->begin();
  for( ; l1OffMap_it != l1JetOfflineJetMap->end(); ++l1OffMap_it ) {
    goodAssocOfflineJets.push_back( *(l1OffMap_it->second) );
    //    cout << "goodAssocOfflineJet et = " << l1OffMap_it->first->et() << endl;

    // Fill the histograms for goodL1Jet vs goodOfflineJet
    const BaseJet * l1JetPtr = l1OffMap_it->first;
    const OfflineJet * offlineJetPtr = l1OffMap_it->second;
    int ipt = (int)(offlineJetPtr->et()/20.);
    if ( ipt>9 ) ipt=9;
    int ieta = (int)(fabs(offlineJetPtr->eta())/1.2);
    if ( ieta>3 ) ieta=3;
    float offlineJetEt = offlineJetPtr->et();
    float l1JetEt = l1JetPtr->et();
    float l1JetEta = l1JetPtr->eta();
    float deta = offlineJetPtr->eta() - l1JetEta;
    float dphi = pi - fabs(fabs(offlineJetPtr->phi()-l1JetPtr->phi())-pi);
    float dr2 = deta*deta+dphi*dphi;
  }
  sort( goodAssocOfflineJets.begin(), goodAssocOfflineJets.end(), EtSort<OfflineJet>() );

  // Associate this new collection the goodGenJets
  AssociatorEt<GenJet, OfflineJet> genOffAssociator( drmin );
  auto_ptr<map<const GenJet*, const OfflineJet*> > genJetOfflineJetMap( genOffAssociator.Associate( goodGenJets, goodAssocOfflineJets ) );

  // Loop on the map and fill the histograms for goodAssocOfflineJet vs genJet
  map<const GenJet*, const OfflineJet*>::const_iterator genOffMap_it = genJetOfflineJetMap->begin();
  for( ; genOffMap_it != genJetOfflineJetMap->end(); ++genOffMap_it ) {
    const GenJet * genJetPtr = genOffMap_it->first;
    const OfflineJet * offlineJetPtr = genOffMap_it->second;
    int ipt = (int)(genJetPtr->et()/20.);
    if ( ipt>9 ) ipt=9; 
    int ieta = (int)(fabs(genJetPtr->eta())/1.2);
    if ( ieta>3 ) ieta=3;
    float genJetEt = genJetPtr->et();
    float genJetEta = genJetPtr->eta();
    float deta = genJetEta - offlineJetPtr->eta();
    float dphi = pi - fabs(fabs(genJetPtr->phi()-offlineJetPtr->phi())-pi);
    float dr2 = deta*deta+dphi*dphi;
  }


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif

}

// ------------ method called once each job just before starting event loop  ------------
// --------------------------------------------------------------------------------------
void MultiplicationAnalyzer::beginJob(const edm::EventSetup&) {
}


// ------------ method called once each job just after ending the event loop  ------------
// ---------------------------------------------------------------------------------------
void MultiplicationAnalyzer::endJob() {

}

// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(MultiplicationAnalyzer);
