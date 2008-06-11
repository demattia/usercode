//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/OfflineProducer.h"

#include "AnalysisExamples/AnalysisClasses/interface/SimpleJet.h"
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaPhi.h"

// Classes to be stored
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/Summary.h"

// For Tracks
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"

// For the offline jets and corrections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Electron.h"
//#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
// Invariant Mass from Calorimeter
#include "DataFormats/BTauReco/interface/TauMassTagInfo.h"

// For the b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"
//#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

// For the reconstructed vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"

// For file output
#include <fstream>
#include <sstream>

#include <cmath>

#include <memory>

//
// constants, enums and typedefs
//

//
// static data member definitions
//

L1Trig OfflineProducer::L1Trigger;

//
// constructors and destructor
//
OfflineProducer::OfflineProducer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  CaloJetAlgorithm( iConfig.getUntrackedParameter<string>( "CaloJetAlgorithm" ) ),
  L2JetCorrectionService( iConfig.getUntrackedParameter<string>( "L2JetCorrectionService" ) ),
  L3JetCorrectionService( iConfig.getUntrackedParameter<string>( "L3JetCorrectionService" ) ),
  METCollection( iConfig.getUntrackedParameter<string>( "METCollection" ) ),
  L2METCollection( iConfig.getUntrackedParameter<string>( "L2METCollection" ) ),
  L3METCollection( iConfig.getUntrackedParameter<string>( "L3METCollection" ) ),
  genParticleCandidates( iConfig.getUntrackedParameter<string>( "genParticleCandidates" ) ),
  trackCountingHighEffJetTags( iConfig.getUntrackedParameter<string>( "trackCountingHighEffJetTags" ) ),
  trackCountingHighPurJetTags( iConfig.getUntrackedParameter<string>( "trackCountingHighPurJetTags" ) ),
  impactParameterTagInfos( iConfig.getUntrackedParameter<string>( "impactParameterTagInfos" ) ),
  recoVertexLabel_( iConfig.getUntrackedParameter<string>( "RecoVertices" ) ),
  paramGlobalMuons_(iConfig.getUntrackedParameter<edm::InputTag>("ParamGlobalMuons") ),
  electronCandidates_(iConfig.getUntrackedParameter<edm::InputTag>("electronCandidates") ),
  electronHcalIsolation_(iConfig.getUntrackedParameter<edm::InputTag>("electronHcalIsolation") ),
  tauTagInfo_(iConfig.getUntrackedParameter<edm::InputTag>("tauTagInfo") ),
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) ),
  cenJets_( iConfig.getParameter<string>( "CenJets" ) ),
  forJets_( iConfig.getParameter<string>( "ForJets" ) ),
  tauJets_( iConfig.getParameter<string>( "TauJets" ) ),
  l1MEt_( iConfig.getParameter<string>( "L1MEt" ) ),
  offlineJets_( iConfig.getParameter<string>( "OfflineJets" ) ),
  offlineMEt_( iConfig.getParameter<string>( "OfflineMEt" ) ),
  simpleTracks_( iConfig.getParameter<string>( "SimpleTracks" ) ),
  baseVertices_( iConfig.getParameter<string>( "OfflineVertices" ) ),
  MCParticles_( iConfig.getParameter<string>( "MCParticles" ) ),
  globalMuons_( iConfig.getParameter<string>( "GlobalMuons" ) ),
  simpleElectrons_( iConfig.getParameter<string>( "SimpleElectrons" ) ),
  simpleTaus_( iConfig.getParameter<string>( "SimpleTaus" ) ),
  summary_( iConfig.getParameter<string>( "Summary" ) ),
  eventType_( iConfig.getParameter<unsigned int>( "EventType" ) ),
  doL1Trig_( iConfig.getParameter<bool>( "doL1Trig" ) )
{
  //now do what ever initialization is needed

  // White background for the canvases
  gROOT->SetStyle("Plain");

  eventcounter_ = 0;
//  PI_ = 3.141593;
  chargedpi_mass_ = 0.13957018;

  // L1
  if ( doL1Trig_ ) {
    produces<anaobj::BaseJetCollection>( cenJets_ );
    produces<anaobj::BaseJetCollection>( forJets_ );
    produces<anaobj::BaseJetCollection>( tauJets_ );
    produces<anaobj::BaseMEt>( l1MEt_ );
  }
  // Offline
  produces<anaobj::OfflineJetCollection>( offlineJets_ );
  produces<anaobj::OfflineMEt>( offlineMEt_ );
  // SimpleTracks
  produces<anaobj::SimpleTrackCollection>( simpleTracks_ );
  // BaseVertex
  produces<anaobj::BaseVertexCollection>( baseVertices_ );
  // MC
  produces<anaobj::MCParticleCollection>( MCParticles_ );
  // Global muons
  produces<anaobj::GlobalMuonCollection>( globalMuons_ );
  // SimpleElectrons
  produces<anaobj::SimpleElectronCollection>( simpleElectrons_ );
  // SimpleTaus
  produces<anaobj::SimpleTauCollection>( simpleTaus_ );
  // Summary
  produces<anaobj::Summary>( summary_ );
}


OfflineProducer::~OfflineProducer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

//  OutputFile->Write();
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
OfflineProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace l1extra;
  using namespace std;
  using namespace anaobj;

  // L1 Calo
  edm::Handle < L1JetParticleCollection > l1eCenJets;
  edm::Handle < L1JetParticleCollection > l1eForJets;
  edm::Handle < L1JetParticleCollection > l1eTauJets;
  edm::Handle < L1EtMissParticle > l1eEtMiss;

  // Offline
  edm::Handle<reco::CaloMETCollection> caloMET;
  edm::Handle<reco::CaloMETCollection> caloL2MET;
  edm::Handle<reco::CaloMETCollection> caloL3MET;
  edm::Handle<reco::CaloJetCollection> caloJets;

  edm::Handle<reco::JetTagCollection> trackCountingHighEff;
  edm::Handle<reco::JetTagCollection> trackCountingHighPur;
  edm::Handle<reco::TrackIPTagInfoCollection> iPtagInfos;

  // Take genParticleCandidates
  edm::Handle< CandidateCollection > MCpartons;

  // ParamGlobalMuons
  edm::Handle< MuonCollection > paramGlobalMuons;

  // ElectronCandidates
  edm::Handle< PixelMatchGsfElectronCollection > electronCandidates;
  // Get the association vector
  edm::Handle< CandViewDoubleAssociations > electronHcalIsolation;

  // TauTagInfos
  edm::Handle<IsolatedTauTagInfoCollection> tauTagInfoHandle;


  // should get rid of this try/catch?
  try {
    edm::InputTag L1CJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eCentralJetsSource");
    edm::InputTag L1FJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eForwardJetsSource");
    edm::InputTag L1TauJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eTauJetsSource");
    edm::InputTag L1EtMissLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eEtMissSource");

    // Level 1
    if ( doL1Trig_ ) {
      iEvent.getByLabel(L1CJetLabel, l1eCenJets);
      iEvent.getByLabel(L1FJetLabel, l1eForJets);
      iEvent.getByLabel(L1TauJetLabel, l1eTauJets);
      iEvent.getByLabel(L1EtMissLabel, l1eEtMiss);
    }
  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the L1 classes?" << std::endl;
    return;
  }

  try {
    // Offline MEt
    iEvent.getByLabel( METCollection, caloMET );
    iEvent.getByLabel( L2METCollection, caloL2MET );
    iEvent.getByLabel( L3METCollection, caloL3MET );
  }
  catch (...) {
    std::cerr << "Could not find METCollection" << std::endl;
    return;
  }
  try {
    // Offline IC5Jets
    iEvent.getByLabel( CaloJetAlgorithm, caloJets );
  }
  catch (...) {
    std::cerr << "Could not find IC5Jet collection" << std::endl;
    return;
  }
  try {
    // Btags
    iEvent.getByLabel( trackCountingHighEffJetTags, trackCountingHighEff );
    iEvent.getByLabel( trackCountingHighPurJetTags, trackCountingHighPur );
    iEvent.getByLabel( impactParameterTagInfos, iPtagInfos );
  }
  catch (...) {
    std::cerr << "Could not find the b-tagging collection" << std::endl;
      return;
  }
  try {
    // Take genParticleCandidates
    iEvent.getByLabel( genParticleCandidates, MCpartons );
  }
  catch (...) {
    std::cerr << "Could not find genParticleCandidates collection" << std::endl;
    return;
  }
  try {
    // Take the ParamGlobalMuons
    iEvent.getByLabel( paramGlobalMuons_, paramGlobalMuons );
  }
  catch (...) {
    std::cerr << "Could not find the ParamGlobalMuons collection" << std::endl;
  }
  try {
    // Take the ElectronCollection
    iEvent.getByLabel( electronCandidates_, electronCandidates );
    iEvent.getByLabel( electronHcalIsolation_, electronHcalIsolation );
  }
  catch (...) {
    std::cerr << "Could not find the ElectronCandidates colletion" << std::endl;
  }
  try {
    // Take the TauTagInfoCollection
    iEvent.getByLabel( tauTagInfo_, tauTagInfoHandle );
  }
  catch (...) {
    std::cerr << "Could not find the TauTagInfo collection" << std::endl;
  }

  // Reco vertex collection
  Handle<reco::VertexCollection> verticesHandle;
  std::vector<reco::Vertex> recoVertices;
  try {
    iEvent.getByLabel( recoVertexLabel_ , verticesHandle );
    recoVertices = *( verticesHandle.product() );
  }
  catch (...) {
    std::cerr << "Could not find the recovertex collection" << std::endl;
    return;
  }

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // ------- //
  // Offline //
  // ------- //

  // MEt
  // ---

  const reco::CaloMET * MET = &( *(caloMET->begin()) );
  const reco::CaloMET * L2MET = &( *(caloL2MET->begin()) );
  const reco::CaloMET * L3MET = &( *(caloL3MET->begin()) );

  // -----------------------------
  // ATTENTION
  // ---------
  // Not including DPhiMin for now
  // -----------------------------
  //  auto_ptr<OfflineMEt> offlineMEt( new OfflineMEt( MET->et(), MET->phi(), MET->sumEt(), MET->mEtSig(), 0. ) );
  auto_ptr<OfflineMEt> offlineMEt( new OfflineMEt( MET->et(), L2MET->et(), L3MET->et(), MET->phi(), MET->sumEt(), MET->mEtSig() ) );

  iEvent.put( offlineMEt, offlineMEt_ );


  // Offline Jets IC5
  // ----------------

  auto_ptr<OfflineJetCollection> vec_IC5Jets_ptr( new OfflineJetCollection );
  auto_ptr<SimpleTrackCollection> vec_AssocTracks_ptr( new SimpleTrackCollection );

  // Step 0: Declare Ref and RefProd
  OfflineJetRefProd offlineJetRefProd = iEvent.getRefBeforePut<OfflineJetCollection>(offlineJets_);
  SimpleTrackRefProd simpleTrackRefProd = iEvent.getRefBeforePut<SimpleTrackCollection>(simpleTracks_);
  // Declare key_type (counter for the Ref)
  OfflineJetRef::key_type offlineJetId = 0;
  SimpleTrackRef::key_type simpleTrackId = 0;

  int caloJetsSize = caloJets->size();

  // B-tagging
  // ---------
  // Brief description of the btagging algorithm (from the workbook page
  // https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging):
  // -------------------------------------------------------------------
  // "Track Counting" algorithm: This is a very simple tag, exploiting the long lifetime of B hadrons.
  // It calculates the signed impact parameter significance of all good tracks, and orders them by decreasing significance.
  // It's b tag discriminator is defined as the significance of the N'th track, where a good value
  // is N = 2 (for high efficiency) or N = 3 (for high purity).
  // -------------------------------------------------------------------

  int tchenum = trackCountingHighEff->size();
  JetTagCollection::const_iterator bTagHighEff_it = trackCountingHighEff->begin();

  int tchpnum = trackCountingHighPur->size();
  JetTagCollection::const_iterator bTagHighPur_it = trackCountingHighPur->begin();

  int iptkinfonum = iPtagInfos->size();

  // Correct offline jets on the fly
  // const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService, iSetup);

  // Use L2+L3 jet corrections (now default in CMSSW)
  // see https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetAnalysis#JetCorrections
  // ---------------------------------------------------------------------------------
  // Corrections on the fly taken from RecoJets/JetAnalyzers/src/L2L3CorJetsExample.cc
  // from the tag jet_corrections_184.
  // ---------------------------------------------------------------------------------
  const JetCorrector* corrector_L2 = JetCorrector::getJetCorrector (L2JetCorrectionService,iSetup);
  const JetCorrector* corrector_L3 = JetCorrector::getJetCorrector (L3JetCorrectionService,iSetup);  

  //   int taginfocounter = 0;

  // ATTENTION: if the same track is shared between more jets, it will appear more times in the SimpleTrackCollection
  // This can be avoided creating a map to check if the track was already inserted.
  // Map with key = (eta, phi) and value = index in the collection
  map<pair<double, double>, int > noDoublesMap;

  CaloJetCollection::const_iterator caloJets_it = caloJets->begin();
  TrackIPTagInfoCollection::const_iterator TkIpTagInfo_it = iPtagInfos->begin();
  for( ; TkIpTagInfo_it != iPtagInfos->end(); ++TkIpTagInfo_it ) {


    //     ++taginfocounter;
    //    double scale = corrector->correction( *caloJets_it );
    //    double corPt = scale*caloJets_it->pt();
    //    double corEt = scale*caloJets_it->et();

    double scale_L2 = corrector_L2->correction(caloJets_it->p4());
    double scale_L3 = corrector_L3->correction(scale_L2*caloJets_it->p4());
    double corPt = scale_L3*scale_L2*caloJets_it->pt();
    double corEt = scale_L3*scale_L2*caloJets_it->et();

    // Collection of tracks associated to the current caloJet
    SimpleTrackCollection vec_SimpleTracks;

//    int tkNum = TkIpTagInfo_it->selectedTracks().size();
//    int tkNum_S1 = 0;
//    int tkNum_S2 = 0;
//    int tkNum_S3 = 0;

//     int probNum_0 = TkIpTagInfo_it->probabilities(0).size();
//     int probNum_1 = TkIpTagInfo_it->probabilities(1).size();

//    double tkSumPt_S1 = 0.;
//    double tkSumPt_S2 = 0.;
//    double tkSumPt_S3 = 0.;
    const TrackRefVector & vec_TkColl = TkIpTagInfo_it->selectedTracks();
    // Take the IP vector (ordered as the selectedTracks vector)
    const vector<TrackIPTagInfo::TrackIPData> & vec_TkIP = TkIpTagInfo_it->impactParameterData();
    // Additional vectors to store subgroups of selectedTracks
//    TrackRefVector vec_TkColl_S1;
//    TrackRefVector vec_TkColl_S2;
//    TrackRefVector vec_TkColl_S3;


    // Jet mass
    double jetMass = 0;
    if ( corEt > corPt ) jetMass = sqrt( corEt*corEt - corPt*corPt );

    // Evaluate tag tracks invariant mass
    double bTagTkInvMass = tagTracksInvariantMass( vec_TkColl );

    // Fill the offline jets collection
    // --------------------------------

    // Take the em and had energies to evaluate the emfraction
    // See CMSSW/RecoJets/JetAlgorithms/src/JetMaker.cc.
    // Ignoring HF contribution
    double EmEnergyInEB = caloJets_it->emEnergyInEB();
    double EmEnergyInEE = caloJets_it->emEnergyInEE();
    //    double EmEnergyInHF = -(caloJets_it->emEnergyInHF());

    double HadEnergyInEB = caloJets_it->hadEnergyInHB();
    double HadEnergyInEE = caloJets_it->hadEnergyInHE();
    //    double HadEnergyInHF = caloJets_it->hadEnergyInHF();

    //     double EmEnergy = EmEnergyInEB + EmEnergyInEE + EmEnergyInHF;
    //     double HadEnergy = HadEnergyInEB + HadEnergyInEE + HadEnergyInHF;
    double EmEnergy = EmEnergyInEB + EmEnergyInEE;
    double HadEnergy = HadEnergyInEB + HadEnergyInEE;
    double TotalEnergy = EmEnergy + HadEnergy;

    double EmFracFromComp = 0.;
    if ( TotalEnergy != 0. ) EmFracFromComp = EmEnergy/(TotalEnergy);

    // Create the new offlineJet here as will be filled with the references to the associated tracks in the next loop
    OfflineJet offlineJet( corEt, 
                           caloJets_it->eta(), 
                           caloJets_it->phi(), 
                           caloJets_it->et(),
                           EmFracFromComp,
                           //  caloJets_it->emEnergyFraction(),
                           caloJets_it->p4(),
                           caloJets_it->vertex(),
                           bTagHighEff_it->second, 
                           bTagHighPur_it->second,
                           jetMass,
                           bTagTkInvMass );

    vector<TrackIPTagInfo::TrackIPData>::const_iterator TkIP_it = vec_TkIP.begin();
    RefVector<TrackCollection>::const_iterator TkColl_it = vec_TkColl.begin();
    for ( ; TkColl_it != vec_TkColl.end(); ++TkColl_it, ++TkIP_it ) {
      // Fill the SimpleTrackCollection

//      if(Et_Jet>jetEtCut && fabs(Tk_IpS2D) < 3.  && eta_Jet < 3 ){

      // Create a new SimpleTrack and set the reference to the corresponding jet
      double tkEta = (*TkColl_it)->eta();
      double tkPhi = (*TkColl_it)->phi();

      pair<map<pair<double, double>, int >::iterator, bool > element(noDoublesMap.insert( make_pair( make_pair(tkEta, tkPhi), simpleTrackId ) ));
      // First is an iterator to the newly inserted member or the already found one
      // The second entry in the map element is the index of the track in the new collection. It
      // coinciedes with simpleTrackId if the element was newly inserted.

      // Store the reference to the track in the offlineJet
      offlineJet.addTkRef( SimpleTrackRef( simpleTrackRefProd, (element.first)->second ) );

      // If a new insertion was made (no correspondence was found in the map, i.e. new track)
      if( element.second ) {
        SimpleTrack simpleTk( (*TkColl_it)->pt(),
                              tkEta,
                              tkPhi,
                              (*TkColl_it)->dz(),
                              (*TkColl_it)->dzError(),
                              TkIP_it->ip2d.significance() );
        simpleTk.addJetRef( OfflineJetRef( offlineJetRefProd, offlineJetId ) );

        vec_AssocTracks_ptr->push_back( simpleTk );
        // Increase the tkCounter only if the track was effectively inserted
        ++simpleTrackId;
      }
      // if it was already associated, add the reference to this jet
      else {
        (*vec_AssocTracks_ptr)[(element.first)->second].addJetRef( OfflineJetRef( offlineJetRefProd, offlineJetId ) );
      }

//      }

    }

    // Passing corrected Et, eta, phi, uncorrected Et, emEnergyFraction
    vec_IC5Jets_ptr->push_back( offlineJet );
    ++offlineJetId;
    ++caloJets_it;
    ++bTagHighEff_it;
    ++bTagHighPur_it;
  }

  if ( tchenum != tchpnum ) cout << "diff 1" << endl;
  if ( tchenum != iptkinfonum ) cout << "diff 2" << endl;
  if ( tchpnum != iptkinfonum ) cout << "diff 3" << endl;

  if ( tchenum != caloJetsSize ) cout << "diff jets 1" << endl;
  if ( tchpnum != caloJetsSize ) cout << "diff jets 1" << endl;
  if ( iptkinfonum != caloJetsSize ) cout << "diff jets 1" << endl;

  // Write ic5 jets
  iEvent.put( vec_IC5Jets_ptr, offlineJets_ );

  // Write the collection of SimpleTracks associated to the jets
  iEvent.put( vec_AssocTracks_ptr, simpleTracks_ );


  // HepMC
  // -----

  auto_ptr<MCParticleCollection> vec_MC_ptr( new MCParticleCollection );

  int MCnum = 0;
  int WtoNuNum = 0;
  int GoodWtoNuNum = 0;
  int bFromHiggsNum = 0;
  for( CandidateCollection::const_iterator MCp = MCpartons->begin(); MCp != MCpartons->end(); ++MCp ) {
    const Candidate* temp_mother = MCp->mother();
    if ( temp_mother != 0 ) {

      int pid = abs( MCp->pdgId() );
      int Mpid = abs( MCp->mother()->pdgId() );

      // If is a top(6), Z(23), W(24) or Higgs(25) or a daughter of those
      if ( pid  == 6 || pid  == 23 || pid  == 24 || pid  == 25 ||
           Mpid == 6 || Mpid == 23 || Mpid == 24 || Mpid == 25 ) {

        vec_MC_ptr->push_back( MCParticle( MCp->pt(), MCp->eta(), MCp->phi(), MCp->mass(), MCp->pdgId(), MCp->mother()->pdgId() ) );
      }
      // Store also the case in which it is a b with mother != b, or a c with mother != b and c
      else if ( ( pid == 5 && Mpid != 5 ) || ( pid == 4 && Mpid != 5 && Mpid != 4) ) {
        vec_MC_ptr->push_back( MCParticle( MCp->pt(), MCp->eta(), MCp->phi(), MCp->mass(), MCp->pdgId(), MCp->mother()->pdgId() ) );
      }

      // Count the number of nu
      if ( ( Mpid == 24 ) && ( pid ==12 || pid == 14 || pid == 16 ) ) {
        ++WtoNuNum;
        // Define analyzable neutrino decays
        // ---------------------------------
        if ( MCp->et() > 20 ) ++GoodWtoNuNum;
      }

      // Count number of b coming from Higgs
      if ( pid == 5 && Mpid == 25 ) {
        ++bFromHiggsNum;
      }
    }
    ++MCnum;
  } // end loop on MC particles

  iEvent.put( vec_MC_ptr, MCParticles_ );

  bool hTobb = false;
  // Higgs -> bbbar
  if ( bFromHiggsNum == 2 ) {
    hTobb = true;
    cout << "Higgs -> bbbar" << endl;
  }

  bool semileptonic = false;
  // Semileptonic event
  if ( WtoNuNum == 1 ) {
    semileptonic = true;
    cout << "Semileptonic event" << endl;
  }
  bool goodSemileptonic = false;
  // Semileptonic good event
  if ( GoodWtoNuNum == 1 ) {
    goodSemileptonic = true;
    cout << "Semileptonic good event" << endl;
  }

  // ParamGlobalMuons
  // ----------------

  // See here: https://twiki.cern.ch/twiki/bin/view/CMS/MuonIsolation#Suggested_cuts
  // for possible isolation cut values.
  // From which:
  // sumPt < 3 in cone 0.3 will give a 96.7 \pm 0.7 % efficiency with about a factor of 10
  // rejection of muon candidates reconstructed in QCD events.
  // Use chi2/ndof < 10, nHits>=8 as a quality cut on reconstructed muon (applied to a
  // global-muon, this cut allows to suppress fake muons from K/pi). 

  auto_ptr<GlobalMuonCollection> vec_glbmuon_ptr( new GlobalMuonCollection );

  MuonCollection::const_iterator muon_it = paramGlobalMuons->begin();
  std::cout << "paramGlobalMuons->size(): " << paramGlobalMuons->size() << std::endl;
  for ( ; muon_it != paramGlobalMuons->end(); ++muon_it ) {
    std::cout << "inside loop on param muons" << std::endl;
    vec_glbmuon_ptr->push_back( GlobalMuon( muon_it->pt(), 
					    muon_it->eta(), 
					    muon_it->phi(), 
					    muon_it->charge(),
					    muon_it->p4(),
					    muon_it->vertex(),
					    muon_it->getCalEnergy().em,
                                            muon_it->getCalEnergy().had, 
					    muon_it->getCalEnergy().ho,
                                            muon_it->getIsolationR03().sumPt, 
                                            muon_it->getIsolationR03().emEt,
                                            muon_it->getIsolationR03().hadEt,
                                            muon_it->getIsolationR03().hoEt,
                                            muon_it->getIsolationR03().nJets,
                                            muon_it->getIsolationR03().nTracks,
					    muon_it->track().get()->normalizedChi2(),
                                            muon_it->track().get()->numberOfValidHits(),
                                            muon_it->track().get()->dxy(),
                                            muon_it->track().get()->dxyError() ) );
  }

  //  if ( eventcounter_/100 == float(eventcounter_)/100. )
    std::cout << "vec_glbmuon_ptr->size(): " << vec_glbmuon_ptr->size() << std::endl;
  iEvent.put( vec_glbmuon_ptr, globalMuons_ );

  // ElectronCandidates
  // ------------------

  auto_ptr<SimpleElectronCollection> vec_simpelec_ptr( new SimpleElectronCollection );

  CandViewDoubleAssociations::const_iterator elec_it = electronHcalIsolation->begin();
  for ( ; elec_it != electronHcalIsolation->end(); ++elec_it ) {
    PixelMatchGsfElectronRef elec_ref = elec_it->first.castTo<PixelMatchGsfElectronRef>();
    const PixelMatchGsfElectron & elec = *elec_ref;
    vec_simpelec_ptr->push_back( SimpleElectron( elec.pt(), 
						 elec.eta(), 
						 elec.phi(),
						 elec.charge(), 
						 elec.p4(),
						 elec.vertex(),
						 elec.et(), 
						 elec.hadronicOverEm(), 
						 elec_it->second,
						 elec.gsfTrack().get()->dxy(),
						 elec.gsfTrack().get()->dxyError() ) );
  }

  iEvent.put( vec_simpelec_ptr, simpleElectrons_ );

  // Tau candidates
  // --------------

  auto_ptr<SimpleTauCollection> vec_simptau_ptr( new SimpleTauCollection );

  // -------------------------------------------------------------- //
  // Taken from Validation/RecoTau/src/TauTagVal.cc                 //
  // -------------------------------------------------------------- //
  // Other important source to look at:                             //
  // ElectroWeakAnalysis/ZTauTau_DoubleTauJet/src/EWKTauAnalyser.cc //
  // -------------------------------------------------------------- //

  const IsolatedTauTagInfoCollection & tauTagInfo = *(tauTagInfoHandle.product());
  IsolatedTauTagInfoCollection::const_iterator tauTagInfo_it = tauTagInfo.begin();
  for ( ; tauTagInfo_it != tauTagInfo.end(); ++tauTagInfo_it ) {

//     double tauJetE = tauTagInfo_it->jet()->energy();
//     double tauJetP = tauTagInfo_it->jet()->p();
//     double tauJetMass = sqrt( pow(tauJetE,2) - pow(tauJetP,2) );

    //get the tracks from the jetTag
    double tauTagTkInvMass = tagTracksInvariantMass( tauTagInfo_it->allTracks() );
    //get the selected tracks used to compute the isolation
    double tauIsoTkInvMass = tagTracksInvariantMass( tauTagInfo_it->selectedTracks() );

    vec_simptau_ptr->push_back( SimpleTau( tauTagInfo_it->jet()->pt(), tauTagInfo_it->jet()->eta(), tauTagInfo_it->jet()->phi(),
                                           tauTagTkInvMass, tauTagInfo_it->allTracks().size(),
                                           tauIsoTkInvMass, tauTagInfo_it->selectedTracks().size() ) ); 
  }

  iEvent.put( vec_simptau_ptr, simpleTaus_ );

  // RecoVertex
  // ----------

  auto_ptr<BaseVertexCollection> vec_baseVertex_ptr( new BaseVertexCollection );

  std::vector<reco::Vertex>::const_iterator recoVerticesIt = recoVertices.begin();
  for(; recoVerticesIt != recoVertices.end(); ++recoVerticesIt){
    vec_baseVertex_ptr->push_back( BaseVertex( recoVerticesIt->z(), recoVerticesIt->zError() ) );
  }

  iEvent.put( vec_baseVertex_ptr, baseVertices_ );

  // Level 1 trigger
  // ---------------
  if ( doL1Trig_ ) {

    std::auto_ptr<BaseJetCollection> v_CenJets_ptr( new BaseJetCollection );
    std::auto_ptr<BaseJetCollection> v_ForJets_ptr( new BaseJetCollection );
    std::auto_ptr<BaseJetCollection> v_TauJets_ptr( new BaseJetCollection );

    // All the jets together fot the L1Trigger
    vector<SimpleJet> vec_TriggerCenJet;
    vector<SimpleJet> vec_TriggerForJet;
    vector<SimpleJet> vec_TriggerTauJet;
    for ( L1JetParticleCollection::const_iterator tcj = l1eCenJets->begin(); tcj != l1eCenJets->end(); ++tcj ) {
      vec_TriggerCenJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
      v_CenJets_ptr->push_back( BaseJet( tcj->et(), tcj->eta(), tcj->phi() ) );
    }
    for ( L1JetParticleCollection::const_iterator tfj = l1eForJets->begin(); tfj != l1eForJets->end(); ++tfj ) {
      vec_TriggerForJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
      v_ForJets_ptr->push_back( BaseJet( tfj->et(), tfj->eta(), tfj->phi() ) );
    }
    // Tau jets
    for ( L1JetParticleCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
      vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
      v_TauJets_ptr->push_back( BaseJet( ttj->et(), ttj->eta(), ttj->phi() ) );
    }

    // etMiss() method returns et() method (at least for version 1_7_0_pre10)
    std::auto_ptr<BaseMEt> L1MEt( new BaseMEt( l1eEtMiss->etMiss(), l1eEtMiss->phi(), 0. ) );

    // Put the L1 jest collections in the event
    iEvent.put( v_CenJets_ptr, cenJets_ );
    iEvent.put( v_ForJets_ptr, forJets_ );
    iEvent.put( v_TauJets_ptr, tauJets_ );
    try {
      iEvent.put( L1MEt, l1MEt_ );
    }
    catch (...) {
      cout << "put crash" << endl;
    }
  } // end if doL1Trig_

  // Summary
  // -------

  // Set hTobb = true for qcd and "other" events
  if ( eventType_ > 7 ) hTobb = true;
  // Set semileptonic = true for all non Higgs events
  if ( eventType_ != 5 && eventType_ != 6 && eventType_ != 7 ) semileptonic = true;

  // For now fill with the number of all caloJets
  auto_ptr<Summary> summary_ptr( new Summary( eventType_, caloJetsSize, hTobb, semileptonic ) );

  iEvent.put( summary_ptr, summary_ );


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void OfflineProducer::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void OfflineProducer::endJob() {
}

/** Taken from CMSSW/RecoTauTag/RecoTau/src/CaloRecoTauAlgorithm.cc
 * Evaluate jet tracks invariant mass
 */
double OfflineProducer::tagTracksInvariantMass( const TrackRefVector & selectedTracks ) const {
  //      int TksQsum = 0;
  // setting invariant mass of the signal Tracks system
  math::XYZTLorentzVector mySignalTksInvariantMass(0.,0.,0.,0.);
  if( (int)(selectedTracks.size()) != 0 ) {
    int mySignalTks_qsum = 0;
    for ( int i=0; i<(int)selectedTracks.size(); ++i ) {
      mySignalTks_qsum += (selectedTracks)[i]->charge();
      math::XYZTLorentzVector mychargedpicand_fromTk_LorentzVect( (selectedTracks)[i]->momentum().x(),
                                                                  (selectedTracks)[i]->momentum().y(),
                                                                  (selectedTracks)[i]->momentum().z(),
                                                                  sqrt(pow((double)(selectedTracks)[i]->momentum().r(),2) + pow(chargedpi_mass_,2)));
      mySignalTksInvariantMass += mychargedpicand_fromTk_LorentzVect;
    }
    //        TksQsum = mySignalTks_qsum;    
  }
  return mySignalTksInvariantMass.mass();
}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflineProducer);
