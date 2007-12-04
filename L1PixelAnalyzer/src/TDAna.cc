//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/TDAna.h"

// Classes to be accessed
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

L1Trig TDAna::L1Trigger;

//
// constructors and destructor
//
TDAna::TDAna(const edm::ParameterSet& iConfig) :
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
  numTkCut_( iConfig.getUntrackedParameter<unsigned int>( "TracksMinimumNum_in_PixelJet" ) ),
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) )
{
  //now do what ever initialization is needed

  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  OutputFile->cd();

  // White background for the canvases
  gROOT->SetStyle("Plain");

  HiVar = new HiVariables( ( conf_.getUntrackedParameter<string>("HiVarName") ).c_str() );

  uncorr_JetPt_IC5_ = new TH1F( "uncorr_JetPt_IC5", "uncorrected JetPt IC5", 100, 0, 200 );
  corr_JetPt_IC5_ = new TH1F( "corr_JetPt_IC5", "corrected JetPt IC5", 100, 0, 200 );
  JetNumber_IC5_ = new TH1F( "JetNumber_IC5", "Number of IC5 jets", 100, 0, 100 );

  MEt_CorrIC5_Pt_ = new TH1F( "MEt_CorrIC5_Pt", "MEt corrected with IC5 jets", 100, 0, 150 );
  MEt_CorrIC5_Phi_ = new TH1F( "MEt_CorrIC5_Phi", "MEt Phi corrected with IC5 jets", 100, -3.15, 3.15 );
  MEt_CorrIC5_SumEt_ = new TH1F( "MEt_CorrIC5_SumEt", "SumEt corrected with IC5 jets", 100, 0, 2500 );
  MEt_CorrIC5_mEtSig_ = new TH1F( "MEt_CorrIC5_mEtSig", "MEt significance corrected with IC5 jets", 100, 0, 10 );

  PixelJet_dz_ = new TH1F( "PixelJet_dz", "PixelJet minimum dz", 100, 0, 10 );
  PixelJet_Num_ = new TH1F( "PixelJet_Num", "Number of pixeljets", 30, 0, 30 );
  PixelJet_Track_Num_ = new TH1F( "PixelJet_Track_Num", "Number of pixeltracks in pixeljets", 15, 0, 15 );

  // Generate histograms for the verteces
  // ------------------------------------
  dz_ = 0.04; // cm
  bins_ = 20;

  dzmax_ = dz_ + dz_*(bins_);

  // Multiple histograms, including mean and stacked histograms
  Multi_Vertex_Dz_ = new MultiTH1F( "Vertex_Dz", "Minimum distance between verteces", 100, 0., 10., bins_, dz_, dzmax_, OutputFile );

  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 100, 0, 3.15 );

}


TDAna::~TDAna()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // Draw the histograms
  //  HiVar->Plot();

  OutputFile->Write();
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
TDAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  // L1 Calo
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  // should get rid of this try/catch?
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
  edm::Handle < GlobalMuonCollection > globalMuons;
  // SimpleElectrons
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  // SimpleTaus
  edm::Handle < SimpleTauCollection > simpleTaus;
  // Summary
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

  // ------- //
  // Offline //
  // ------- //

  // MEt
  // ---

  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  //  Associator<CaloMET, CaloJet> associator( 0.5 );
  //  std::auto_ptr<std::map<const reco::CaloMET*, const reco::CaloJet*> > AssocMap( associator.Associate( *caloMET, *caloJets ) );

  // Take the first (and it should be the only one) element in the map and get the eta of the closest jet (in DeltaR)
  //   double JetPhi = (*AssocMap).begin()->second->phi();

  // Do not use DR association. The MEt has no z component. Use DeltaPhi association.

  // HiVariables
  // -----------

  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Count IC5 jets with Et>=30GeV and |eta| < 3.0
  int goodIc5Jets = 0;

  // Offline cuts
  bool offline = false;

  if ( caloJets->size() != 0 ) {

    std::vector<double> vec_DPhi;
    vec_DPhi.reserve(caloJets->size());

    vector<SimpleJet> vec_calojet; 
    vec_calojet.reserve(caloJets->size());
    // Correct offline jets on the fly
    for( OfflineJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal ) {
      vec_calojet.push_back( SimpleJet( cal->et(), cal->eta(), cal->phi() ) );
      uncorr_JetPt_IC5_->Fill( cal->uncorrEt() );
      corr_JetPt_IC5_->Fill( cal->et() );

      // Evaluate DeltaPhi between MET and calo-jets
      // Consider only high-luminosity--good-jets
      //       if ( cal->et() >= 40. && fabs( cal->eta() ) < 3.0 ) {
      //         vec_DPhi.push_back( DeltaPhi( MET_phi, cal->phi() ) );
      //       }
      if ( cal->et() >= 30. && fabs( cal->eta() ) < 3.0 ) {
        ++goodIc5Jets;
      }
    }

    // Offline cut
    // -----------
    if ( goodIc5Jets >= 5 ) offline = true;


    JetNumber_IC5_->Fill( vec_calojet.size() );

    HiVar->Fill( vec_calojet );

    // Minimum DeltaPhi
    std::sort( vec_DPhi.begin(), vec_DPhi.end() );
    DPhimin_->Fill( *(vec_DPhi.begin()) );

  }
  else {
    std::cout << "ATTENTION: Jet collection empty" << std::endl;
  }

  // Take genParticleCandidates
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  // Take the mothers
  MCParticleCollection::const_iterator MCp = MCpartons->begin();
  for ( ; MCp != MCpartons->end(); ++MCp ) {

#ifdef DEBUG
    std::cout << "For parton number = " << counter << std::endl;
    std::cout << "status = " << MCp->status() << std::endl;
    std::cout << "pdgId = " << MCp->pdgId() << std::endl;
    std::cout << "Et = " << MCp->et() << std::endl;
    std::cout << "Eta = " << MCp->eta() << std::endl;
    std::cout << "Phi = " << MCp->phi() << std::endl;
    std::cout << "Number of mothers = " << MCp->numberOfMothers() << std::endl;
    std::cout << "first mother = " << MCp->mother() << std::endl;
    std::cout << "Mother pdgId = " << MCp->mother()->pdgId() << std::endl;
#endif // DEBUG
  }

  // Pixel jets
  // ----------

  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );

  const SimplePixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // Count the number of pixeljets with at least numTkCut_ tracks
  SimplePixelJetCollection::const_iterator pj_it = pixeljets.begin();
  int goodPjNum = 0;
  for ( ; pj_it != pixeljets.end(); ++pj_it ) {
    if ( pj_it->tkNum() >= int(numTkCut_) ) ++goodPjNum;
  }

  // Pixel trigger requiring at least 2 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)

#ifdef DEBUG
    std::cout << "numTkCut = " << numTkCut_ << std::endl;
#endif

    L1PixelTrig<SimplePixelJet> PJtrig(2, 0.4, numTkCut_);

    PJtrig.Fill( pixeljets );

#ifdef DEBUG
    std::cout << "Pixel trigger response = " << PJtrig.Response() << std::endl;
#endif

    // Level 1 trigger
    // ---------------
    // All the jets together fot the L1Trigger
    vector<SimpleJet> vec_TriggerCenJet;
    vector<SimpleJet> vec_TriggerForJet;
    vector<SimpleJet> vec_TriggerTauJet;
    for ( BaseJetCollection::const_iterator tcj = l1eCenJets->begin(); tcj != l1eCenJets->end(); ++tcj ) {
      vec_TriggerCenJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
    }
    int fjcount = 0;
    for ( BaseJetCollection::const_iterator tfj = l1eForJets->begin(); tfj != l1eForJets->end(); ++tfj ) {
      vec_TriggerForJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
#ifdef DEBUG
      std::cout << "ForwardJet Et["<<fjcount<<"] = " << tfj->et() << std::endl;
      std::cout << "ForwardJet Eta["<<fjcount<<"] = " << tfj->eta() << std::endl;
      std::cout << "ForwardJet Phi["<<fjcount<<"] = " << tfj->phi() << std::endl;
#endif
      ++fjcount;
    }
    // Tau jets
    for ( BaseJetCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
      vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
    }

    // Multijet
    // --------
    // Central
    bool response_cen = false;
    L1Trigger.Fill( vec_TriggerCenJet );
    response_cen = L1Trigger.Response();
    // Forward
    bool response_for = false;
    L1Trigger.Fill( vec_TriggerForJet );
    response_for = L1Trigger.Response();
    // Tau
    bool response_tau = false;
    L1Trigger.Fill( vec_TriggerTauJet );
    response_tau = L1Trigger.Response();
    // Full and no-forward
    bool response = ( response_cen || response_tau || response_for );
    bool response_nofor = ( response_cen || response_tau );

    // MEt + Jet
    // ---------
    // Central
    bool response_MEtJet_cen = false;
    sort( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
    reverse( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
    if ( vec_TriggerCenJet.size() != 0 ) {
      if ( (vec_TriggerCenJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
	response_MEtJet_cen = true;
      }
    }
    // Tau
    bool response_MEtJet_tau = false;
    sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
    reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
    if ( vec_TriggerTauJet.size() != 0 ) {
      if ( (vec_TriggerTauJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
	response_MEtJet_tau = true;
      }
    }
    // Forward
    bool response_MEtJet_for = false;
    sort( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
    reverse( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
    if ( vec_TriggerForJet.size() != 0 ) {
      if ( (vec_TriggerForJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
	response_MEtJet_for = true;
      }
    }
    // Full and no-forward
    bool response_MEtJet = ( response_MEtJet_cen || response_MEtJet_tau || response_MEtJet_for );
    bool response_MEtJet_nofor = ( response_MEtJet_cen || response_MEtJet_tau );

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif
  }

  // ------------ method called once each job just before starting event loop  ------------
  void TDAna::beginJob(const edm::EventSetup&) {
  }

  // void TDAna::meanHistoSetup( TH1F * meanhisto_ptr, vector<TH1F *> vec_histo ) {
  //   meanhisto_ptr->SetBinContent( numdz+1, mean[numdz]->GetMean() );
  //   meanhisto_ptr->SetBinError( numdz+1, Vertex_Dz_[numdz]->GetMeanError() );
  //   meanhisto_ptr->GetXaxis()->CenterLabels();
  //   meanhisto_ptr->GetXaxis()->SetNdivisions(Vertex_Dz_Mean_->GetSize()-2, false);
  // }

  // ------------ method called once each job just after ending the event loop  ------------
  void TDAna::endJob() {

    Multi_Vertex_Dz_->Write();


  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(TDAna);
