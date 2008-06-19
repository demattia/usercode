//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/OfflineAnalyzer.h"

// Classes to be accessed
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
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

L1Trig OfflineAnalyzer::L1Trigger;

//
// constructors and destructor
//
OfflineAnalyzer::OfflineAnalyzer(const edm::ParameterSet& iConfig) :
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
  minDz_( iConfig.getUntrackedParameter<double>( "MinDz" ) ),
  maxDz_( iConfig.getUntrackedParameter<double>( "MaxDz" ) ),
  doTrigger_( iConfig.getUntrackedParameter<bool>( "doTrigger" ) ),
  extendedInfo_( iConfig.getUntrackedParameter<bool>( "extendedInfo" ) ),
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) )
{
  //now do what ever initialization is needed

  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,"RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  OutputFile->cd();

  // White background for the canvases
  gROOT->SetStyle("Plain");

  HiVar = new HiVariables( ( conf_.getUntrackedParameter<string>("HiVarName") ).c_str() );

  uncorr_JetPt_IC5_ = new TH1F( "uncorr_JetPt_IC5", "uncorrected JetPt IC5", 100, 0, 200 );
  corr_JetPt_IC5_ = new TH1F( "corr_JetPt_IC5", "corrected JetPt IC5", 100, 0, 200 );
  JetNumber_IC5_ = new TH1F( "JetNumber_IC5", "Number of IC5 jets", 100, 0, 100 );

  // number of correction types
  corrTypeNum_ = 3;

  TString corrTypeName[] = { "uncorr", "corrL2", "corrL3" };
  TString corrTypeTitle[] = { "uncorrected", "L2 corrected", "L3 corrected" };
  TString varName("MEt");

  for (int i=0; i<corrTypeNum_; ++i) {
    MEtPt_.push_back( new TH1F( corrTypeName[i]+varName, corrTypeTitle[i]+" "+varName, 100, 0, 150 ) );
    MEtPhi_.push_back( new TH1F( corrTypeName[i]+varName+"Phi", corrTypeTitle[i]+" "+varName+" phi", 100, -3.15, 3.15 ) );
    MEtSumEt_.push_back( new TH1F( corrTypeName[i]+varName+"SumEt", corrTypeTitle[i]+" "+varName+" sumEt", 100, 0, 2500 ) );
    MEtmEtSig_.push_back( new TH1F( corrTypeName[i]+varName+"mEtSig", corrTypeTitle[i]+" "+varName+"missing Et significance", 100, 0, 10 ) );
  }

  // Generate histograms for the verteces
  // ------------------------------------
  dz_ = 0.1;
  bins_ = 30;

  dzmax_ = dz_ + dz_*(bins_);

  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 100, 0, 3.15 );

  // Trigger efficiency counters
  // Multijet
  Eff_ = 0;
  Eff_et1_ = 0;
  Eff_et2_ = 0;
  Eff_et3_ = 0;
  Eff_et4_ = 0;
  Eff_cen_ = 0;
  Eff_cen_et1_ = 0;
  Eff_cen_et2_ = 0;
  Eff_cen_et3_ = 0;
  Eff_cen_et4_ = 0;
  Eff_tau_ = 0;
  Eff_tau_et1_ = 0;
  Eff_tau_et2_ = 0;
  Eff_tau_et3_ = 0;
  Eff_tau_et4_ = 0;
  Eff_for_ = 0;
  Eff_for_et1_ = 0;
  Eff_for_et2_ = 0;
  Eff_for_et3_ = 0;
  Eff_for_et4_ = 0;
  Eff_nofor_ = 0;
  Eff_nofor_et1_ = 0;
  Eff_nofor_et2_ = 0;
  Eff_nofor_et3_ = 0;
  Eff_nofor_et4_ = 0;
  // MEt+Jet
  Eff_MEtJet_ = 0;
  Eff_MEtJet_cen_ = 0;
  Eff_MEtJet_tau_ = 0;
  Eff_MEtJet_for_ = 0;
  Eff_MEtJet_nofor_ = 0;
  // Tau
  Eff_tautrig_ = 0;
  Eff_tautrig_single_ = 0;
  Eff_tautrig_ditau_ = 0;

  // MultiJet || MEtJet
  Eff_Multi_Or_MEtJet_ = 0;
  Eff_Multi_Or_MEtJet_nofor_ = 0;

  // Offline
  offlineEffMultijet_ = 0;
  offlineEffMEtJet_ = 0;
  offlineEffTauTrig_ = 0;
  offlineEff_Multi_Or_MEtJet_ = 0;
  offlineEff_Multi_Or_MEtJet_nofor_ = 0;

  eventcounter_ = 0;
  //  PI_ = 3.141593;

  numgoodpjeff_ = new int[20];
  numgoodpjeff_3_ = new int[20];
  numgoodpjeff_4_ = new int[20];
  numgoodpjeff_5_ = new int[20];
  numgoodpjeff_6_ = new int[20];

  for (int i=0; i<20; ++i) {
    numgoodpjeff_[i] = 0;
    numgoodpjeff_3_[i] = 0;
    numgoodpjeff_4_[i] = 0;
    numgoodpjeff_5_[i] = 0;
    numgoodpjeff_6_[i] = 0;
  }
}


OfflineAnalyzer::~OfflineAnalyzer()
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
OfflineAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  MEtPt_[0]->Fill( caloMET->et() );
  MEtPhi_[0]->Fill( caloMET->phi() );
  MEtSumEt_[0]->Fill( caloMET->sumEt() );
  MEtmEtSig_[0]->Fill( caloMET->mEtSig() );

  MEtPt_[1]->Fill( caloMET->corrL2et() );
  MEtPhi_[1]->Fill( caloMET->corrL2phi() );
  MEtSumEt_[1]->Fill( caloMET->corrL2sumEt() );
  MEtmEtSig_[1]->Fill( caloMET->corrL2mEtSig() );

  MEtPt_[2]->Fill( caloMET->corrL3et() );
  MEtPhi_[2]->Fill( caloMET->corrL3phi() );
  MEtSumEt_[2]->Fill( caloMET->corrL3sumEt() );
  MEtmEtSig_[2]->Fill( caloMET->corrL3mEtSig() );

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
    int i = 0;
    for( OfflineJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal, ++i ) {
      vec_calojet.push_back( SimpleJet( cal->et(), cal->eta(), cal->phi() ) );
      uncorr_JetPt_IC5_->Fill( cal->uncorrEt() );
      corr_JetPt_IC5_->Fill( cal->et() );

      // Evaluate DeltaPhi between MET and calo-jets
      // Consider only high-luminosity--good-jets
      //       if ( cal->et() >= 40. && fabs( cal->eta() ) < 3.0 ) {
      //         vec_DPhi.push_back( DeltaPhi( MET_phi, cal->phi() ) );
      //       }

      // cout << "L3 jet["<<i<<"] et = " << cal->et() << endl;
      // cout << "L3 jet["<<i<<"] eta = " << cal->eta() << endl;

      if ( cal->et() >= 25. && fabs( cal->eta() ) < 3.0 ) {
        ++goodIc5Jets;
      }

      // cout << "good jets count = " << goodIc5Jets << endl;
    }

    // Offline cut
    // -----------
    // cout << "L3 corrected Missing Et significance = " << caloMET->corrL3mEtSig() << endl;

    if ( goodIc5Jets >= 5 && caloMET->corrL3mEtSig() > 3. ) offline = true;

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
  bool response_cen_et1 = false;
  bool response_cen_et2 = false;
  bool response_cen_et3 = false;
  bool response_cen_et4 = false;
  L1Trigger.Fill( vec_TriggerCenJet );
  response_cen = L1Trigger.Response();
  response_cen_et1 = L1Trigger.Et1Cut();
  response_cen_et2 = L1Trigger.Et2Cut();
  response_cen_et3 = L1Trigger.Et3Cut();
  response_cen_et4 = L1Trigger.Et4Cut();
  // Forward
  bool response_for = false;
  bool response_for_et1 = false;
  bool response_for_et2 = false;
  bool response_for_et3 = false;
  bool response_for_et4 = false;
  L1Trigger.Fill( vec_TriggerForJet );
  response_for = L1Trigger.Response();
  response_for_et1 = L1Trigger.Et1Cut();
  response_for_et2 = L1Trigger.Et2Cut();
  response_for_et3 = L1Trigger.Et3Cut();
  response_for_et4 = L1Trigger.Et4Cut();
  // Tau
  bool response_tau = false;
  bool response_tau_et1 = false;
  bool response_tau_et2 = false;
  bool response_tau_et3 = false;
  bool response_tau_et4 = false;
  L1Trigger.Fill( vec_TriggerTauJet );
  response_tau = L1Trigger.Response();
  response_tau_et1 = L1Trigger.Et1Cut();
  response_tau_et2 = L1Trigger.Et2Cut();
  response_tau_et3 = L1Trigger.Et3Cut();
  response_tau_et4 = L1Trigger.Et4Cut();
  // Full and no-forward
  bool response = ( response_cen || response_tau || response_for );
  bool response_Et1 = ( response_cen_et1 || response_tau_et1 || response_for_et1 );
  bool response_Et2 = ( response_cen_et2 || response_tau_et2 || response_for_et2 );
  bool response_Et3 = ( response_cen_et3 || response_tau_et3 || response_for_et3 );
  bool response_Et4 = ( response_cen_et4 || response_tau_et4 || response_for_et4 );
  bool response_nofor = ( response_cen || response_tau );
  bool response_Et1_nofor = ( response_cen_et1 || response_tau_et1 );
  bool response_Et2_nofor = ( response_cen_et2 || response_tau_et2 );
  bool response_Et3_nofor = ( response_cen_et3 || response_tau_et3 );
  bool response_Et4_nofor = ( response_cen_et4 || response_tau_et4 );

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

  // Tau trigger
  // -----------
  bool response_tautrig = false;
  bool response_tautrig_single = false;
  bool response_tautrig_ditau = false;
  // Already sorted for the previous case
  //   sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  //   reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );

  // Single tau trigger
  if ( vec_TriggerTauJet.size() != 0 ) {
    if ( vec_TriggerTauJet[0].pt() >= 150. ) {
      response_tautrig_single = true;
    }
  }
  // Di-tau trigger
  if ( vec_TriggerTauJet.size() >= 1 ) {
    if ( vec_TriggerTauJet[1].pt() >= 80. ) {
      response_tautrig_ditau = true;
    }
  }
  if ( response_tautrig_single || response_tautrig_ditau ) response_tautrig = true;

#ifdef DEBUG
  std::cout << "L1Trigger response cen = " << response_cen << std::endl;
  std::cout << "L1Trigger response for = " << response_for << std::endl;
  std::cout << "L1Trigger response tau = " << response_tau << std::endl;
  std::cout << "L1Trigger response = " << response << std::endl;

  std::cout << "Number of L1Jets central = " << vec_TriggerJetCen.size() << std::endl;
  std::cout << "Number of L1Jets forward = " << vec_TriggerJetFor.size() << std::endl;
  std::cout << "Number of L1Jets tau = " << vec_TriggerJetTau.size() << std::endl;
#endif

  // Count the trigger efficiency
  // ----------------------------

  // Multijet
  if ( response_cen ) ++Eff_cen_;
  if ( response_cen_et1 ) ++Eff_cen_et1_;
  if ( response_cen_et2 ) ++Eff_cen_et2_;
  if ( response_cen_et3 ) ++Eff_cen_et3_;
  if ( response_cen_et4 ) ++Eff_cen_et4_;

  if ( response_tau ) ++Eff_tau_;
  if ( response_tau_et1 ) ++Eff_tau_et1_;
  if ( response_tau_et2 ) ++Eff_tau_et2_;
  if ( response_tau_et3 ) ++Eff_tau_et3_;
  if ( response_tau_et4 ) ++Eff_tau_et4_;

  if ( response_for ) ++Eff_for_;
  if ( response_for_et1 ) ++Eff_for_et1_;
  if ( response_for_et2 ) ++Eff_for_et2_;
  if ( response_for_et3 ) ++Eff_for_et3_;
  if ( response_for_et4 ) ++Eff_for_et4_;

  if ( response ) ++Eff_;
  if ( response_Et1 ) ++Eff_et1_;
  if ( response_Et2 ) ++Eff_et2_;
  if ( response_Et3 ) ++Eff_et3_;
  if ( response_Et4 ) ++Eff_et4_;

  if ( response_nofor ) ++Eff_nofor_;
  if ( response_Et1_nofor ) ++Eff_nofor_et1_;
  if ( response_Et2_nofor ) ++Eff_nofor_et2_;
  if ( response_Et3_nofor ) ++Eff_nofor_et3_;
  if ( response_Et4_nofor ) ++Eff_nofor_et4_;

  // MEt+jet
  if ( response_MEtJet_cen ) ++Eff_MEtJet_cen_;
  if ( response_MEtJet_tau ) ++Eff_MEtJet_tau_;
  if ( response_MEtJet_for ) ++Eff_MEtJet_for_;
  if ( response_MEtJet ) ++Eff_MEtJet_;
  if ( response_MEtJet_nofor ) ++Eff_MEtJet_nofor_;

  // Tau trigger
  if ( response_tautrig ) ++Eff_tautrig_;
  if ( response_tautrig_single ) ++Eff_tautrig_single_;
  if ( response_tautrig_ditau ) ++Eff_tautrig_ditau_;

  // multijet || MEt+jet
  if ( response || response_MEtJet ) ++Eff_Multi_Or_MEtJet_;
  if ( response_nofor || response_MEtJet_nofor ) ++Eff_Multi_Or_MEtJet_nofor_;

  // Offline
  if ( response && offline ) ++offlineEffMultijet_;
  if ( response_MEtJet && offline ) ++offlineEffMEtJet_;
  if ( response_tautrig && offline ) ++offlineEffTauTrig_;
  if ( (response || response_MEtJet) && offline ) ++offlineEff_Multi_Or_MEtJet_;
  if ( (response_nofor || response_MEtJet_nofor) && offline ) ++offlineEff_Multi_Or_MEtJet_nofor_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void OfflineAnalyzer::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void OfflineAnalyzer::endJob() {

  OutputFile->cd();

  uncorr_JetPt_IC5_->Write();
  corr_JetPt_IC5_->Write();
  JetNumber_IC5_->Write();

  for (int i=0; i<corrTypeNum_; ++i) {
    MEtPt_[i]->Write();
    MEtPhi_[i]->Write();
    MEtSumEt_[i]->Write();
    MEtmEtSig_[i]->Write();
  }

  ofstream Effoutputfile( OutputEffFileName.c_str() );

  Effoutputfile << "Multijet trigger" << endl;
  Effoutputfile << "----------------" << endl;
  Effoutputfile << "Eff multijet = " << float(Eff_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et1 = " << float(Eff_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et2 = " << float(Eff_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et3 = " << float(Eff_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et4 = " << float(Eff_et4_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff multijet no-forward = " <<  float(Eff_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et1 no-forward = " << float(Eff_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et2 no-forward = " << float(Eff_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et3 no-forward = " << float(Eff_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et4 no-forward = " << float(Eff_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff central jet = " << float(Eff_cen_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et1 = " << float(Eff_cen_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et2 = " << float(Eff_cen_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et3 = " << float(Eff_cen_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et4 = " << float(Eff_cen_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff forward jet = " << float(Eff_for_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et1 = " << float(Eff_for_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et2 = " << float(Eff_for_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et3 = " << float(Eff_for_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et4 = " << float(Eff_for_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff tau jet = " << float(Eff_tau_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et1 = " << float(Eff_tau_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et2 = " << float(Eff_tau_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et3 = " << float(Eff_tau_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et4 = " << float(Eff_tau_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "MEt + Jet trigger" << endl;
  Effoutputfile << "-----------------" << endl;
  Effoutputfile << "Eff MEt + Jet = " << float(Eff_MEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + Jet no-forward = " << float(Eff_MEtJet_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + central jet = " << float(Eff_MEtJet_cen_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + tau jet = " << float(Eff_MEtJet_tau_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + forward jet = " << float(Eff_MEtJet_for_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "Tau trigger" << endl;
  Effoutputfile << "-----------" << endl;
  Effoutputfile << "Eff tau trigger = " << float(Eff_tautrig_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff single tau trigger = " << float(Eff_tautrig_single_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff ditau trigger = " << float(Eff_tautrig_ditau_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "MultiJet || MEtJet" << endl;
  Effoutputfile << "------------------" << endl;
  Effoutputfile << "Eff MultiJet || MEtJet = " << float(Eff_Multi_Or_MEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MultiJet || MEtJet no-forward = " << float(Eff_Multi_Or_MEtJet_nofor_)/float(eventcounter_) << endl;

  Effoutputfile << "Offline efficiency" << endl;
  Effoutputfile << "------------------" << endl;
  Effoutputfile << "offline efficiency after multijet trigger = " << float(offlineEffMultijet_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after MEt + Jet trigger = " << float(offlineEffMEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after Tau trigger = " << float(offlineEffTauTrig_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after (MultiJet || MEtJet) trigger = " << float(offlineEff_Multi_Or_MEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after (MultiJet || MEtJet) no-forward trigger = " << float(offlineEff_Multi_Or_MEtJet_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "Total events = " << eventcounter_ << endl;

  Effoutputfile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflineAnalyzer);
