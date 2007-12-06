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
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"

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
  eventcounter_=0;
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


  dR_    = 0.05;
  dRmax_ = dR_ + dR_*(bins_);
  // multiple b tag discriminator histograms
  Multi_bDiscrHighEff   = new MultiTH1F( "bDiscrHighEff",   "high efficiency discriminator for b jet",   200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  Multi_bDiscrHighPur   = new MultiTH1F( "bDiscrHighPur",   "high purity discriminator for b jet",       200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  Multi_nobDiscrHighEff = new MultiTH1F( "nobDiscrHighEff", "high efficiency discriminator for no b jet",200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  Multi_nobDiscrHighPur = new MultiTH1F( "nobDiscrHighPur", "high purity discriminator for no b jet",    200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  Multi_tagMassS1       = new MultiTH1F( "tagMassS1",       "tag mass S1 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_tagMassS2       = new MultiTH1F( "tagMassS2",       "tag mass S2 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_tagMassS3       = new MultiTH1F( "tagMassS3",       "tag mass S3 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_nobtagMassS1    = new MultiTH1F( "nobtagMassS1",    "tag mass S1for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_nobtagMassS2    = new MultiTH1F( "nobtagMassS2",    "tag mass S2for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_nobtagMassS3    = new MultiTH1F( "nobtagMassS3",    "tag mass S3for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  Multi_jetEtVSbPt      = new MultiTProfile("jetEtVSbPt",   "jet Et as function as b parton Pt",   500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );
  Multi_jetEtVSnobPt    = new MultiTProfile("jetEtVSnobPt", "jet Et as function as no b parton Pt",500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );

  bHighEff   = new TH1F("bHighEff",  "high efficiency discriminator for b jet",   300,-50.,100.);
  bHighPur   = new TH1F("bHighPur",  "high purity discriminator for b jet",       300,-50.,100.);
  nobHighEff = new TH1F("nobHighEff","high efficiency discriminator for no b jet",300,-50.,100.);
  nobHighPur = new TH1F("nobHighPur","high purity discriminator for no b jet",    300,-50.,100.);

deltaR                        = new TH1F("deltaR","deltaR between matched jet and b-parton",250,0.,1.0);
jetEtVSbParton                = new TH2F("jetEtVSbParton", "jetEtVSbParton", 500,0.,500.,500,0.,500. );
jetUncorrEtVSbParton          = new TH2F("jetUncorrEtVSbParton", "jetUncorrEtVSbParton", 500,0.,500.,500,0.,500. );
tagTkMassS1                   = new TH1F("tagTkMassS1","tagTkMassS1",150,0.,15.);                   
tagTkMassS2                   = new TH1F("tagTkMassS2","tagTkMassS2",150,0.,15.);
tagTkMassS3                   = new TH1F("tagTkMassS3","tagTkMassS3",150,0.,15.);

mPdiVsDiscHighEff             = new TH2F("mPid", "mPidVsDiscHighEff", 34,-7.,26.,200,-100.,100. );
uncorrEtVsDiscHighEff         = new TH2F("uncorrEt", "uncorrEtVsDiscHighEff", 300,0.,300.,200,-100.,100. );
emEnergyFractionVsDiscHighEff = new TH2F("emEnergyFraction", "emEnergyFractionVSDiscHighEff",10,0.,1.,200,-100.,100. );
jetMassVsDiscHighEff          = new TH2F("jetMass", "jetMassVSDiscHighEff",100,0.,100.,200,-100.,100.);
tkNumS1VsDiscHighEff          = new TH2F("tkNumS1", "tkNumS1VSDiscHighEff",30,0.,30.,200,-100.,100.);
tkSumPtS1VsDiscHighEff        = new TH2F("tkSumPtS1","tkSumPtS1VSDiscHighEff",650,0.,650.,200,-100.,100.);
tagTkMassS1VsDiscHighEff      = new TH2F("tagTkMassS1VS","tagTkMassS1VSDiscHighEff",150,0.,15.,200,-100.,100.);
tkNumS2VsDiscHighEff          = new TH2F("tkNumS2", "tkNumS2VSDiscHighEff",30,0.,30.,200,-100.,100.);
tkSumPtS2VsDiscHighEff        = new TH2F("tkSumPtS2","tkSumPtS2VSDiscHighEff",650,0.,650.,200,-100.,100.);
tagTkMassS2VsDiscHighEff      = new TH2F("tagTkMassS2VS","tagTkMassS2VSDiscHighEff",150,0.,15.,200,-100.,100.);
tkNumS3VsDiscHighEff          = new TH2F("tkNumS3", "tkNumS3VSDiscHighEff",30,0.,30.,200,-100.,100.);
tkSumPtS3VsDiscHighEff        = new TH2F("tkSumPtS3","tkSumPtS3VSDiscHighEff",650,0.,650.,200,-100.,100.);
tagTkMassS3VsDiscHighEff      = new TH2F("tagTkMassS3VS","tagTkMassS3VSDiscHighEff",150,0.,15.,200,-100.,100.);


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

  //  AssociatorEt<CaloMET, CaloJet> associator( 0.5 );
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

  vector<OfflineJet> goodIc5JetVec;

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
	goodIc5JetVec.push_back(*cal);
      //         vec_DPhi.push_back( DeltaPhi( MET_phi, cal->phi() ) );
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

//  int Zdaughters=0;
//  int Zcounter=0;
//  vector<int> HiggsDau;
  for ( ; MCp != MCpartons->end(); ++MCp ) {
//    if(MCp->mPid()==25) HiggsDau.push_back(MCp->pid());
//    if(MCp->mPid()==23) ++Zdaughters;
//    if(MCp->pid()==23)  ++Zcounter;

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
  // CHECK for Higgs BR: Higgs inclusive
//  if(Zcounter!=0 || Zdaughters !=0) {
//    cout << "number of Z mother: " << Zcounter << " <--> number of Z daughters: " << Zdaughters << endl;
//    int iter=0;
//    for(vector<int>::iterator HiggsDau_it = HiggsDau.begin(); HiggsDau_it != HiggsDau.end(); ++HiggsDau_it, ++iter)
//      cout << "HiggsDaughters[" << iter << "]: " << *HiggsDau_it << endl;
//  }


  // Pixel jets
  // ----------

  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );

  const SimplePixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // Count the number of pixeljets with at leas numTkCut_ tracks
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

    // if #b-tag
    for( int numdR = 0; numdR < bins_; ++numdR){
      double dRcut = dR_ + dR_*(numdR);
      AssociatorEt<OfflineJet, MCParticle > associator( dRcut );      
//      std::auto_ptr<std::map<const OfflineJet*, 
//	                     const MCParticle*> > assocMap( associator.Associate( *caloJets, *MCpartons ) );
      std::auto_ptr<std::map<const OfflineJet*, 
	                     const MCParticle*> > assocMap( associator.Associate( goodIc5JetVec, *MCpartons ) );
      
      std::map<const OfflineJet*, 
	       const MCParticle*>::const_iterator assocMap_it = assocMap->begin();
      for ( ; assocMap_it != assocMap->end(); ++assocMap_it ) {
	float discHighPur = assocMap_it->first->discriminatorHighPur();
	float discHighEff = assocMap_it->first->discriminatorHighEff();
	float tagMassS1   = assocMap_it->first->tagTkMassS1();
	float tagMassS2   = assocMap_it->first->tagTkMassS2();
	float tagMassS3   = assocMap_it->first->tagTkMassS3();
	float jetEt       = assocMap_it->first->et();
	float partonPt    = assocMap_it->second->pt();

	if( fabs( assocMap_it->second->pid() ) == 5 ){
	  Multi_bDiscrHighEff->Fill(discHighEff,numdR);
	  Multi_bDiscrHighPur->Fill(discHighPur,numdR);  
	  if(tagMassS1 > 0.2) Multi_tagMassS1->Fill(tagMassS1,numdR);
	  if(tagMassS2 > 0.2) Multi_tagMassS2->Fill(tagMassS2,numdR);
	  if(tagMassS3 > 0.2) Multi_tagMassS3->Fill(tagMassS3,numdR);
	  Multi_jetEtVSbPt->Fill(partonPt,jetEt,numdR);
	  if(dRcut == 1.0){
	    double jetEta = assocMap_it->first->eta();
	    double jetPhi = assocMap_it->first->phi();
	    double bEta   = assocMap_it->second->eta();
	    double bPhi   = assocMap_it->second->phi();
	    double dR     = DeltaR(jetEta,jetPhi,bEta,bPhi);
	    deltaR->Fill(dR);
	  }
	}else{
	  Multi_nobDiscrHighEff->Fill(discHighEff,numdR);
	  Multi_nobDiscrHighPur->Fill(discHighPur,numdR);
//	  if(tagMassS1 > 0.2) Multi_nobtagMassS1->Fill(tagMassS1,numdR);
//	  if(tagMassS2 > 0.2) Multi_nobtagMassS2->Fill(tagMassS2,numdR);
//	  if(tagMassS3 > 0.2) Multi_nobtagMassS3->Fill(tagMassS3,numdR);
	  Multi_nobtagMassS1->Fill(tagMassS1,numdR);
	  Multi_nobtagMassS2->Fill(tagMassS2,numdR);
	  Multi_nobtagMassS3->Fill(tagMassS3,numdR);
	  Multi_jetEtVSnobPt->Fill(partonPt,jetEt,numdR);
	}
      }
    }

    AssociatorEt<OfflineJet, MCParticle > associator( 0.3 );
    std::auto_ptr<std::map<const OfflineJet*, 
                           const MCParticle*> > assocMap( associator.Associate( *caloJets, *MCpartons ) );

    std::map<const OfflineJet*, 
             const MCParticle*>::const_iterator assocMap_it = assocMap->begin();
    for ( ; assocMap_it != assocMap->end(); ++assocMap_it ) {
      if( fabs( assocMap_it->second->pid() ) == 5 ){

	jetEtVSbParton->Fill(      assocMap_it->second->pt(),assocMap_it->first->et());
	jetUncorrEtVSbParton->Fill(assocMap_it->second->pt(),assocMap_it->first->uncorrEt());

        float tagMassS1 = assocMap_it->first->tagTkMassS1();
        float tagMassS2 = assocMap_it->first->tagTkMassS2();
        float tagMassS3 = assocMap_it->first->tagTkMassS3();
	if(tagMassS1 > 0.2) tagTkMassS1->Fill(tagMassS1);
	if(tagMassS2 > 0.2) tagTkMassS2->Fill(tagMassS2);
	if(tagMassS3 > 0.2) tagTkMassS3->Fill(tagMassS3);

	float discHighPur = assocMap_it->first->discriminatorHighPur();
	float discHighEff = assocMap_it->first->discriminatorHighEff();
	mPdiVsDiscHighEff->Fill(                          assocMap_it->second->mPid(),            discHighEff );
        uncorrEtVsDiscHighEff->Fill(                      assocMap_it->first->uncorrEt(),         discHighEff );
        emEnergyFractionVsDiscHighEff->Fill(              assocMap_it->first->emEnergyFraction(), discHighEff );
	jetMassVsDiscHighEff->Fill(                       assocMap_it->first->jetMass(),          discHighEff );
	tkNumS1VsDiscHighEff->Fill(                       assocMap_it->first->tkNumS1(),          discHighEff );
	tkSumPtS1VsDiscHighEff->Fill(                     assocMap_it->first->tkSumPtS1(),        discHighEff );
	if(tagMassS1 != 0) tagTkMassS1VsDiscHighEff->Fill(assocMap_it->first->tagTkMassS1(),      discHighEff );
	tkNumS2VsDiscHighEff->Fill(                       assocMap_it->first->tkNumS2(),          discHighEff );
	tkSumPtS2VsDiscHighEff->Fill(                     assocMap_it->first->tkSumPtS2(),        discHighEff );
	if(tagMassS2 != 0) tagTkMassS2VsDiscHighEff->Fill(assocMap_it->first->tagTkMassS2(),      discHighEff );
	tkNumS3VsDiscHighEff->Fill(                       assocMap_it->first->tkNumS3(),          discHighEff );
	tkSumPtS3VsDiscHighEff->Fill(                     assocMap_it->first->tkSumPtS3(),        discHighEff );
	if(tagMassS3 != 0) tagTkMassS3VsDiscHighEff->Fill(assocMap_it->first->tagTkMassS3(),      discHighEff );

      } else{
	nobHighEff->Fill( assocMap_it->first->discriminatorHighEff() );
	nobHighPur->Fill( assocMap_it->first->discriminatorHighPur() );
      }
   }

  // Take the first (and it should be the only one) element in the map and 
  // get the eta of the closest jet (in DeltaR)
  //   double JetPhi = (*AssocMap).begin()->second->phi();

  // Do not use DR association. The MEt has no z component. Use DeltaPhi association.



   // if goodjet >=4,>=5,>=6,>=7  goodIc5jets
 
   // if #b-tag>=2,>=3

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
    
    vector<TH1F*> bDiscrHighEffMultiHistos   = Multi_bDiscrHighEff->multiHistos();
    vector<TH1F*> bDiscrHighPurMultiHistos   = Multi_bDiscrHighPur->multiHistos();
    vector<TH1F*> nobDiscrHighEffMultiHistos = Multi_nobDiscrHighEff->multiHistos();
    vector<TH1F*> nobDiscrHighPurMultiHistos = Multi_nobDiscrHighPur->multiHistos();
    vector<TH1F*> tagMassS1MultiHistos       = Multi_tagMassS1->multiHistos();
    vector<TH1F*> tagMassS2MultiHistos       = Multi_tagMassS2->multiHistos();
    vector<TH1F*> tagMassS3MultiHistos       = Multi_tagMassS3->multiHistos();
    vector<TH1F*> nobtagMassS1MultiHistos    = Multi_nobtagMassS1->multiHistos();
    vector<TH1F*> nobtagMassS2MultiHistos    = Multi_nobtagMassS2->multiHistos();
    vector<TH1F*> nobtagMassS3MultiHistos    = Multi_nobtagMassS3->multiHistos();
    vector<TProfile*> jetEtVSbPtMultiProfiles = Multi_jetEtVSbPt->multiProfiles();
    vector<TProfile*> jetEtVSnobPtMultiProfiles = Multi_jetEtVSnobPt->multiProfiles();

    vector<TH1F*>::iterator bDiscrHighEffMultiHistos_it   = bDiscrHighEffMultiHistos.begin();
    vector<TH1F*>::iterator nobDiscrHighEffMultiHistos_it = nobDiscrHighEffMultiHistos.begin();
    vector<TH1F*>::iterator bDiscrHighPurMultiHistos_it   = bDiscrHighPurMultiHistos.begin();
    vector<TH1F*>::iterator nobDiscrHighPurMultiHistos_it = nobDiscrHighPurMultiHistos.begin();
    vector<TH1F*>::iterator tagMassS1MultiHistos_it       = tagMassS1MultiHistos.begin();
    vector<TH1F*>::iterator tagMassS2MultiHistos_it       = tagMassS2MultiHistos.begin();
    vector<TH1F*>::iterator tagMassS3MultiHistos_it       = tagMassS3MultiHistos.begin();
    vector<TH1F*>::iterator nobtagMassS1MultiHistos_it    = nobtagMassS1MultiHistos.begin();
    vector<TH1F*>::iterator nobtagMassS2MultiHistos_it    = nobtagMassS2MultiHistos.begin();
    vector<TH1F*>::iterator nobtagMassS3MultiHistos_it    = nobtagMassS3MultiHistos.begin();
    vector<TProfile*>::iterator jetEtVSbPtMultiProfiles_it   = jetEtVSbPtMultiProfiles.begin();
    vector<TProfile*>::iterator jetEtVSnobPtMultiProfiles_it = jetEtVSnobPtMultiProfiles.begin();
    TCanvas* canvasEff1        = new TCanvas( "discrHighEffVsDeltaRCanvas_1", "b tag discriminator in high efficiency VS delta R", 1000, 800 );
    canvasEff1->Divide(2,2);
    TCanvas* canvasEff2        = new TCanvas( "discrHighEffVsDeltaRCanvas_2", "b tag discriminator in high efficiency VS delta R", 1000, 800 );
    canvasEff2->Divide(2,2);
    TCanvas* canvasPur1        = new TCanvas( "discrHighPurVsDeltaRCanvas_1", "b tag discriminator in high purity VS delta R", 1000, 800 );
    canvasPur1->Divide(2,2);
    TCanvas* canvasPur2        = new TCanvas( "discrHighPurVsDeltaRCanvas_2", "b tag discriminator in high purity VS delta R", 1000, 800 );
    canvasPur2->Divide(2,2);
    TCanvas* canvasTagMass     = new TCanvas( "tagMassVsDeltaRCanvas", "tag mass VS delta R", 1000, 800 );
    canvasTagMass->Divide(2,2);
    TCanvas* canvasTagMassS1   = new TCanvas( "tagMassS1VsDeltaRCanvas", "tag mass S1 VS delta R", 1000, 800 );
    canvasTagMassS1->Divide(2,2);
    TCanvas* canvasTagMassS2   = new TCanvas( "tagMassS2VsDeltaRCanvas", "tag mass S2 VS delta R", 1000, 800 );
    canvasTagMassS2->Divide(2,2);
    TCanvas* canvasTagMassS3   = new TCanvas( "tagMassS3VsDeltaRCanvas", "tag mass S3 VS delta R", 1000, 800 );
    canvasTagMassS3->Divide(2,2);
    TCanvas* canvasJetVsParton = new TCanvas( "jetVsParton", "jet Et VS parton Pt", 1000, 800 );
    canvasJetVsParton->Divide(2,2);
    int padCounter=0;
    for(;bDiscrHighEffMultiHistos_it!=bDiscrHighEffMultiHistos.end();
	  ++bDiscrHighEffMultiHistos_it,
	  ++nobDiscrHighEffMultiHistos_it,
	  ++bDiscrHighPurMultiHistos_it,
	  ++nobDiscrHighPurMultiHistos_it,
	  ++tagMassS1MultiHistos_it,
	  ++tagMassS2MultiHistos_it,
	  ++tagMassS3MultiHistos_it,
	  ++nobtagMassS1MultiHistos_it,
	  ++nobtagMassS2MultiHistos_it,
	  ++nobtagMassS3MultiHistos_it,
	  ++jetEtVSbPtMultiProfiles_it,
	  ++jetEtVSnobPtMultiProfiles_it,
	  ++padCounter){

      stringstream dRstringstream;
      dRstringstream<< dR_+dR_*(padCounter);
      TString title = "discrHighEffVsDeltaRStack";
      TString dRstring(dRstringstream.str()); 
      THStack* stack = new THStack(title + dRstring,title + dRstring);
      TLegend* legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
      Double_t integral_b = (*bDiscrHighEffMultiHistos_it)->GetEntries();
      Double_t integral_nob = (*nobDiscrHighEffMultiHistos_it)->GetEntries();
      if ( integral_b != 0 && integral_nob != 0 ) {
	(*bDiscrHighEffMultiHistos_it)->Scale(1./integral_b);
	(*nobDiscrHighEffMultiHistos_it)->Scale(1./integral_nob);
      }
      (*nobDiscrHighEffMultiHistos_it)->SetLineColor(kRed);
      stack->Add(*bDiscrHighEffMultiHistos_it);
      stack->Add(*nobDiscrHighEffMultiHistos_it);
      legend->AddEntry(*bDiscrHighEffMultiHistos_it, (*bDiscrHighEffMultiHistos_it)->GetName(), "l");
      legend->AddEntry(*nobDiscrHighEffMultiHistos_it, (*nobDiscrHighEffMultiHistos_it)->GetName(), "l");
      if(padCounter>1 && padCounter<6 ){
	canvasEff1->cd(padCounter-1);
	stack->Draw("nostack");
	legend->Draw();
      }
      if(padCounter>5 && padCounter<10){
	canvasEff2->cd(padCounter-5);
	stack->Draw("nostack");
	legend->Draw();
      }

      
      title = "discrHighPurVsDeltaRStack";
      stack = new THStack(title + dRstring,title + dRstring);
      legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
      integral_b = (*bDiscrHighPurMultiHistos_it)->Integral();
      integral_nob = (*nobDiscrHighPurMultiHistos_it)->Integral();
      if ( integral_b != 0 && integral_nob != 0 ) {
	(*bDiscrHighPurMultiHistos_it)->Scale(1./integral_b);
	(*nobDiscrHighPurMultiHistos_it)->Scale(1./integral_nob);
      }
      (*nobDiscrHighPurMultiHistos_it)->SetLineColor(kRed);
      stack->Add(*bDiscrHighPurMultiHistos_it);
      stack->Add(*nobDiscrHighPurMultiHistos_it);
      legend->AddEntry(*bDiscrHighPurMultiHistos_it, (*bDiscrHighPurMultiHistos_it)->GetName(), "l");
      legend->AddEntry(*nobDiscrHighPurMultiHistos_it, (*nobDiscrHighPurMultiHistos_it)->GetName(), "l");
      if(padCounter>1 && padCounter<6 ){
	canvasPur1->cd(padCounter-1);
	stack->Draw("nostack");
	legend->Draw();
      }
      if(padCounter>5 && padCounter<10){ 
	canvasPur2->cd(padCounter-5);
	stack->Draw("nostack");
	legend->Draw();
      }

      title = "tagMassVsDeltaRStack";
      stack = new THStack(title + dRstring,title + dRstring);
      legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
      Double_t integral_s1 = (*tagMassS1MultiHistos_it)->GetEntries();
      Double_t integral_s2 = (*tagMassS2MultiHistos_it)->GetEntries();
      Double_t integral_s3 = (*tagMassS3MultiHistos_it)->GetEntries();
      if ( integral_s1 != 0 && integral_s2 != 0 && integral_s3 != 0) {
	(*tagMassS1MultiHistos_it)->Scale(1./integral_s1);
	(*tagMassS2MultiHistos_it)->Scale(1./integral_s2);
	(*tagMassS3MultiHistos_it)->Scale(1./integral_s3);
      }
      (*tagMassS2MultiHistos_it)->SetLineColor(kBlue);
      (*tagMassS3MultiHistos_it)->SetLineColor(kGreen);
      stack->Add(*tagMassS1MultiHistos_it);
      stack->Add(*tagMassS2MultiHistos_it);
      stack->Add(*tagMassS3MultiHistos_it);
      legend->AddEntry(*tagMassS1MultiHistos_it, (*tagMassS1MultiHistos_it)->GetName(), "l");
      legend->AddEntry(*tagMassS2MultiHistos_it, (*tagMassS2MultiHistos_it)->GetName(), "l");
      legend->AddEntry(*tagMassS3MultiHistos_it, (*tagMassS3MultiHistos_it)->GetName(), "l");
      if(padCounter>1 && padCounter<6){
	canvasTagMass->cd(padCounter-1);
	stack->Draw("nostack");
	legend->Draw();
      }

      title = "tagMassS1VsDeltaRStack";
      stack = new THStack(title + dRstring,title + dRstring);
      legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );
      Double_t integral_nobs1 = (*nobtagMassS1MultiHistos_it)->GetEntries();
      if ( integral_s1 != 0 && integral_nos1 ) {
	(*nobtagMassS1MultiHistos_it)->Scale(1./integral_nos1);
      }
      (*nobtagMassS1MultiHistos_it)->SetLineColor(kRed);
      stack->Add(*tagMassS1MultiHistos_it);
      stack->Add(*tagMassS2MultiHistos_it);
      stack->Add(*tagMassS3MultiHistos_it);
      legend->AddEntry(*tagMassS1MultiHistos_it, (*tagMassS1MultiHistos_it)->GetName(), "l");
      legend->AddEntry(*tagMassS2MultiHistos_it, (*tagMassS2MultiHistos_it)->GetName(), "l");
      legend->AddEntry(*tagMassS3MultiHistos_it, (*tagMassS3MultiHistos_it)->GetName(), "l");
      if(padCounter>1 && padCounter<6){
	canvasTagMass->cd(padCounter-1);
	stack->Draw("nostack");
	legend->Draw();
      }

    }
    canvasEff1->Write();
    canvasEff2->Write();
    canvasPur1->Write();
    canvasPur2->Write();
    canvasTagMass->Write();

    Multi_bDiscrHighEff->Write();
    Multi_bDiscrHighPur->Write();
    Multi_nobDiscrHighEff->Write();
    Multi_nobDiscrHighPur->Write();
    Multi_tagMassS1->Write();
    Multi_tagMassS2->Write();
    Multi_tagMassS3->Write();
    Multi_jetEtVSbPt->Write();
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(TDAna);
