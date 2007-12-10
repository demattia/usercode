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
  MultidR_bDiscrHighEff   = new MultiTH1F( "multidRbDiscrHighEff",   "high efficiency discriminator for b jet",   200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  MultidR_bDiscrHighPur   = new MultiTH1F( "multidRbDiscrHighPur",   "high purity discriminator for b jet",       200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobDiscrHighEff = new MultiTH1F( "multidRnobDiscrHighEff", "high efficiency discriminator for no b jet",200,-20.,80., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobDiscrHighPur = new MultiTH1F( "multidRnobDiscrHighPur", "high purity discriminator for no b jet",    200,-20.,80., bins_, dR_, dRmax_, OutputFile );

  MultidR_tagMassS1       = new MultiTH1F( "multidRtagMassS1",       "tag mass S1 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS2       = new MultiTH1F( "multidRtagMassS2",       "tag mass S2 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS3       = new MultiTH1F( "multidRtagMassS3",       "tag mass S3 for b jet",                     150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS1    = new MultiTH1F( "multidRnobtagMassS1",    "tag mass S1for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS2    = new MultiTH1F( "multidRnobtagMassS2",    "tag mass S2for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS3    = new MultiTH1F( "multidRnobtagMassS3",    "tag mass S3for no b jet",                   150,  0.,15., bins_, dR_, dRmax_, OutputFile );

  MultidR_tagMassS1_loose   =new MultiTH1F( "multidRtagMassS1_loose",   "tag mass S1 for b jet loose b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS2_loose   =new MultiTH1F( "multidRtagMassS2_loose",   "tag mass S2 for b jet loose b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS3_loose   =new MultiTH1F( "multidRtagMassS3_loose",   "tag mass S3 for b jet loose b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS1_loose=new MultiTH1F( "multidRnobtagMassS1_loose","tag mass S1for no b jet loose b-tagging cut", 150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS2_loose=new MultiTH1F( "multidRnobtagMassS2_loose","tag mass S2for no b jet loose b-tagging cut", 150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS3_loose=new MultiTH1F( "multidRnobtagMassS3_loose","tag mass S3for no b jet loose b-tagging cut", 150,  0.,15., bins_, dR_, dRmax_, OutputFile );

  MultidR_tagMassS1_medium   =new MultiTH1F("multidRtagMassS1_medium",   "tag mass S1 for b jet medium b-tagging cut",  150,0.,15., bins_,dR_,dRmax_,OutputFile );
  MultidR_tagMassS2_medium   =new MultiTH1F("multidRtagMassS2_medium",   "tag mass S2 for b jet medium b-tagging cut",  150,0.,15., bins_,dR_,dRmax_,OutputFile );
  MultidR_tagMassS3_medium   =new MultiTH1F("multidRtagMassS3_medium",   "tag mass S3 for b jet medium b-tagging cut",  150,0.,15., bins_,dR_,dRmax_,OutputFile );
  MultidR_nobtagMassS1_medium=new MultiTH1F("multidRnobtagMassS1_medium","tag mass S1for no b jet medium b-tagging cut",150,0.,15., bins_,dR_,dRmax_,OutputFile );
  MultidR_nobtagMassS2_medium=new MultiTH1F("multidRnobtagMassS2_medium","tag mass S2for no b jet medium b-tagging cut",150,0.,15., bins_,dR_,dRmax_,OutputFile );
  MultidR_nobtagMassS3_medium=new MultiTH1F("multidRnobtagMassS3_medium","tag mass S3for no b jet medium b-tagging cut",150,0.,15., bins_,dR_,dRmax_,OutputFile );

  MultidR_tagMassS1_tight   =new MultiTH1F( "multidRtagMassS1_tight",   "tag mass S1 for b jet tight b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS2_tight   =new MultiTH1F( "multidRtagMassS2_tight",   "tag mass S2 for b jet tight b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_tagMassS3_tight   =new MultiTH1F( "multidRtagMassS3_tight",   "tag mass S3 for b jet tight b-tagging cut",  150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS1_tight=new MultiTH1F( "multidRnobtagMassS1_tight","tag mass S1for no b jet tight b-tagging cut",150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS2_tight=new MultiTH1F( "multidRnobtagMassS2_tight","tag mass S2for no b jet tight b-tagging cut",150,  0.,15., bins_, dR_, dRmax_, OutputFile );
  MultidR_nobtagMassS3_tight=new MultiTH1F( "multidRnobtagMassS3_tight","tag mass S3for no b jet tight b-tagging cut",150,  0.,15., bins_, dR_, dRmax_, OutputFile );

  MultidR_jetEtVSbPt      = new MultiTProfile("multidRjetEtVSbPt",   "jet Et as function as b parton Pt",   500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );
  MultidR_jetEtVSnobPt    = new MultiTProfile("multidRjetEtVSnobPt", "jet Et as function as no b parton Pt",500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );
  MultidR_jetUncorrEtVSbPt= new MultiTProfile("multidRjetUncorrEtVSbPt","jet uncorrected Et as function as b parton Pt",   500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );
  MultidR_jetUncorrEtVSnobPt= new MultiTProfile("multidRjetUncorrEtVSnobPt","jet uncorrected Et as function as no b parton Pt",500,0.,500.,0.,500., bins_, dR_, dRmax_, OutputFile );

  jets_   = 4;
  jetMin_ = 4;
  jetMax_ = jetMin_+jets_;
  MultijetNum_tagMassS1       = new MultiTH1F( "numJettagMassS1",       "tag mass S1 for b jet",                     150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS2       = new MultiTH1F( "numJettagMassS2",       "tag mass S2 for b jet",                     150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS3       = new MultiTH1F( "numJettagMassS3",       "tag mass S3 for b jet",                     150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS1    = new MultiTH1F( "numJetnobtagMassS1",    "tag mass S1for no b jet",                   150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS2    = new MultiTH1F( "numJetnobtagMassS2",    "tag mass S2for no b jet",                   150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS3    = new MultiTH1F( "numJetnobtagMassS3",    "tag mass S3for no b jet",                   150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );

  MultijetNum_tagMassS1_loose   =new MultiTH1F( "numJettagMassS1_loose",   "tag mass S1 for b jet loose b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS2_loose   =new MultiTH1F( "numJettagMassS2_loose",   "tag mass S2 for b jet loose b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS3_loose   =new MultiTH1F( "numJettagMassS3_loose",   "tag mass S3 for b jet loose b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS1_loose=new MultiTH1F( "numJetnobtagMassS1_loose","tag mass S1for no b jet loose b-tagging cut", 150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS2_loose=new MultiTH1F( "numJetnobtagMassS2_loose","tag mass S2for no b jet loose b-tagging cut", 150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS3_loose=new MultiTH1F( "numJetnobtagMassS3_loose","tag mass S3for no b jet loose b-tagging cut", 150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );

  MultijetNum_tagMassS1_medium   =new MultiTH1F("numJettagMassS1_medium",   "tag mass S1 for b jet medium b-tagging cut",  150,0.,15., jets_,jetMin_,jetMax_,OutputFile );
  MultijetNum_tagMassS2_medium   =new MultiTH1F("numJettagMassS2_medium",   "tag mass S2 for b jet medium b-tagging cut",  150,0.,15., jets_,jetMin_,jetMax_,OutputFile );
  MultijetNum_tagMassS3_medium   =new MultiTH1F("numJettagMassS3_medium",   "tag mass S3 for b jet medium b-tagging cut",  150,0.,15., jets_,jetMin_,jetMax_,OutputFile );
  MultijetNum_nobtagMassS1_medium=new MultiTH1F("numJetnobtagMassS1_medium","tag mass S1for no b jet medium b-tagging cut",150,0.,15., jets_,jetMin_,jetMax_,OutputFile );
  MultijetNum_nobtagMassS2_medium=new MultiTH1F("numJetnobtagMassS2_medium","tag mass S2for no b jet medium b-tagging cut",150,0.,15., jets_,jetMin_,jetMax_,OutputFile );
  MultijetNum_nobtagMassS3_medium=new MultiTH1F("numJetnobtagMassS3_medium","tag mass S3for no b jet medium b-tagging cut",150,0.,15., jets_,jetMin_,jetMax_,OutputFile );

  MultijetNum_tagMassS1_tight   =new MultiTH1F( "numJettagMassS1_tight",   "tag mass S1 for b jet tight b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS2_tight   =new MultiTH1F( "numJettagMassS2_tight",   "tag mass S2 for b jet tight b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_tagMassS3_tight   =new MultiTH1F( "numJettagMassS3_tight",   "tag mass S3 for b jet tight b-tagging cut",  150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS1_tight=new MultiTH1F( "numJetnobtagMassS1_tight","tag mass S1for no b jet tight b-tagging cut",150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS2_tight=new MultiTH1F( "numJetnobtagMassS2_tight","tag mass S2for no b jet tight b-tagging cut",150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );
  MultijetNum_nobtagMassS3_tight=new MultiTH1F( "numJetnobtagMassS3_tight","tag mass S3for no b jet tight b-tagging cut",150,  0.,15., jets_, jetMin_, jetMax_, OutputFile );

  deltaR = new TH1F("deltaR","deltaR between matched jet and b-parton",250,0.,1.0);


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

    // if deltaR cut [0.05;1.0]
    for( int numdR = 0; numdR < bins_; ++numdR){
      double dRcut = dR_ + dR_*(numdR);
      AssociatorEt<OfflineJet, MCParticle > associator( dRcut );      
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
	float jetUncorrEt = assocMap_it->first->uncorrEt();
	float partonPt    = assocMap_it->second->pt();

	if( fabs( assocMap_it->second->pid() ) == 5 ){
	  MultidR_bDiscrHighEff->Fill(discHighEff,numdR);
	  MultidR_bDiscrHighPur->Fill(discHighPur,numdR);  
	  if(tagMassS1 > 0.2) MultidR_tagMassS1->Fill(tagMassS1,numdR);
	  if(tagMassS2 > 0.2) MultidR_tagMassS2->Fill(tagMassS2,numdR);
	  if(tagMassS3 > 0.2) MultidR_tagMassS3->Fill(tagMassS3,numdR);

	  if(discHighEff > 2.3){
	    if(tagMassS1 > 0.2) MultidR_tagMassS1_loose->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_tagMassS2_loose->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_tagMassS3_loose->Fill(tagMassS3,numdR);
	  }
	  if(discHighEff > 5.3){
	    if(tagMassS1 > 0.2) MultidR_tagMassS1_medium->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_tagMassS2_medium->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_tagMassS3_medium->Fill(tagMassS3,numdR);
	  }
	  if(discHighPur > 4.8){
	    if(tagMassS1 > 0.2) MultidR_tagMassS1_tight->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_tagMassS2_tight->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_tagMassS3_tight->Fill(tagMassS3,numdR);
	  }
	  MultidR_jetEtVSbPt->Fill(partonPt,jetEt,numdR);
	  MultidR_jetUncorrEtVSbPt->Fill(partonPt,jetUncorrEt,numdR);
	  if(dRcut == 1.0){
	    double jetEta = assocMap_it->first->eta();
	    double jetPhi = assocMap_it->first->phi();
	    double bEta   = assocMap_it->second->eta();
	    double bPhi   = assocMap_it->second->phi();
	    double dR     = DeltaR(jetEta,jetPhi,bEta,bPhi);
	    deltaR->Fill(dR);
	  }
	}else{
	  MultidR_nobDiscrHighEff->Fill(discHighEff,numdR);
	  MultidR_nobDiscrHighPur->Fill(discHighPur,numdR);
	  if(tagMassS1 > 0.2) MultidR_nobtagMassS1->Fill(tagMassS1,numdR);
	  if(tagMassS2 > 0.2) MultidR_nobtagMassS2->Fill(tagMassS2,numdR);
	  if(tagMassS3 > 0.2) MultidR_nobtagMassS3->Fill(tagMassS3,numdR);

	  if(discHighEff > 2.3){
	    if(tagMassS1 > 0.2) MultidR_nobtagMassS1_loose->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_nobtagMassS2_loose->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_nobtagMassS3_loose->Fill(tagMassS3,numdR);
	  }
	  if(discHighEff > 5.3){
	    if(tagMassS1 > 0.2) MultidR_nobtagMassS1_medium->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_nobtagMassS2_medium->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_nobtagMassS3_medium->Fill(tagMassS3,numdR);
	  }
	  if(discHighPur > 4.8){
	    if(tagMassS1 > 0.2) MultidR_nobtagMassS1_tight->Fill(tagMassS1,numdR);
	    if(tagMassS2 > 0.2) MultidR_nobtagMassS2_tight->Fill(tagMassS2,numdR);
	    if(tagMassS3 > 0.2) MultidR_nobtagMassS3_tight->Fill(tagMassS3,numdR);
	  }
	  MultidR_jetEtVSnobPt->Fill(partonPt,jetEt,numdR);
	  MultidR_jetUncorrEtVSnobPt->Fill(partonPt,jetUncorrEt,numdR);
	}
      }
    }

    // if( goodIc5JetVec.size() >=4 ) {,>=5,>=6,>=7  goodIc5jets
    for( int numGoodIc5Jet = 0; numGoodIc5Jet < jets_; ++numGoodIc5Jet){
      int goodIc5Jetcut = jetMin_+numGoodIc5Jet;
      if( goodIc5JetVec.size() >= goodIc5Jetcut ) {
	AssociatorEt<OfflineJet, MCParticle > associator( 0.3 );
	std::auto_ptr<std::map<const OfflineJet*, 
  	                       const MCParticle*> > assocMap( associator.Associate( goodIc5JetVec, *MCpartons ) );
	
	std::map<const OfflineJet*, 
	         const MCParticle*>::const_iterator assocMap_it = assocMap->begin();
	for ( ; assocMap_it != assocMap->end(); ++assocMap_it ) {
	  // jet variables
	  float discHighPur      = assocMap_it->first->discriminatorHighPur();
	  float discHighEff      = assocMap_it->first->discriminatorHighEff();
	  float tagMassS1        = assocMap_it->first->tagTkMassS1();
	  float tagMassS2        = assocMap_it->first->tagTkMassS2();
	  float tagMassS3        = assocMap_it->first->tagTkMassS3();
	  float jetEt            = assocMap_it->first->et();
	  float jetUncorrEt      = assocMap_it->first->uncorrEt();
	  float emEnergyFraction = assocMap_it->first->emEnergyFraction();
	  float jetMass          = assocMap_it->first->jetMass();
	  int   tagNumTkS1       = assocMap_it->first->tkNumS1();
	  int   tagNumTkS2       = assocMap_it->first->tkNumS2();
	  int   tagNumTkS3       = assocMap_it->first->tkNumS3();
	  float tagSumPtTkS1     = assocMap_it->first->tkSumPtS1();
	  float tagSumPtTkS2     = assocMap_it->first->tkSumPtS2();
	  float tagSumPtTkS3     = assocMap_it->first->tkSumPtS3();

	  // MCparton variables
	  int   partonPid   = assocMap_it->second->mPid();
	  float partonPt    = assocMap_it->second->pt();

	  if( fabs( assocMap_it->second->pid() ) == 5 ){
	    if(tagMassS1 > 0.2) MultijetNum_tagMassS1->Fill(tagMassS1,numGoodIc5Jet);
	    if(tagMassS2 > 0.2) MultijetNum_tagMassS2->Fill(tagMassS2,numGoodIc5Jet);
	    if(tagMassS3 > 0.2) MultijetNum_tagMassS3->Fill(tagMassS3,numGoodIc5Jet);
	    
	    if(discHighEff > 2.3){
	      if(tagMassS1 > 0.2) MultijetNum_tagMassS1_loose->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_tagMassS2_loose->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_tagMassS3_loose->Fill(tagMassS3,numGoodIc5Jet);
	    }
	    if(discHighEff > 5.3){
	      if(tagMassS1 > 0.2) MultijetNum_tagMassS1_medium->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_tagMassS2_medium->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_tagMassS3_medium->Fill(tagMassS3,numGoodIc5Jet);
	    }
	    if(discHighPur > 4.8){
	      if(tagMassS1 > 0.2) MultijetNum_tagMassS1_tight->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_tagMassS2_tight->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_tagMassS3_tight->Fill(tagMassS3,numGoodIc5Jet);
	    }
	  }else{
	    if(tagMassS1 > 0.2) MultijetNum_nobtagMassS1->Fill(tagMassS1,numGoodIc5Jet);
	    if(tagMassS2 > 0.2) MultijetNum_nobtagMassS2->Fill(tagMassS2,numGoodIc5Jet);
	    if(tagMassS3 > 0.2) MultijetNum_nobtagMassS3->Fill(tagMassS3,numGoodIc5Jet);
	    
	    if(discHighEff > 2.3){
	      if(tagMassS1 > 0.2) MultijetNum_nobtagMassS1_loose->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_nobtagMassS2_loose->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_nobtagMassS3_loose->Fill(tagMassS3,numGoodIc5Jet);
	    }
	    if(discHighEff > 5.3){
	      if(tagMassS1 > 0.2) MultijetNum_nobtagMassS1_medium->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_nobtagMassS2_medium->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_nobtagMassS3_medium->Fill(tagMassS3,numGoodIc5Jet);
	    }
	    if(discHighPur > 4.8){
	      if(tagMassS1 > 0.2) MultijetNum_nobtagMassS1_tight->Fill(tagMassS1,numGoodIc5Jet);
	      if(tagMassS2 > 0.2) MultijetNum_nobtagMassS2_tight->Fill(tagMassS2,numGoodIc5Jet);
	      if(tagMassS3 > 0.2) MultijetNum_nobtagMassS3_tight->Fill(tagMassS3,numGoodIc5Jet);
	    }
	    
	  }
	}
      }
    }

  // Take the first (and it should be the only one) element in the map and 
  // get the eta of the closest jet (in DeltaR)
  //   double JetPhi = (*AssocMap).begin()->second->phi();

  // Do not use DR association. The MEt has no z component. Use DeltaPhi association.


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
    
  MultiStack<MultiTH1F,TH1F>* canvasdREff        = new MultiStack<MultiTH1F,TH1F>( "discrHighEffVsDeltaRCanvas", 
										 "b tag discriminator in high efficiency VS delta R", 
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRPur        = new MultiStack<MultiTH1F,TH1F>( "discrHighPurVsDeltaRCanvas", 
										 "b tag discriminator in high purity VS delta R",
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS1  = new MultiStack<MultiTH1F,TH1F>( "tagMassS1VsDeltaRCanvas", 
										 "tag mass S1 VS delta R",
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS2  = new MultiStack<MultiTH1F,TH1F>( "tagMassS2VsDeltaRCanvas",
										 "tag mass S2 VS delta R",
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS3  = new MultiStack<MultiTH1F,TH1F>( "tagMassS3VsDeltaRCanvas", 
										 "tag mass S3 VS delta R",
										 1000, 800, 4, 3, 11, true );

  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS1_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_looseVsDeltaRCanvas",
										      "tag mass S1 for loose b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS2_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_looseVsDeltaRCanvas",
										      "tag mass S2 for loose b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS3_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_looseVsDeltaRCanvas",
										      "tag mass S3 for loose b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS1_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_mediumVsDeltaRCanvas",
										       "tag mass S1 for medium b-tagging cut VS delta R", 
										       1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS2_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_mediumVsDeltaRCanvas",
										       "tag mass S2 for medium b-tagging cut VS delta R", 
										       1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS3_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_mediumVsDeltaRCanvas",
										       "tag mass S3 for medium b-tagging cut VS delta R", 
										       1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS1_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_tightVsDeltaRCanvas",
										      "tag mass S1 for tight b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS2_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_tightVsDeltaRCanvas",
										      "tag mass S2 for tight b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasdRTagMassS3_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_tightVsDeltaRCanvas",
										      "tag mass S3 for tight b-tagging cut VS delta R", 
										      1000, 800, 4, 3, 11, true );
  MultiStack<MultiTProfile,TProfile>* canvasdRJetVsParton= new MultiStack<MultiTProfile,TProfile>( "jetVsParton",
												 "jet Et VS parton Pt", 
												 1000, 800, 4, 3, 11, false );
  MultiStack<MultiTProfile,TProfile>* canvasdRUncorrJetVsParton= new MultiStack<MultiTProfile,TProfile>( "uncorrJetVsParton",
												       "jet uncorrected Et VS parton Pt", 
												       1000, 800, 4, 3, 11, false );

  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS1  = new MultiStack<MultiTH1F,TH1F>( "tagMassS1VsDeltaRCanvas", 
										 "tag mass S1 VS delta R",
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS2  = new MultiStack<MultiTH1F,TH1F>( "tagMassS2VsDeltaRCanvas",
										 "tag mass S2 VS delta R",
										 1000, 800, 4, 3, 11, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS3  = new MultiStack<MultiTH1F,TH1F>( "tagMassS3VsDeltaRCanvas", 
										 "tag mass S3 VS delta R",
										 1000, 800, 4, 3, 11, true );

  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS1_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_looseVsJetNumCanvas",
										      "tag mass S1 for loose b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS2_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_looseVsJetNumCanvas",
										      "tag mass S2 for loose b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS3_loose = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_looseVsJetNumCanvas",
										      "tag mass S3 for loose b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS1_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_mediumVsJetNumCanvas",
										       "tag mass S1 for medium b-tagging cut VS good jet number cut ", 
										       1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS2_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_mediumVsJetNumCanvas",
										       "tag mass S2 for medium b-tagging cut VS good jet number cut ", 
										       1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS3_medium = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_mediumVsJetNumCanvas",
										       "tag mass S3 for medium b-tagging cut VS good jet number cut ", 
										       1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS1_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS1_tightVsJetNumCanvas",
										      "tag mass S1 for tight b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS2_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS2_tightVsJetNumCanvas",
										      "tag mass S2 for tight b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );
  MultiStack<MultiTH1F,TH1F>* canvasjetNumTagMassS3_tight = new MultiStack<MultiTH1F,TH1F>( "tagMassS3_tightVsJetNumCanvas",
										      "tag mass S3 for tight b-tagging cut VS good jet number cut ", 
										      1000, 800, 4, 1, 5, true );

  canvasdREff->Fill(MultidR_bDiscrHighEff,MultidR_nobDiscrHighEff);
  canvasdRPur->Fill(MultidR_bDiscrHighPur,MultidR_nobDiscrHighPur);
  canvasdRTagMassS1->Fill(MultidR_tagMassS1,MultidR_nobtagMassS1);
  canvasdRTagMassS1_loose->Fill(MultidR_tagMassS1_loose,MultidR_nobtagMassS1_loose);
  canvasdRTagMassS1_medium->Fill(MultidR_tagMassS1_medium,MultidR_nobtagMassS1_medium);
  canvasdRTagMassS1_tight->Fill(MultidR_tagMassS1_tight,MultidR_nobtagMassS1_tight);
  canvasdRTagMassS2->Fill(MultidR_tagMassS2,MultidR_nobtagMassS2);
  canvasdRTagMassS2_loose->Fill(MultidR_tagMassS2_loose,MultidR_nobtagMassS2_loose);
  canvasdRTagMassS2_medium->Fill(MultidR_tagMassS2_medium,MultidR_nobtagMassS2_medium);
  canvasdRTagMassS2_tight->Fill(MultidR_tagMassS2_tight,MultidR_nobtagMassS2_tight);
  canvasdRTagMassS3->Fill(MultidR_tagMassS3,MultidR_nobtagMassS3);
  canvasdRTagMassS3_loose->Fill(MultidR_tagMassS3_loose,MultidR_nobtagMassS3_loose);
  canvasdRTagMassS3_medium->Fill(MultidR_tagMassS3_medium,MultidR_nobtagMassS3_medium);
  canvasdRTagMassS3_tight->Fill(MultidR_tagMassS3_tight,MultidR_nobtagMassS3_tight);
  canvasdRJetVsParton->Fill(MultidR_jetEtVSbPt,MultidR_jetEtVSnobPt);
  canvasdRUncorrJetVsParton->Fill(MultidR_jetUncorrEtVSbPt,MultidR_jetUncorrEtVSnobPt);

  canvasjetNumTagMassS1->Fill(MultijetNum_tagMassS1,MultijetNum_nobtagMassS1);
  canvasjetNumTagMassS1_loose->Fill(MultijetNum_tagMassS1_loose,MultijetNum_nobtagMassS1_loose);
  canvasjetNumTagMassS1_medium->Fill(MultijetNum_tagMassS1_medium,MultijetNum_nobtagMassS1_medium);
  canvasjetNumTagMassS1_tight->Fill(MultijetNum_tagMassS1_tight,MultijetNum_nobtagMassS1_tight);
  canvasjetNumTagMassS2->Fill(MultijetNum_tagMassS2,MultijetNum_nobtagMassS2);
  canvasjetNumTagMassS2_loose->Fill(MultijetNum_tagMassS2_loose,MultijetNum_nobtagMassS2_loose);
  canvasjetNumTagMassS2_medium->Fill(MultijetNum_tagMassS2_medium,MultijetNum_nobtagMassS2_medium);
  canvasjetNumTagMassS2_tight->Fill(MultijetNum_tagMassS2_tight,MultijetNum_nobtagMassS2_tight);
  canvasjetNumTagMassS3->Fill(MultijetNum_tagMassS3,MultijetNum_nobtagMassS3);
  canvasjetNumTagMassS3_loose->Fill(MultijetNum_tagMassS3_loose,MultijetNum_nobtagMassS3_loose);
  canvasjetNumTagMassS3_medium->Fill(MultijetNum_tagMassS3_medium,MultijetNum_nobtagMassS3_medium);
  canvasjetNumTagMassS3_tight->Fill(MultijetNum_tagMassS3_tight,MultijetNum_nobtagMassS3_tight);

  canvasdREff->Write();
  canvasdRPur->Write();
  canvasdRTagMassS1->Write();
  canvasdRTagMassS1_loose->Write();
  canvasdRTagMassS1_medium->Write();
  canvasdRTagMassS1_tight->Write();
  canvasdRTagMassS2->Write();
  canvasdRTagMassS2_loose->Write();
  canvasdRTagMassS2_medium->Write();
  canvasdRTagMassS2_tight->Write();
  canvasdRTagMassS3->Write();
  canvasdRTagMassS3_loose->Write();
  canvasdRTagMassS3_medium->Write();
  canvasdRTagMassS3_tight->Write();
  canvasdRJetVsParton->Write();
  canvasdRUncorrJetVsParton->Write();

  /*
  canvasjetNumTagMassS1->Write();
  canvasjetNumTagMassS1_loose->Write();
  canvasjetNumTagMassS1_medium->Write();
  canvasjetNumTagMassS1_tight->Write();
  canvasjetNumTagMassS2->Write();
  canvasjetNumTagMassS2_loose->Write();
  canvasjetNumTagMassS2_medium->Write();
  canvasjetNumTagMassS2_tight->Write();
  canvasjetNumTagMassS3->Write();
  canvasjetNumTagMassS3_loose->Write();
  canvasjetNumTagMassS3_medium->Write();
  canvasjetNumTagMassS3_tight->Write();
  */

  MultidR_bDiscrHighEff->Write();
  MultidR_bDiscrHighPur->Write();
  MultidR_nobDiscrHighEff->Write();
  MultidR_nobDiscrHighPur->Write();
  MultidR_tagMassS1->Write();
  MultidR_tagMassS2->Write();
  MultidR_tagMassS3->Write();
  MultidR_tagMassS1_loose->Write();
  MultidR_tagMassS2_loose->Write();
  MultidR_tagMassS3_loose->Write();
  MultidR_tagMassS1_medium->Write();
  MultidR_tagMassS2_medium->Write();
  MultidR_tagMassS3_medium->Write();
  MultidR_tagMassS1_tight->Write();
  MultidR_tagMassS2_tight->Write();
  MultidR_tagMassS3_tight->Write();
  MultidR_nobtagMassS1->Write();
  MultidR_nobtagMassS2->Write();
  MultidR_nobtagMassS3->Write();
  MultidR_nobtagMassS1_loose->Write();
  MultidR_nobtagMassS2_loose->Write();
  MultidR_nobtagMassS3_loose->Write();
  MultidR_nobtagMassS1_medium->Write();
  MultidR_nobtagMassS2_medium->Write();
  MultidR_nobtagMassS3_medium->Write();
  MultidR_nobtagMassS1_tight->Write();
  MultidR_nobtagMassS2_tight->Write();
  MultidR_nobtagMassS3_tight->Write();
  MultidR_jetEtVSbPt->Write();
  MultidR_jetEtVSnobPt->Write();
  MultidR_jetUncorrEtVSbPt->Write();
  MultidR_jetUncorrEtVSnobPt->Write();

  /*
  MultijetNum_tagMassS1->Write();
  MultijetNum_tagMassS2->Write();
  MultijetNum_tagMassS3->Write();
  MultijetNum_tagMassS1_loose->Write();
  MultijetNum_tagMassS2_loose->Write();
  MultijetNum_tagMassS3_loose->Write();
  MultijetNum_tagMassS1_medium->Write();
  MultijetNum_tagMassS2_medium->Write();
  MultijetNum_tagMassS3_medium->Write();
  MultijetNum_tagMassS1_tight->Write();
  MultijetNum_tagMassS2_tight->Write();
  MultijetNum_tagMassS3_tight->Write();
  MultijetNum_nobtagMassS1->Write();
  MultijetNum_nobtagMassS2->Write();
  MultijetNum_nobtagMassS3->Write();
  MultijetNum_nobtagMassS1_loose->Write();
  MultijetNum_nobtagMassS2_loose->Write();
  MultijetNum_nobtagMassS3_loose->Write();
  MultijetNum_nobtagMassS1_medium->Write();
  MultijetNum_nobtagMassS2_medium->Write();
  MultijetNum_nobtagMassS3_medium->Write();
  MultijetNum_nobtagMassS1_tight->Write();
  MultijetNum_nobtagMassS2_tight->Write();
  MultijetNum_nobtagMassS3_tight->Write();
  */
}

//define this as a plug-in
DEFINE_FWK_MODULE(TDAna);
