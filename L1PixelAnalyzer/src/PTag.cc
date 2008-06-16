// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/PTag.h"

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

L1Trig PTag::L1Trigger;

// Constructors and destructor
// ---------------------------
PTag::PTag(const edm::ParameterSet& iConfig) :
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
  // ----------------------------------------------------------
  Multi_Vertex_Dz_ = new MultiTH1F( "Vertex_Dz", "Minimum distance between verteces", 
				    50, 0., 10., bins_, dz_, dzmax_, OutputFile );
  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 50, 0, 3.15 );

  
  NJets_ = new TH1F ( "NJets", "Number of selected jets", 20, 0, 20 );
  UncorrHt_ = new TH1F ( "UncorrHt", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHt_ = new TH1F ( "CorrHt", "Ht with corrected jets", 50, 0, 4000 );
  GoodHt_ = new TH1F ( "GoodHt", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2_ = new TH1F ( "GoodHt2", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEt_ = new TH1F ( "UncorrSumEt", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEt_ = new TH1F ( "CorrSumEt", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEt_ = new TH1F ( "GoodSumEt", "SumEt with selected jets", 50, 0, 4000 );
  MEt_ = new TH1F ( "MEt", "MEt", 50, 0, 500 );
  MEtSig_ = new TH1F ( "MEtSig", "MEtSig", 50, 0, 20 );
  MEtSigNew_ = new TH1F ( "MEtSigNew", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPM_ = new TH1F ( "MEtDPM", "MEtDPM", 50, 0, 3.2 );
  MEtDP1_ = new TH1F ( "MEtDP1", "MEtDP1", 50, 0, 3.2 );
  MEtDP2_ = new TH1F ( "MEtDP2", "MEtDP2", 50, 0, 3.2 );
  MEtDP3_ = new TH1F ( "MEtDP3", "MEtDP3", 50, 0, 3.2 );
  UncorrMEtSig_ = new TH1F ( "UncorrMEtSig", "MEtSig", 50, 0, 20 );
  CorrMEtSig_ = new TH1F ( "CorrMEtSig", "MEtSig", 50, 0, 20 );
  M3best_ = new TH1F ( "M3best", "M3best", 50, 0, 300 );
  Mwbest_ = new TH1F ( "Mwbest", "Mwbest", 50, 0, 200 );
  Chi2mass_ = new TH1F ( "Chi2mass", "Chi2mass", 50, 0, 50 );
  M45best_ = new TH1F ( "M45best", "M45best", 50, 0, 300 );
  Chi2ext_ = new TH1F ( "Chi2ext", "Chi2ext", 50, 0, 50 );
  MEx_SumEt_ = new TH2F ( "MEx_SumEt", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12_ = new TH1F ( "DP12", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbb_ = new TH1F ( "DPbb", "Delta phi bb", 50, 0., 3.2 );
  M_others_ = new TH1F ( "M_others", "Mass of other jets", 50, 0., 1500. );
  Mbbnoh_ = new TH1F ( "Mbbnoh", "Mass of b pair not from h", 50, 0., 400 );
  DPbbnoh_ = new TH1F ( "DPbbnoh", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8_ = new TH1F ( "M8", "Mass 8 jets", 50, 0., 2000.) ;
  C8_ = new TH1F ( "C8", "Centrality 8 jets", 50, 0., 1.2 );

  NJetsS_ = new TH1F ( "NJetsS", "Number of selected jets", 20, 0, 20 );
  UncorrHtS_ = new TH1F ( "UncorrHtS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtS_ = new TH1F ( "CorrHtS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtS_ = new TH1F ( "GoodHtS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2S_ = new TH1F ( "GoodHt2S", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtS_ = new TH1F ( "UncorrSumEtS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtS_ = new TH1F ( "CorrSumEtS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtS_ = new TH1F ( "GoodSumEtS", "SumEt with selected jets", 50, 0, 4000 );
  MEtS_ = new TH1F ( "MEtS", "MEt", 50, 0, 500 );
  MEtSigS_ = new TH1F ( "MEtSigS", "MEtSig", 50, 0, 20 );
  MEtSigNewS_ = new TH1F ( "MEtSigNewS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMS_ = new TH1F ( "MEtDPMS", "MEtDPM", 50, 0, 3.2 );
  MEtDP1S_ = new TH1F ( "MEtDP1S", "MEtDP1S", 50, 0, 3.2 );
  MEtDP2S_ = new TH1F ( "MEtDP2S", "MEtDP2S", 50, 0, 3.2 );
  MEtDP3S_ = new TH1F ( "MEtDP3S", "MEtDP3S", 50, 0, 3.2 );
  UncorrMEtSigS_ = new TH1F ( "UncorrMEtSigS", "MEtSig", 50, 0, 20 );
  CorrMEtSigS_ = new TH1F ( "CorrMEtSigS", "MEtSig", 50, 0, 20 );
  M3bestS_ = new TH1F ( "M3bestS", "M3bestS", 50, 0, 300 );
  MwbestS_ = new TH1F ( "MwbestS", "MwbestS", 50, 0, 200 );
  Chi2massS_ = new TH1F ( "Chi2massS", "Chi2massS", 50, 0, 50 );
  M45bestS_ = new TH1F ( "M45bestS", "M45best", 50, 0, 300 );
  Chi2extS_ = new TH1F ( "Chi2extS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtS_ = new TH2F ( "MEx_SumEtS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12S_ = new TH1F ( "DP12S", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbS_ = new TH1F ( "DPbbS", "Delta phi bb", 50, 0., 3.2 );
  M_othersS_ = new TH1F ( "M_othersS", "Mass of other jets", 50, 0., 1500. );
  MbbnohS_ = new TH1F ( "MbbnohS", "Mass of b pair not from h", 50, 0., 400 );
  DPbbnohS_ = new TH1F ( "DPbbnohS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8S_ = new TH1F ( "M8S", "Mass 8 jets", 50, 0., 2000.) ;
  C8S_ = new TH1F ( "C8S", "Centrality 8 jets", 50, 0., 1.2 );

  NJetsSS_ = new TH1F ( "NJetsSS", "Number of selected jets", 20, 0, 20 );
  UncorrHtSS_ = new TH1F ( "UncorrHtSS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSS_ = new TH1F ( "CorrHtSS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSS_ = new TH1F ( "GoodHtSS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SS_ = new TH1F ( "GoodHt2SS", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSS_ = new TH1F ( "UncorrSumEtSS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSS_ = new TH1F ( "CorrSumEtSS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSS_ = new TH1F ( "GoodSumEtSS", "SumEt with selected jets", 50, 0, 4000 );
  MEtSS_ = new TH1F ( "MEtSS", "MEt", 50, 0, 500 );
  MEtSigSS_ = new TH1F ( "MEtSigSS", "MEtSig", 50, 0, 20 );
  MEtSigNewSS_ = new TH1F ( "MEtSigNewSS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSS_ = new TH1F ( "MEtDPMSS", "MEtDPMSS", 50, 0, 3.2 );
  MEtDP1SS_ = new TH1F ( "MEtDP1SS", "MEtDP1SS", 50, 0, 3.2 );
  MEtDP2SS_ = new TH1F ( "MEtDP2SS", "MEtDP2SS", 50, 0, 3.2 );
  MEtDP3SS_ = new TH1F ( "MEtDP3SS", "MEtDP3SS", 50, 0, 3.2 );
  UncorrMEtSigSS_ = new TH1F ( "UncorrMEtSigSS", "UncorrMEtSigSS", 50, 0, 20 );
  CorrMEtSigSS_ = new TH1F ( "CorrMEtSigSS", "CorrMEtSigSS", 50, 0, 20 );
  M3bestSS_ = new TH1F ( "M3bestSS", "M3bestSS", 50, 0, 300 );
  MwbestSS_ = new TH1F ( "MwbestSS", "MwbestSS", 50, 0, 200 );
  Chi2massSS_ = new TH1F ( "Chi2massSS", "Chi2massSS", 50, 0, 50 );
  M45bestSS_ = new TH1F ( "M45bestSS", "M45best", 50, 0, 300 );
  Chi2extSS_ = new TH1F ( "Chi2extSS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSS_ = new TH2F ( "MEx_SumEtSS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SS_ = new TH1F ( "DP12SS", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSS_ = new TH1F ( "DPbbSS", "Delta phi bb", 50, 0., 3.2 );
  M_othersSS_ = new TH1F ( "M_othersSS", "Mass of other jets", 50, 0., 1500. );
  MbbnohSS_ = new TH1F ( "MbbnohSS", "Mass of b pair not from h", 50, 0., 400 );
  DPbbnohSS_ = new TH1F ( "DPbbnohSS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SS_ = new TH1F ( "M8SS", "Mass 8 jets", 50, 0., 2000.) ;
  C8SS_ = new TH1F ( "C8SS", "Centrality 8 jets", 50, 0., 1.2 );

  NJetsSSS_ = new TH1F ( "NJetsSSS", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSS_ = new TH1F ( "UncorrHtSSS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSS_ = new TH1F ( "CorrHtSSS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSS_ = new TH1F ( "GoodHtSSS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSS_ = new TH1F ( "GoodHt2SSS", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSS_ = new TH1F ( "UncorrSumEtSSS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSS_ = new TH1F ( "CorrSumEtSSS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSS_ = new TH1F ( "GoodSumEtSSS", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSS_ = new TH1F ( "MEtSSS", "MEt", 50, 0, 500 );
  MEtSigSSS_ = new TH1F ( "MEtSigSSS", "MEtSig", 50, 0, 20 );
  MEtSigNewSSS_ = new TH1F ( "MEtSigNewSSS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSS_ = new TH1F ( "MEtDPMSSS", "MEtDPMSSS", 50, 0, 3.2 );
  MEtDP1SSS_ = new TH1F ( "MEtDP1SSS", "MEtDP1SSS", 50, 0, 3.2 );
  MEtDP2SSS_ = new TH1F ( "MEtDP2SSS", "MEtDP2SSS", 50, 0, 3.2 );
  MEtDP3SSS_ = new TH1F ( "MEtDP3SSS", "MEtDP3SSS", 50, 0, 3.2 );
  UncorrMEtSigSSS_ = new TH1F ( "UncorrMEtSigSSS", "UncorrMEtSigSSS", 50, 0, 20 );
  CorrMEtSigSSS_ = new TH1F ( "CorrMEtSigSSS", "CorrMEtSigSSS", 50, 0, 20 );
  M3bestSSS_ = new TH1F ( "M3bestSSS", "M3bestSSS", 50, 0, 300 );
  MwbestSSS_ = new TH1F ( "MwbestSSS", "MwbestSSS", 50, 0, 200 );
  Chi2massSSS_ = new TH1F ( "Chi2massSSS", "Chi2massSSS", 50, 0, 50 );
  M45bestSSS_ = new TH1F ( "M45bestSSS", "M45best", 50, 0, 300 );
  Chi2extSSS_ = new TH1F ( "Chi2extSSS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSS_ = new TH2F ( "MEx_SumEtSSS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSS_ = new TH1F ( "DP12SSS", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSS_ = new TH1F ( "DPbbSSS", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSS_ = new TH1F ( "M_othersSSS", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSS_ = new TH1F ( "MbbnohSSS", "Mass of b pair not from h", 50, 0., 400 );
  DPbbnohSSS_ = new TH1F ( "DPbbnohSSS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SSS_ = new TH1F ( "M8SSS", "Mass 8 jets", 50, 0., 2000.) ;
  C8SSS_ = new TH1F ( "C8SSS", "Centrality 8 jets", 50, 0., 1.2 );

  UncorrMEt_SumEt_ = new TH2F ( "UncorrMEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEt_ = new TH2F ( "CorrMEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEt_ = new TH2F ( "MEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  UncorrMEt_SumEtC_ = new TH2F ( "UncorrMEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEtC_ = new TH2F ( "CorrMEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEtC_ = new TH2F ( "MEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  UncorrMEt_SumEtJ_ = new TH2F ( "UncorrMEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEtJ_ = new TH2F ( "CorrMEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEtJ_ = new TH2F ( "MEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  

  dR_    = 0.05;
  dRmax_ = dR_ + dR_*(bins_);

  loose_  = 2.3;
  medium_ = 5.3;
  tight_  = 4.8;

}


PTag::~PTag()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  // Draw the histograms
  //  HiVar->Plot();
  // -------------------

  OutputFile->Write();

}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void PTag::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  // Count IC5 jets with Et>=30 GeV and |eta| < 3.0
  // ----------------------------------------------
  double goodIc5Jets = 0; // Number of selected jets
  double NTagsL = 0;      // Number of loose b-tags
  double NTags = 0;       // Number of b-tags in the selected jets
  double NTagsT = 0;      // Number of tight b-tags

  double UncorrSumEx=0;
  double UncorrSumEy=0;
  double UncorrSumEt=0;  // SumEt with all jets before corrections
  double UncorrHt=0;     // Ht with all jets before corrections

  double CorrSumEx=0;
  double CorrSumEy=0;
  double CorrSumEt=0;    // SumEt with all jets after corrections 
  double CorrHt=0;       // Ht with all jets after corrections

  double GoodSumEt=0;    // SumEt with all selected jets
  double GoodHt=0;       // Ht with all selected jets
  double GoodHt2=0;      // Ht with all selected jets and standard MET

  double dpmin=20;  // DP min MET / selected jets
  double dp1st=20;  // DP MET/leading jet
  double et1st=0;   // Et of leading jet among those selected by cuts in et, eta
  double dp2nd=20;  // DP MET/2nd leading jet
  double et2nd=0;   // Et of second leading jet
  double dp3rd=20;  // DP MET/3rd leading jet
  double et3rd=0;   // Et of third leading jet

  double m3best=-999; // Massa del tripletto piu' vicino a 172 GeV
  double mwbest=-999;    // Massa del doppietto piu' vicino a 80.4 GeV tra i tre nel tripletto definito sopra
  double chi2=1000;   // Chiquadro costruito con le due masse qui sopra
  double m45best=-999;
  double chi2ext=1000;
  double dp12 = 0;        // Delta phi jet 1 jet 2
  double phi1st = 0;      // phi 1st jet (used to compute dp12)
  double phi2nd = 0;      // phi 2nd jet (       "            )
  double dpbb = 0;    // delta phi bb
  double m_others = -999; // mass of jets not part of best triplet and pair
  double mbbnohmax = -999; // mass of b pair with highest mass excluding b-jets from h decay
  double dpbbnohmax = 0; // angle of b pair with highest mass excluding b-jets from h decay

  // Offline cuts
  // ------------
  bool NJetsCut = false;   // Whether event passes NJet cut
  bool MEtSigCut = false;  // Whether event pass MET significance cut
  bool NTagsCut = false;

  double m8=0;
  double set=0;
  double c8=0;

  // good jets array
  // ---------------
  int iJ=0;
  double JEt[100];
  double Jeta[100];
  double Jphi[100];
  double JHET[100];
  double JHPT[100];
  double NHSJ = 8; // number of jets from hard subprocess  

  vector<OfflineJet> goodIc5JetVec;
  if ( caloJets->size() != 0 ) {
    std::vector<double> vec_DPhi;
    vec_DPhi.reserve(caloJets->size());
    vector<SimpleJet> vec_calojet; 
    vec_calojet.reserve(caloJets->size());

    // Correct offline jets on the fly
    // -------------------------------
    for( OfflineJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal ) {
      vec_calojet.push_back( SimpleJet( cal->et(), cal->eta(), cal->phi() ) );
      uncorr_JetPt_IC5_->Fill( cal->uncorrEt() );
      corr_JetPt_IC5_->Fill( cal->et() );

      UncorrSumEx+=cal->uncorrEt()*cos(cal->phi());
      UncorrSumEy+=cal->uncorrEt()*sin(cal->phi());
      UncorrSumEt+=cal->uncorrEt();
      UncorrHt+=cal->uncorrEt();
      CorrSumEx+=cal->et()*cos(cal->phi());
      CorrSumEy+=cal->et()*sin(cal->phi());
      CorrSumEt+=cal->et();
      CorrHt+=cal->et();

      if ( cal->et() >= 30. && fabs( cal->eta() ) < 2.5 ) {
        goodIc5Jets=goodIc5Jets+1;
	if ( cal->discriminatorHighEff() > loose_ ) NTagsL++;
	if ( cal->discriminatorHighEff() > medium_ ) NTags++;
	if ( cal->discriminatorHighEff() > tight_ ) NTagsT++;
	GoodHt+=cal->et();
	GoodHt2+=cal->et();
	GoodSumEt+=cal->et();
	double dp = 3.1415926-fabs(fabs(cal->phi()-caloMET->phi())-3.1415926);
	if ( dp<dpmin ) dpmin=dp;
	if ( cal->et()>et1st ) {
	  dp1st=dp;
	  et1st=cal->et();
	  phi1st=cal->phi();
	} else if ( cal->et()>et2nd ) {
	  dp2nd=dp;
	  et2nd=cal->et();
	  phi2nd=cal->phi();
	} else if ( cal->et()>et3rd ) {
	  dp3rd=dp;
	  et3rd=cal->et();
	}
	goodIc5JetVec.push_back(*cal);

	// Store information on good jets in arrays
	// ----------------------------------------
	if ( iJ<100 ) {
	  JEt[iJ]=cal->et();
	  Jeta[iJ]=cal->eta();
	  Jphi[iJ]=cal->phi();
	  JHET[iJ]=cal->discriminatorHighEff();
	  JHPT[iJ]=cal->discriminatorHighPur();
	  iJ++;
	}
      }
    } // End loop on cal jet collection

    // Sort good jets by Et
    // --------------------
    for ( int times=0; times<iJ; times++ ) {
      for ( int j=0; j<iJ-1; j++ ) {
	int k = iJ-j-1;
	if ( JEt[k]>JEt[k-1] ) {
	  double a1 = JEt[k];
	  double a2 = Jeta[k];
	  double a3 = Jphi[k];
	  double a4 = JHET[k];
	  double a5 = JHPT[k];
	  JEt[k]=JEt[k-1];
	  Jeta[k]=Jeta[k-1];
	  Jphi[k]=Jphi[k-1];
	  JHET[k]=JHET[k-1];
	  JHPT[k]=JHPT[k-1];
	  JEt[k-1]=a1;
	  Jeta[k-1]=a2;
	  Jphi[k-1]=a3;
	  JHET[k-1]=a4;
	  JHPT[k-1]=a5;

	}
      }
    }

    dp12 = 3.1415926-fabs(fabs(phi1st-phi2nd)-3.1415926);
    JetNumber_IC5_->Fill( goodIc5Jets );
    HiVar->Fill( vec_calojet );
    
    // Compute best top and W masses and determine simple chi2
    // -------------------------------------------------------
    double px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,e1,e2,e3;
    double px4,px5,py4,py5,pz4,pz5,e4,e5;
    double px,py,pz,e;
    double spx,spy,spz,se;
    double m3,m12,m13,m23;
    double m45;
    double mw=0;
    double bestWhasb=0;
    double bestbnob=0;
    int it1=-1;
    int it2=-1;
    int it3=-1;
    int ih1=-1;
    int ih2=-1;
    for ( int i=0; i<iJ-2 && i<NHSJ-2; i++ ) {
      for ( int j=i+1; j<iJ-1 && j<NHSJ-1; j++ ) {
	for ( int k=j+1; k<iJ && k<NHSJ; k++ ) { // limit to 8 the number of jets in hard subprocess

	  // Require that the triplet owns exactly one b-tag
	  // -----------------------------------------------
	  int nb=0;
	  if ( JHET[i]>loose_ ) nb++; 
	  if ( JHET[j]>loose_ ) nb++;
	  if ( JHET[k]>loose_ ) nb++;
	  if ( nb>=1 ) {
	    bestWhasb=0;
	    px1=JEt[i]*cos(Jphi[i]);
	    px2=JEt[j]*cos(Jphi[j]);
	    px3=JEt[k]*cos(Jphi[k]);
	    py1=JEt[i]*sin(Jphi[i]);
	    py2=JEt[j]*sin(Jphi[j]);
	    py3=JEt[k]*sin(Jphi[k]);
	    pz1=JEt[i]/tan(2*atan(exp(-Jeta[i])));
	    pz2=JEt[j]/tan(2*atan(exp(-Jeta[j])));
	    pz3=JEt[k]/tan(2*atan(exp(-Jeta[k])));
	    e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
	    e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
	    e3 =sqrt(px3*px3+py3*py3+pz3*pz3);
	    m3=(e1+e2+e3)*(e1+e2+e3)-(px1+px2+px3)*(px1+px2+px3)-
	      (py1+py2+py3)*(py1+py2+py3)-(pz1+pz2+pz3)*(pz1+pz2+pz3);
	    if ( m3>0 ) m3=sqrt(m3);
	    m12=(e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2);
	    if ( m12>0 ) m12=sqrt(m12);
	    m13=(e1+e3)*(e1+e3)-(px1+px3)*(px1+px3)-(py1+py3)*(py1+py3)-(pz1+pz3)*(pz1+pz3);
	    if ( m13>0 ) m13=sqrt(m13);
	    m23=(e2+e3)*(e2+e3)-(px2+px3)*(px2+px3)-(py2+py3)*(py2+py3)-(pz2+pz3)*(pz2+pz3);
	    if ( m23>0 ) m23=sqrt(m23);
	    
	    // Simple non-selecting way to use info of b-tag in W selection
	    // ------------------------------------------------------------
	    if ( fabs(m12-80.4)<fabs(m13-80.4) && fabs(m12-80.4)<fabs(m23-80.4) ) {
	      mw=m12;
	      if ( JHET[i]>loose_ || JHET[j]>loose_ ) bestWhasb=1;
	      if ( JHET[k]<=loose_ ) bestbnob=1;
	    }
	    if ( fabs(m13-80.4)<fabs(m12-80.4) && fabs(m13-80.4)<fabs(m23-80.4) ) { 
	      mw=m13;
	      if ( JHET[i]>loose_ || JHET[k]>loose_ ) bestWhasb=1;
	      if ( JHET[j]<=loose_ ) bestbnob=1;
	    }
	    if ( fabs(m23-80.4)<fabs(m13-80.4) && fabs(m23-80.4)<fabs(m12-80.4) ) {
	      mw=m23;
	      if ( JHET[j]>loose_ || JHET[k]>loose_ ) bestWhasb=1;
	      if ( JHET[i]<=loose_ ) bestbnob=1;
	    }
	    if ( fabs(m3-172)<fabs(m3best-172) ) {
	      m3best=m3;
	      mwbest=mw;
	      chi2=sqrt((m3best-172)*(m3best-172)/400+(mwbest-80.4)*(mwbest-80.4)/100)+
		3*bestWhasb+3*bestbnob;
	      it1=i;
	      it2=j;
	      it3=k;
	    }		       
	  }	    
	}
      }
    }

    // Search for a H signal in the remaining jets
    // -------------------------------------------
    if ( it1+it2+it3>0 ) {
      for ( int ii=0; ii<iJ-1 && ii<NHSJ-1; ii++ ) {
	if ( ii!=it1 && ii!=it2 && ii!=it3 ) {
	  for ( int jj=ii+1; jj<iJ && jj<NHSJ; jj++ ) { // limit to NHSJ the number of jets from h.s.
	    if ( jj!=it1 && jj!=it2 && jj!=it3 ) {
		
	      // Demand that both b-jets from H be tagged
	      // ----------------------------------------
	      if ( JHET[ii]>loose_ && JHET[jj]>loose_ ) {
		px4=JEt[ii]*cos(Jphi[ii]);
		px5=JEt[jj]*cos(Jphi[jj]);
		py4=JEt[ii]*sin(Jphi[ii]);
		py5=JEt[jj]*sin(Jphi[jj]);
		pz4=JEt[ii]/tan(2*atan(exp(-Jeta[ii])));
		pz5=JEt[jj]/tan(2*atan(exp(-Jeta[jj])));
		e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
		e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
		m45=(e4+e5)*(e4+e5)-(px4+px5)*(px4+px5)-(py4+py5)*(py4+py5)-(pz4+pz5)*(pz4+pz5);
		if ( m45>0 ) m45=sqrt(m45);
		if ( fabs(m45-120)<fabs(m45best-120) ) {
		  m45best=m45;
		  chi2ext=sqrt(chi2*chi2+(m45best-120)*(m45best-120)/200);
		  dpbb=3.145926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
		  ih1=ii;
		  ih2=jj;
		}
	      }
	    }
	  }
	}
      }

      // Compute mass of jets not selected in triplet of top and doublet of H decay
      // --------------------------------------------------------------------------
      if ( ih1+ih2>0 ) {
	spx=0;
	spy=0;
	spz=0;
	se=0;
	for ( int kk=0; kk<iJ && kk<NHSJ ; kk++ ) { // limit to NHSJ the number of jets from hard subprocess
	  if ( kk!= ih1 && kk!=ih2 && kk!=it1 && kk!=it2 && kk!=it3 ) { 
	    px=JEt[kk]*cos(Jphi[kk]);
	    py=JEt[kk]*sin(Jphi[kk]);
	    pz=JEt[kk]/tan(2*atan(exp(-Jeta[kk])));
	    spx+=px;
	    spy+=py;
	    spz+=pz;
	    se+=sqrt(px*px+py*py+pz*pz);
	  }
	  m_others=se*se-spx*spx-spy*spy-spz*spz;
	  if ( m_others>0 ) m_others=sqrt(m_others);
	}
	
	// Now compute mass and angle of two b-jets not chosen as H decay
	// --------------------------------------------------------------
	spx=0;
	spy=0;
	spz=0;
	se=0;
	double mbbnoh;
	for ( int k1=0; k1<iJ-1 && k1<NHSJ-1; k1++ ) {
	  if ( k1!=ih1 && k1!=ih2 && JHET[k1]>loose_ ) {
	    for ( int k2=k1+1; k2<iJ && k2<NHSJ; k2++ ) { // limit to 8 the number of jets from h.s.
	      if ( k2!=ih1 && k2!=ih2 && JHET[k2]>loose_ ) {
		px1=JEt[k1]*cos(Jphi[k1]);
		py1=JEt[k1]*sin(Jphi[k1]);
		pz1=JEt[k1]/tan(2*atan(exp(-Jeta[k1])));
		px2=JEt[k2]*cos(Jphi[k2]);
		py2=JEt[k2]*sin(Jphi[k2]);
		pz2=JEt[k2]/tan(2*atan(exp(-Jeta[k2])));
		spx=px1+px2;
		spy=py1+py2;
		spz=pz1+pz2;
		se=sqrt(px1*px1+py1*py1+pz1*pz1)+sqrt(px2*px2+py2*py2+pz2*pz2);
	      }
	      mbbnoh=se*se-spx*spx-spy*spy-spz*spz;
	      if ( mbbnoh>0 ) mbbnoh=sqrt(mbbnoh);
	      if ( mbbnoh>mbbnohmax ) {
		mbbnohmax=mbbnoh;
		dpbbnohmax=3.1415926-fabs(fabs(Jphi[k1]-Jphi[k2])-3.1415926);
	      }
	    }
	  }      
	}
      }      
    }      

    // Now compute mass of first 8 jets and their centrality
    // -----------------------------------------------------
    spx=0;
    spy=0;
    spz=0;
    se=0;
    set=0;
    for ( int k=0; k<iJ && k<NHSJ; k++ ) {
      px=JEt[k]*cos(Jphi[k]);
      py=JEt[k]*sin(Jphi[k]);
      pz=JEt[k]/tan(2*atan(exp(-Jeta[k])));
      spx+=px;
      spy+=py;
      spz+=pz;
      se+=sqrt(px*px+py*py+pz*pz);
      set+=JEt[k];
    }
    m8=se*se-spx*spx-spy*spy-spz*spz;
    if ( m8>0 ) { 
      m8=sqrt(m8);
      c8=set/m8;
    }

  } // end if caloJETS.size()
  else {
    std::cout << "ATTENTION: Jet collection empty" << std::endl;
  }



//   // Take genParticleCandidates
//   edm::Handle < MCParticleCollection > MCpartons;
//   iEvent.getByLabel( MCParticleLabel_, MCpartons );

//   // Take the mothers
//   MCParticleCollection::const_iterator MCp = MCpartons->begin();

//   //  int Zdaughters=0;
//   //  int Zcounter=0;
//   //  vector<int> HiggsDau;
//   for ( ; MCp != MCpartons->end(); ++MCp ) {
//     //    if(MCp->mPid()==25) HiggsDau.push_back(MCp->pid());
//     //    if(MCp->mPid()==23) ++Zdaughters;
//     //    if(MCp->pid()==23)  ++Zcounter;

// #ifdef DEBUG
//     std::cout << "For parton number = " << counter << std::endl;
//     std::cout << "status = " << MCp->status() << std::endl;
//     std::cout << "pdgId = " << MCp->pdgId() << std::endl;
//     std::cout << "Et = " << MCp->et() << std::endl;
//     std::cout << "Eta = " << MCp->eta() << std::endl;
//     std::cout << "Phi = " << MCp->phi() << std::endl;
//     std::cout << "Number of mothers = " << MCp->numberOfMothers() << std::endl;
//     std::cout << "first mother = " << MCp->mother() << std::endl;
//     std::cout << "Mother pdgId = " << MCp->mother()->pdgId() << std::endl;
// #endif // DEBUG
//   }
//   // CHECK for Higgs BR: Higgs inclusive
//   //  if(Zcounter!=0 || Zdaughters !=0) {
//   //    cout << "number of Z mother: " << Zcounter << " <--> number of Z daughters: " << Zdaughters << endl;
//   //    int iter=0;
//   //    for(vector<int>::iterator HiggsDau_it = HiggsDau.begin(); HiggsDau_it != HiggsDau.end(); ++HiggsDau_it, ++iter)
//   //      cout << "HiggsDaughters[" << iter << "]: " << *HiggsDau_it << endl;
//   //  }


  // Pixel jets
  // ----------
  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );
  const SimplePixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // Count the number of pixeljets with at leas numTkCut_ tracks
  // -----------------------------------------------------------
  SimplePixelJetCollection::const_iterator pj_it = pixeljets.begin();
  int goodPjNum = 0;
  for ( ; pj_it != pixeljets.end(); ++pj_it ) {
    if ( pj_it->tkNum() >= int(numTkCut_) ) ++goodPjNum;
  }

  // Pixel trigger requiring at least 2 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)
  // -----------------------------------------------------------------------------
  
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
  // All the jets together for the L1Trigger
  // ---------------------------------------
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
  // --------
  for ( BaseJetCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
    vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
  }
  
  // Multijet
  // --------
  
  // Central
  // -------
  bool response_cen = false;
  L1Trigger.Fill( vec_TriggerCenJet );
  response_cen = L1Trigger.Response();
  
  // Forward
  // -------
  bool response_for = false;
  L1Trigger.Fill( vec_TriggerForJet );
  response_for = L1Trigger.Response();
  
  // Tau
  // ---
  bool response_tau = false;
  L1Trigger.Fill( vec_TriggerTauJet );
  response_tau = L1Trigger.Response();
  
  // Full and no-forward
  // -------------------
  bool response = ( response_cen || response_tau || response_for );
  bool response_nofor = ( response_cen || response_tau );
  
  // MEt + Jet
  // ---------
  // Central
  // -------
  bool response_MEtJet_cen = false;
  sort( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
  reverse( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
  if ( vec_TriggerCenJet.size() != 0 ) {
    if ( (vec_TriggerCenJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_cen = true;
    }
  }
  // Tau
  // ---
  bool response_MEtJet_tau = false;
  sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  if ( vec_TriggerTauJet.size() != 0 ) {
    if ( (vec_TriggerTauJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_tau = true;
    }
  }
  // Forward
  // -------
  bool response_MEtJet_for = false;
  sort( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
  reverse( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
  if ( vec_TriggerForJet.size() != 0 ) {
    if ( (vec_TriggerForJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_for = true;
    }
  }
  // Full and no-forward
  // -------------------
  bool response_MEtJet = ( response_MEtJet_cen || response_MEtJet_tau || response_MEtJet_for );
  bool response_MEtJet_nofor = ( response_MEtJet_cen || response_MEtJet_tau );
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif


  /////////////////////////////////////////////////////////////////////////////////////
  // Now fill histograms
  // -------------------
  // Offline cuts
  // ------------
  if ( iJ >= 5 ) NJetsCut = true;
  if ( caloMET->mEtSig()>3 ) MEtSigCut = true;
  if ( NTags>=4 ) NTagsCut = true;
  
  // Fill kinematical histograms
  // ---------------------------
  double UncorrMEt =sqrt(UncorrSumEx*UncorrSumEx+UncorrSumEy*UncorrSumEy);
  double CorrMEt = sqrt(CorrSumEx*CorrSumEx+CorrSumEy*CorrSumEy);
  UncorrHt+=UncorrMEt;
  CorrHt+=CorrMEt;
  GoodHt+=CorrMEt;
  GoodHt2+=caloMET->et();
  
  double UncorrMEtSig=0;
  if ( UncorrSumEt>0 ) UncorrMEtSig = UncorrMEt/sqrt(UncorrSumEt);
  double CorrMEtSig = 0;
  if ( CorrSumEt>0 ) CorrMEtSig = CorrMEt/sqrt(CorrSumEt);  
  
  // Missing Et significance computed with tuned resolution function
  // ---------------------------------------------------------------
  double metsignew=0;
  if ( caloMET->sumEt()>0 ) metsignew = caloMET->et()/(sqrt(2)*0.8033*pow(caloMET->sumEt(),0.5004));

  UncorrMEt_SumEt_->Fill(caloMET->sumEt(),UncorrMEt);
  CorrMEt_SumEt_->Fill(caloMET->sumEt(),CorrMEt);
  MEt_SumEt_->Fill(caloMET->sumEt(),caloMET->et());
  UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt);
  CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt);
  MEt_SumEtC_->Fill(CorrSumEt,caloMET->et());
  UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt);
  CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt);
  MEt_SumEtJ_->Fill(GoodSumEt,caloMET->et());

  // Apply trigger requirement
  // -------------------------
  if ( response && NJetsCut && MEtSigCut && NTagsL>=2 ) {
    NJets_->Fill(goodIc5Jets);
    UncorrHt_->Fill(UncorrHt);
    CorrHt_->Fill(CorrHt);
    GoodHt_->Fill(GoodHt);
    GoodHt2_->Fill(GoodHt2);
    UncorrSumEt_->Fill(UncorrSumEt);
    CorrSumEt_->Fill(CorrSumEt);
    GoodSumEt_->Fill(GoodSumEt);
    MEt_->Fill(caloMET->et());
    MEtSig_->Fill(caloMET->mEtSig());
    MEtSigNew_->Fill(metsignew);
    MEtDPM_->Fill(dpmin);
    MEtDP1_->Fill(dp1st);
    MEtDP2_->Fill(dp2nd);
    MEtDP3_->Fill(dp3rd);
    UncorrMEtSig_->Fill(UncorrMEtSig);
    CorrMEtSig_->Fill(CorrMEtSig);
    M3best_->Fill(m3best);
    Mwbest_->Fill(mwbest);
    Chi2mass_->Fill(chi2);
    M45best_->Fill(m45best);
    Chi2ext_->Fill(chi2ext);
    MEx_SumEt_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()));
    MEx_SumEt_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()));
    DP12_->Fill(dp12);
    DPbb_->Fill(dpbb);
    M_others_->Fill(m_others);
    Mbbnoh_->Fill(mbbnohmax);
    DPbbnoh_->Fill(dpbbnohmax);
    M8_->Fill(m8);
    C8_->Fill(c8);
    
    if ( NTagsL>=3 ) {
      NJetsS_->Fill(iJ);
      UncorrHtS_->Fill(UncorrHt);
      CorrHtS_->Fill(CorrHt);
      GoodHtS_->Fill(GoodHt);
      GoodHt2S_->Fill(GoodHt2);
      UncorrSumEtS_->Fill(UncorrSumEt);
      CorrSumEtS_->Fill(CorrSumEt);
      GoodSumEtS_->Fill(GoodSumEt);
      MEtS_->Fill(caloMET->et());
      MEtSigS_->Fill(caloMET->mEtSig());
      MEtSigNewS_->Fill(metsignew);
      MEtDPMS_->Fill(dpmin);
      MEtDP1S_->Fill(dp1st);
      MEtDP2S_->Fill(dp2nd);
      MEtDP3S_->Fill(dp3rd);
      UncorrMEtSigS_->Fill(UncorrMEtSig);
      CorrMEtSigS_->Fill(CorrMEtSig);
      M3bestS_->Fill(m3best);
      MwbestS_->Fill(mwbest);
      Chi2massS_->Fill(chi2);
      M45bestS_->Fill(m45best);
      Chi2extS_->Fill(chi2ext);
      MEx_SumEtS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()));
      MEx_SumEtS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()));
      DP12S_->Fill(dp12);
      DPbbS_->Fill(dpbb);
      M_othersS_->Fill(m_others);
      MbbnohS_->Fill(mbbnohmax);
      DPbbnohS_->Fill(dpbbnohmax);
      M8S_->Fill(m8);
      C8S_->Fill(c8);
    }
      
    if ( NTagsCut ) {
      NJetsSS_->Fill(iJ);
      UncorrHtSS_->Fill(UncorrHt);
      CorrHtSS_->Fill(CorrHt);
      GoodHtSS_->Fill(GoodHt);
      GoodHt2SS_->Fill(GoodHt2);
      UncorrSumEtSS_->Fill(UncorrSumEt);
      CorrSumEtSS_->Fill(CorrSumEt);
      GoodSumEtSS_->Fill(GoodSumEt);
      MEtSS_->Fill(caloMET->et());
      MEtSigSS_->Fill(caloMET->mEtSig());
      MEtSigNewSS_->Fill(metsignew);
      MEtDPMSS_->Fill(dpmin);
      MEtDP1SS_->Fill(dp1st);
      MEtDP2SS_->Fill(dp2nd);
      MEtDP3SS_->Fill(dp3rd);
      UncorrMEtSigSS_->Fill(UncorrMEtSig);
      CorrMEtSigSS_->Fill(CorrMEtSig);
      M3bestSS_->Fill(m3best);
      MwbestSS_->Fill(mwbest);
      Chi2massSS_->Fill(chi2);
      M45bestSS_->Fill(m45best);
      Chi2extSS_->Fill(chi2ext);
      MEx_SumEtSS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()));
      MEx_SumEtSS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()));
      DP12SS_->Fill(dp12);
      DPbbSS_->Fill(dpbb);
      M_othersSS_->Fill(m_others);      
      MbbnohSS_->Fill(mbbnohmax);
      DPbbnohSS_->Fill(dpbbnohmax);
      M8SS_->Fill(m8);
      C8SS_->Fill(c8);
    }

    if ( NTagsT>=4 ) {
      NJetsSSS_->Fill(iJ);
      UncorrHtSSS_->Fill(UncorrHt);
      CorrHtSSS_->Fill(CorrHt);
      GoodHtSSS_->Fill(GoodHt);
      GoodHt2SSS_->Fill(GoodHt2);
      UncorrSumEtSSS_->Fill(UncorrSumEt);
      CorrSumEtSSS_->Fill(CorrSumEt);
      GoodSumEtSSS_->Fill(GoodSumEt);
      MEtSSS_->Fill(caloMET->et());
      MEtSigSSS_->Fill(caloMET->mEtSig());
      MEtSigNewSSS_->Fill(metsignew);
      MEtDPMSSS_->Fill(dpmin);
      MEtDP1SSS_->Fill(dp1st);
      MEtDP2SSS_->Fill(dp2nd);
      MEtDP3SSS_->Fill(dp3rd);
      UncorrMEtSigSSS_->Fill(UncorrMEtSig);
      CorrMEtSigSSS_->Fill(CorrMEtSig);
      M3bestSSS_->Fill(m3best);
      MwbestSSS_->Fill(mwbest);
      Chi2massSSS_->Fill(chi2);
      M45bestSSS_->Fill(m45best);
      Chi2extSSS_->Fill(chi2ext);
      MEx_SumEtSSS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()));
      MEx_SumEtSSS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()));
      DP12SSS_->Fill(dp12);
      DPbbSSS_->Fill(dpbb);
      M_othersSSS_->Fill(m_others);      
      MbbnohSSS_->Fill(mbbnohmax);
      DPbbnohSSS_->Fill(dpbbnohmax);
      M8SSS_->Fill(m8);
      C8SSS_->Fill(c8);
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
// --------------------------------------------------------------------------------------
void PTag::beginJob(const edm::EventSetup&) {
}


// ------------ method called once each job just after ending the event loop  ------------
// ---------------------------------------------------------------------------------------
void PTag::endJob() {

  Multi_Vertex_Dz_->Write();

  // Histograms for events passing trigger
  // -------------------------------------
  NJets_->Write();
  UncorrHt_->Write();
  CorrHt_->Write();
  GoodHt_->Write();
  GoodHt2_->Write();
  UncorrSumEt_->Write();
  CorrSumEt_->Write();
  GoodSumEt_->Write();
  MEt_->Write();
  MEtSig_->Write();
  MEtSigNew_->Write();
  MEtDPM_->Write();
  MEtDP1_->Write();
  MEtDP2_->Write();
  MEtDP3_->Write();
  UncorrMEtSig_->Write();
  CorrMEtSig_->Write();
  M3best_->Write();
  Mwbest_->Write();
  Chi2mass_->Write();
  M45best_->Write();
  Chi2ext_->Write();
  MEx_SumEt_->Write();
  DP12_->Write();
  DPbb_->Write();
  M_others_->Write();
  Mbbnoh_->Write();
  DPbbnoh_->Write();
  M8_->Write();
  C8_->Write();
  
  // Histograms for events passing trigger and NJet cut
  // --------------------------------------------------
  NJetsS_->Write();
  UncorrHtS_->Write();
  CorrHtS_->Write();
  GoodHtS_->Write();
  GoodHt2S_->Write();
  UncorrSumEtS_->Write();
  CorrSumEtS_->Write();
  GoodSumEtS_->Write();
  MEtS_->Write();
  MEtSigS_->Write();
  MEtSigNewS_->Write();
  MEtDPMS_->Write();
  MEtDP1S_->Write();
  MEtDP2S_->Write();
  MEtDP3S_->Write();
  UncorrMEtSigS_->Write();
  CorrMEtSigS_->Write();
  M3bestS_->Write();
  MwbestS_->Write();
  Chi2massS_->Write();
  M45bestS_->Write();
  Chi2extS_->Write();
  MEx_SumEtS_->Write();
  DP12S_->Write();
  DPbbS_->Write();
  M_othersS_->Write();
  MbbnohS_->Write();
  DPbbnohS_->Write();
  M8S_->Write();
  C8S_->Write();

  // Histograms for events passing trigger, NJet, and METS cut
  // ---------------------------------------------------------
  NJetsSS_->Write();
  UncorrHtSS_->Write();
  CorrHtSS_->Write();
  GoodHtSS_->Write();
  GoodHt2SS_->Write();
  UncorrSumEtSS_->Write();
  CorrSumEtSS_->Write();
  GoodSumEtSS_->Write();
  MEtSS_->Write();
  MEtSigSS_->Write();
  MEtSigNewSS_->Write();
  MEtDPMSS_->Write();
  MEtDP1SS_->Write();
  MEtDP2SS_->Write();
  MEtDP3SS_->Write();
  UncorrMEtSigSS_->Write();
  CorrMEtSigSS_->Write();
  M3bestSS_->Write();
  MwbestSS_->Write();
  Chi2massSS_->Write();
  M45bestSS_->Write();
  Chi2extSS_->Write();
  MEx_SumEtSS_->Write();
  DP12SS_->Write();
  DPbbSS_->Write();
  M_othersSS_->Write();
  MbbnohSS_->Write();
  DPbbnohSS_->Write();
  M8SS_->Write();
  C8SS_->Write();

  // Histograms for events passing trigger, NJet, and METS cut
  // ---------------------------------------------------------
  NJetsSSS_->Write();
  UncorrHtSSS_->Write();
  CorrHtSSS_->Write();
  GoodHtSSS_->Write();
  GoodHt2SSS_->Write();
  UncorrSumEtSSS_->Write();
  CorrSumEtSSS_->Write();
  GoodSumEtSSS_->Write();
  MEtSSS_->Write();
  MEtSigSSS_->Write();
  MEtSigNewSSS_->Write();
  MEtDPMSSS_->Write();
  MEtDP1SSS_->Write();
  MEtDP2SSS_->Write();
  MEtDP3SSS_->Write();
  UncorrMEtSigSSS_->Write();
  CorrMEtSigSSS_->Write();
  M3bestSSS_->Write();
  MwbestSSS_->Write();
  Chi2massSSS_->Write();
  M45bestSSS_->Write();
  Chi2extSSS_->Write();
  MEx_SumEtSSS_->Write();
  DP12SSS_->Write();
  DPbbSSS_->Write();
  M_othersSSS_->Write();
  MbbnohSSS_->Write();
  DPbbnohSSS_->Write();
  M8SSS_->Write();
  C8SSS_->Write();
  
  // Profile plots of MET vs SumEt
  // -----------------------------
  UncorrMEt_SumEt_->Write();
  CorrMEt_SumEt_->Write();
  MEt_SumEt_->Write();
  UncorrMEt_SumEtC_->Write();
  CorrMEt_SumEtC_->Write();
  MEt_SumEtC_->Write();
  UncorrMEt_SumEtJ_->Write();
  CorrMEt_SumEtJ_->Write();
  MEt_SumEtJ_->Write();

}

// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(PTag);
