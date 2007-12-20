// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/TDAna.h"

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
#include <string>

// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------

L1Trig TDAna::L1Trigger;

// Constructors and destructor
// ---------------------------
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
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) ),
  QCD_( iConfig.getUntrackedParameter<bool> ( "QCD" ) )
{

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd();

  // Histograms
  // ----------
  NJets_ = new TH1D ( "NJets", "Number of selected jets", 20, 0, 20 );
  UncorrHt_ = new TH1D ( "UncorrHt", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHt_ = new TH1D ( "CorrHt", "Ht with corrected jets", 50, 0, 4000 );
  GoodHt_ = new TH1D ( "GoodHt", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2_ = new TH1D ( "GoodHt2", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEt_ = new TH1D ( "UncorrSumEt", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEt_ = new TH1D ( "CorrSumEt", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEt_ = new TH1D ( "GoodSumEt", "SumEt with selected jets", 50, 0, 4000 );
  MEt_ = new TH1D ( "MEt", "MEt", 50, 0, 500 );
  MEtSig_ = new TH1D ( "MEtSig", "MEtSig", 50, 0, 20 );
  MEtSigNew_ = new TH1D ( "MEtSigNew", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPM_ = new TH1D ( "MEtDPM", "MEtDPM", 50, 0, 3.2 );
  MEtDP1_ = new TH1D ( "MEtDP1", "MEtDP1", 50, 0, 3.2 );
  MEtDP2_ = new TH1D ( "MEtDP2", "MEtDP2", 50, 0, 3.2 );
  MEtDP3_ = new TH1D ( "MEtDP3", "MEtDP3", 50, 0, 3.2 );
  UncorrMEtSig_ = new TH1D ( "UncorrMEtSig", "MEtSig", 50, 0, 20 );
  CorrMEtSig_ = new TH1D ( "CorrMEtSig", "MEtSig", 50, 0, 20 );
  M3best_ = new TH1D ( "M3best", "M3best", 50, 0, 300 );
  Mwbest_ = new TH1D ( "Mwbest", "Mwbest", 50, 0, 200 );
  Chi2mass_ = new TH1D ( "Chi2mass", "Chi2mass", 50, 0, 50 );
  M45best_ = new TH1D ( "M45best", "M45best", 50, 0, 300 );
  Chi2ext_ = new TH1D ( "Chi2ext", "Chi2ext", 50, 0, 50 );
  MEx_SumEt_ = new TH2D ( "MEx_SumEt", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12_ = new TH1D ( "DP12", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbb_ = new TH1D ( "DPbb", "Delta phi bb", 50, 0., 3.2 );
  M_others_ = new TH1D ( "M_others", "Mass of other jets", 50, 0., 1500. );
  Mbbnoh_ = new TH1D ( "Mbbnoh", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnoh_ = new TH1D ( "DPbbnoh", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8_ = new TH1D ( "M8", "Mass 8 jets", 50, 0., 2000.) ;
  C8_ = new TH1D ( "C8", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestall_ = new TH1D ( "M45bestall", "M45best", 50, 0, 300 );
  Chi2extall_ = new TH1D ( "Chi2extall", "Chi2ext", 50, 0, 50 );
  DPbball_ = new TH1D ( "DPbball", "Delta phi bb", 50, 0., 3.2 );

  NJetsS_ = new TH1D ( "NJetsS", "Number of selected jets", 20, 0, 20 );
  UncorrHtS_ = new TH1D ( "UncorrHtS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtS_ = new TH1D ( "CorrHtS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtS_ = new TH1D ( "GoodHtS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2S_ = new TH1D ( "GoodHt2S", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtS_ = new TH1D ( "UncorrSumEtS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtS_ = new TH1D ( "CorrSumEtS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtS_ = new TH1D ( "GoodSumEtS", "SumEt with selected jets", 50, 0, 4000 );
  MEtS_ = new TH1D ( "MEtS", "MEt", 50, 0, 500 );
  MEtSigS_ = new TH1D ( "MEtSigS", "MEtSig", 50, 0, 20 );
  MEtSigNewS_ = new TH1D ( "MEtSigNewS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMS_ = new TH1D ( "MEtDPMS", "MEtDPM", 50, 0, 3.2 );
  MEtDP1S_ = new TH1D ( "MEtDP1S", "MEtDP1S", 50, 0, 3.2 );
  MEtDP2S_ = new TH1D ( "MEtDP2S", "MEtDP2S", 50, 0, 3.2 );
  MEtDP3S_ = new TH1D ( "MEtDP3S", "MEtDP3S", 50, 0, 3.2 );
  UncorrMEtSigS_ = new TH1D ( "UncorrMEtSigS", "MEtSig", 50, 0, 20 );
  CorrMEtSigS_ = new TH1D ( "CorrMEtSigS", "MEtSig", 50, 0, 20 );
  M3bestS_ = new TH1D ( "M3bestS", "M3bestS", 50, 0, 300 );
  MwbestS_ = new TH1D ( "MwbestS", "MwbestS", 50, 0, 200 );
  Chi2massS_ = new TH1D ( "Chi2massS", "Chi2massS", 50, 0, 50 );
  M45bestS_ = new TH1D ( "M45bestS", "M45best", 50, 0, 300 );
  Chi2extS_ = new TH1D ( "Chi2extS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtS_ = new TH2D ( "MEx_SumEtS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12S_ = new TH1D ( "DP12S", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbS_ = new TH1D ( "DPbbS", "Delta phi bb", 50, 0., 3.2 );
  M_othersS_ = new TH1D ( "M_othersS", "Mass of other jets", 50, 0., 1500. );
  MbbnohS_ = new TH1D ( "MbbnohS", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohS_ = new TH1D ( "DPbbnohS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8S_ = new TH1D ( "M8S", "Mass 8 jets", 50, 0., 2000.) ;
  C8S_ = new TH1D ( "C8S", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallS_ = new TH1D ( "M45bestallS", "M45best", 50, 0, 300 );
  Chi2extallS_ = new TH1D ( "Chi2extallS", "Chi2ext", 50, 0, 50 );
  DPbballS_ = new TH1D ( "DPbballS", "Delta phi bb", 50, 0., 3.2 );

  NJetsSS_ = new TH1D ( "NJetsSS", "Number of selected jets", 20, 0, 20 );
  UncorrHtSS_ = new TH1D ( "UncorrHtSS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSS_ = new TH1D ( "CorrHtSS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSS_ = new TH1D ( "GoodHtSS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SS_ = new TH1D ( "GoodHt2SS", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSS_ = new TH1D ( "UncorrSumEtSS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSS_ = new TH1D ( "CorrSumEtSS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSS_ = new TH1D ( "GoodSumEtSS", "SumEt with selected jets", 50, 0, 4000 );
  MEtSS_ = new TH1D ( "MEtSS", "MEt", 50, 0, 500 );
  MEtSigSS_ = new TH1D ( "MEtSigSS", "MEtSig", 50, 0, 20 );
  MEtSigNewSS_ = new TH1D ( "MEtSigNewSS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSS_ = new TH1D ( "MEtDPMSS", "MEtDPMSS", 50, 0, 3.2 );
  MEtDP1SS_ = new TH1D ( "MEtDP1SS", "MEtDP1SS", 50, 0, 3.2 );
  MEtDP2SS_ = new TH1D ( "MEtDP2SS", "MEtDP2SS", 50, 0, 3.2 );
  MEtDP3SS_ = new TH1D ( "MEtDP3SS", "MEtDP3SS", 50, 0, 3.2 );
  UncorrMEtSigSS_ = new TH1D ( "UncorrMEtSigSS", "UncorrMEtSigSS", 50, 0, 20 );
  CorrMEtSigSS_ = new TH1D ( "CorrMEtSigSS", "CorrMEtSigSS", 50, 0, 20 );
  M3bestSS_ = new TH1D ( "M3bestSS", "M3bestSS", 50, 0, 300 );
  MwbestSS_ = new TH1D ( "MwbestSS", "MwbestSS", 50, 0, 200 );
  Chi2massSS_ = new TH1D ( "Chi2massSS", "Chi2massSS", 50, 0, 50 );
  M45bestSS_ = new TH1D ( "M45bestSS", "M45best", 50, 0, 300 );
  Chi2extSS_ = new TH1D ( "Chi2extSS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSS_ = new TH2D ( "MEx_SumEtSS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SS_ = new TH1D ( "DP12SS", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSS_ = new TH1D ( "DPbbSS", "Delta phi bb", 50, 0., 3.2 );
  M_othersSS_ = new TH1D ( "M_othersSS", "Mass of other jets", 50, 0., 1500. );
  MbbnohSS_ = new TH1D ( "MbbnohSS", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSS_ = new TH1D ( "DPbbnohSS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SS_ = new TH1D ( "M8SS", "Mass 8 jets", 50, 0., 2000.) ;
  C8SS_ = new TH1D ( "C8SS", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSS_ = new TH1D ( "M45bestallSS", "M45best", 50, 0, 300 );
  Chi2extallSS_ = new TH1D ( "Chi2extallSS", "Chi2ext", 50, 0, 50 );
  DPbballSS_ = new TH1D ( "DPbballSS", "Delta phi bb", 50, 0., 3.2 );

  NJetsSSS_ = new TH1D ( "NJetsSSS", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSS_ = new TH1D ( "UncorrHtSSS", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSS_ = new TH1D ( "CorrHtSSS", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSS_ = new TH1D ( "GoodHtSSS", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSS_ = new TH1D ( "GoodHt2SSS", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSS_ = new TH1D ( "UncorrSumEtSSS", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSS_ = new TH1D ( "CorrSumEtSSS", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSS_ = new TH1D ( "GoodSumEtSSS", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSS_ = new TH1D ( "MEtSSS", "MEt", 50, 0, 500 );
  MEtSigSSS_ = new TH1D ( "MEtSigSSS", "MEtSig", 50, 0, 20 );
  MEtSigNewSSS_ = new TH1D ( "MEtSigNewSSS", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSS_ = new TH1D ( "MEtDPMSSS", "MEtDPMSSS", 50, 0, 3.2 );
  MEtDP1SSS_ = new TH1D ( "MEtDP1SSS", "MEtDP1SSS", 50, 0, 3.2 );
  MEtDP2SSS_ = new TH1D ( "MEtDP2SSS", "MEtDP2SSS", 50, 0, 3.2 );
  MEtDP3SSS_ = new TH1D ( "MEtDP3SSS", "MEtDP3SSS", 50, 0, 3.2 );
  UncorrMEtSigSSS_ = new TH1D ( "UncorrMEtSigSSS", "UncorrMEtSigSSS", 50, 0, 20 );
  CorrMEtSigSSS_ = new TH1D ( "CorrMEtSigSSS", "CorrMEtSigSSS", 50, 0, 20 );
  M3bestSSS_ = new TH1D ( "M3bestSSS", "M3bestSSS", 50, 0, 300 );
  MwbestSSS_ = new TH1D ( "MwbestSSS", "MwbestSSS", 50, 0, 200 );
  Chi2massSSS_ = new TH1D ( "Chi2massSSS", "Chi2massSSS", 50, 0, 50 );
  M45bestSSS_ = new TH1D ( "M45bestSSS", "M45best", 50, 0, 300 );
  Chi2extSSS_ = new TH1D ( "Chi2extSSS", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSS_ = new TH2D ( "MEx_SumEtSSS", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSS_ = new TH1D ( "DP12SSS", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSS_ = new TH1D ( "DPbbSSS", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSS_ = new TH1D ( "M_othersSSS", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSS_ = new TH1D ( "MbbnohSSS", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSSS_ = new TH1D ( "DPbbnohSSS", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SSS_ = new TH1D ( "M8SSS", "Mass 8 jets", 50, 0., 2000.) ;
  C8SSS_ = new TH1D ( "C8SSS", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSS_ = new TH1D ( "M45bestallSSS", "M45best", 50, 0, 300 );
  Chi2extallSSS_ = new TH1D ( "Chi2extallSSS", "Chi2ext", 50, 0, 50 );
  DPbballSSS_ = new TH1D ( "DPbballSSS", "Delta phi bb", 50, 0., 3.2 );
  N4NJSSS_ = new TH1D ( "N4NJSSS", "N of 4HEL tags vs N jets", 20, 0, 20 );
  E4NJSSS_ = new TH1D ( "E4NJSSS", "Efficiency of 4HEL tags vs N jets", 20, 0, 20 );

  NJetsW_ = new TH1D ( "NJetsW", "Number of selected jets", 20, 0, 20 );
  UncorrHtW_ = new TH1D ( "UncorrHtW", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtW_ = new TH1D ( "CorrHtW", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtW_ = new TH1D ( "GoodHtW", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2W_ = new TH1D ( "GoodHt2W", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtW_ = new TH1D ( "UncorrSumEtW", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtW_ = new TH1D ( "CorrSumEtW", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtW_ = new TH1D ( "GoodSumEtW", "SumEt with selected jets", 50, 0, 4000 );
  MEtW_ = new TH1D ( "MEtW", "MEt", 50, 0, 500 );
  MEtSigW_ = new TH1D ( "MEtSigW", "MEtSig", 50, 0, 20 );
  MEtSigNewW_ = new TH1D ( "MEtSigNewW", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMW_ = new TH1D ( "MEtDPMW", "MEtDPM", 50, 0, 3.2 );
  MEtDP1W_ = new TH1D ( "MEtDP1W", "MEtDP1", 50, 0, 3.2 );
  MEtDP2W_ = new TH1D ( "MEtDP2W", "MEtDP2", 50, 0, 3.2 );
  MEtDP3W_ = new TH1D ( "MEtDP3W", "MEtDP3", 50, 0, 3.2 );
  UncorrMEtSigW_ = new TH1D ( "UncorrMEtSigW", "MEtSig", 50, 0, 20 );
  CorrMEtSigW_ = new TH1D ( "CorrMEtSigW", "MEtSig", 50, 0, 20 );
  M3bestW_ = new TH1D ( "M3bestW", "M3best", 50, 0, 300 );
  MwbestW_ = new TH1D ( "MwbestW", "Mwbest", 50, 0, 200 );
  Chi2massW_ = new TH1D ( "Chi2massW", "Chi2mass", 50, 0, 50 );
  M45bestW_ = new TH1D ( "M45bestW", "M45best", 50, 0, 300 );
  Chi2extW_ = new TH1D ( "Chi2extW", "Chi2ext", 50, 0, 50 );
  MEx_SumEtW_ = new TH2D ( "MEx_SumEtW", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12W_ = new TH1D ( "DP12W", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbW_ = new TH1D ( "DPbbW", "Delta phi bb", 50, 0., 3.2 );
  M_othersW_ = new TH1D ( "M_othersW", "Mass of other jets", 50, 0., 1500. );
  MbbnohW_ = new TH1D ( "MbbnohW", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohW_ = new TH1D ( "DPbbnohW", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8W_ = new TH1D ( "M8W", "Mass 8 jets", 50, 0., 2000.) ;
  C8W_ = new TH1D ( "C8W", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallW_ = new TH1D ( "M45bestallW", "M45best", 50, 0, 300 );
  Chi2extallW_ = new TH1D ( "Chi2extallW", "Chi2ext", 50, 0, 50 );
  DPbballW_ = new TH1D ( "DPbballW", "Delta phi bb", 50, 0., 3.2 );

  NJetsSW_ = new TH1D ( "NJetsSW", "Number of selected jets", 20, 0, 20 );
  UncorrHtSW_ = new TH1D ( "UncorrHtSW", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSW_ = new TH1D ( "CorrHtSW", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSW_ = new TH1D ( "GoodHtSW", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SW_ = new TH1D ( "GoodHt2SW", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSW_ = new TH1D ( "UncorrSumEtSW", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSW_ = new TH1D ( "CorrSumEtSW", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSW_ = new TH1D ( "GoodSumEtSW", "SumEt with selected jets", 50, 0, 4000 );
  MEtSW_ = new TH1D ( "MEtSW", "MEt", 50, 0, 500 );
  MEtSigSW_ = new TH1D ( "MEtSigSW", "MEtSig", 50, 0, 20 );
  MEtSigNewSW_ = new TH1D ( "MEtSigNewSW", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSW_ = new TH1D ( "MEtDPMSW", "MEtDPM", 50, 0, 3.2 );
  MEtDP1SW_ = new TH1D ( "MEtDP1SW", "MEtDP1S", 50, 0, 3.2 );
  MEtDP2SW_ = new TH1D ( "MEtDP2SW", "MEtDP2S", 50, 0, 3.2 );
  MEtDP3SW_ = new TH1D ( "MEtDP3SW", "MEtDP3S", 50, 0, 3.2 );
  UncorrMEtSigSW_ = new TH1D ( "UncorrMEtSigSW", "MEtSig", 50, 0, 20 );
  CorrMEtSigSW_ = new TH1D ( "CorrMEtSigSW", "MEtSig", 50, 0, 20 );
  M3bestSW_ = new TH1D ( "M3bestSW", "M3bestS", 50, 0, 300 );
  MwbestSW_ = new TH1D ( "MwbestSW", "MwbestS", 50, 0, 200 );
  Chi2massSW_ = new TH1D ( "Chi2massSW", "Chi2massS", 50, 0, 50 );
  M45bestSW_ = new TH1D ( "M45bestSW", "M45best", 50, 0, 300 );
  Chi2extSW_ = new TH1D ( "Chi2extSW", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSW_ = new TH2D ( "MEx_SumEtSW", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SW_ = new TH1D ( "DP12SW", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSW_ = new TH1D ( "DPbbSW", "Delta phi bb", 50, 0., 3.2 );
  M_othersSW_ = new TH1D ( "M_othersSW", "Mass of other jets", 50, 0., 1500. );
  MbbnohSW_ = new TH1D ( "MbbnohSW", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSW_ = new TH1D ( "DPbbnohSW", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SW_ = new TH1D ( "M8SW", "Mass 8 jets", 50, 0., 2000.) ;
  C8SW_ = new TH1D ( "C8SW", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSW_ = new TH1D ( "M45bestallSW", "M45best", 50, 0, 300 );
  Chi2extallSW_ = new TH1D ( "Chi2extallSW", "Chi2ext", 50, 0, 50 );
  DPbballSW_ = new TH1D ( "DPbballSW", "Delta phi bb", 50, 0., 3.2 );

  NJetsSSW_ = new TH1D ( "NJetsSSW", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSW_ = new TH1D ( "UncorrHtSSW", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSW_ = new TH1D ( "CorrHtSSW", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSW_ = new TH1D ( "GoodHtSSW", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSW_ = new TH1D ( "GoodHt2SSW", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSW_ = new TH1D ( "UncorrSumEtSSW", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSW_ = new TH1D ( "CorrSumEtSSW", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSW_ = new TH1D ( "GoodSumEtSSW", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSW_ = new TH1D ( "MEtSSW", "MEt", 50, 0, 500 );
  MEtSigSSW_ = new TH1D ( "MEtSigSSW", "MEtSig", 50, 0, 20 );
  MEtSigNewSSW_ = new TH1D ( "MEtSigNewSSW", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSW_ = new TH1D ( "MEtDPMSSW", "MEtDPMSS", 50, 0, 3.2 );
  MEtDP1SSW_ = new TH1D ( "MEtDP1SSW", "MEtDP1SS", 50, 0, 3.2 );
  MEtDP2SSW_ = new TH1D ( "MEtDP2SSW", "MEtDP2SS", 50, 0, 3.2 );
  MEtDP3SSW_ = new TH1D ( "MEtDP3SSW", "MEtDP3SS", 50, 0, 3.2 );
  UncorrMEtSigSSW_ = new TH1D ( "UncorrMEtSigSSW", "UncorrMEtSigSS", 50, 0, 20 );
  CorrMEtSigSSW_ = new TH1D ( "CorrMEtSigSSW", "CorrMEtSigSS", 50, 0, 20 );
  M3bestSSW_ = new TH1D ( "M3bestSSW", "M3bestSS", 50, 0, 300 );
  MwbestSSW_ = new TH1D ( "MwbestSSW", "MwbestSS", 50, 0, 200 );
  Chi2massSSW_ = new TH1D ( "Chi2massSSW", "Chi2massSS", 50, 0, 50 );
  M45bestSSW_ = new TH1D ( "M45bestSSW", "M45best", 50, 0, 300 );
  Chi2extSSW_ = new TH1D ( "Chi2extSSW", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSW_ = new TH2D ( "MEx_SumEtSSW", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSW_ = new TH1D ( "DP12SSW", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSW_ = new TH1D ( "DPbbSSW", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSW_ = new TH1D ( "M_othersSSW", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSW_ = new TH1D ( "MbbnohSSW", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSSW_ = new TH1D ( "DPbbnohSSW", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SSW_ = new TH1D ( "M8SSW", "Mass 8 jets", 50, 0., 2000.) ;
  C8SSW_ = new TH1D ( "C8SSW", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSW_ = new TH1D ( "M45bestallSSW", "M45best", 50, 0, 300 );
  Chi2extallSSW_ = new TH1D ( "Chi2extallSSW", "Chi2ext", 50, 0, 50 );
  DPbballSSW_ = new TH1D ( "DPbballSSW", "Delta phi bb", 50, 0., 3.2 );

  NJetsSSSW_ = new TH1D ( "NJetsSSSW", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSSW_ = new TH1D ( "UncorrHtSSSW", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSSW_ = new TH1D ( "CorrHtSSSW", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSSW_ = new TH1D ( "GoodHtSSSW", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSSW_ = new TH1D ( "GoodHt2SSSW", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSSW_ = new TH1D ( "UncorrSumEtSSSW", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSSW_ = new TH1D ( "CorrSumEtSSSW", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSSW_ = new TH1D ( "GoodSumEtSSSW", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSSW_ = new TH1D ( "MEtSSSW", "MEt", 50, 0, 500 );
  MEtSigSSSW_ = new TH1D ( "MEtSigSSSW", "MEtSig", 50, 0, 20 );
  MEtSigNewSSSW_ = new TH1D ( "MEtSigNewSSSW", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSSW_ = new TH1D ( "MEtDPMSSSW", "MEtDPMSSS", 50, 0, 3.2 );
  MEtDP1SSSW_ = new TH1D ( "MEtDP1SSSW", "MEtDP1SSS", 50, 0, 3.2 );
  MEtDP2SSSW_ = new TH1D ( "MEtDP2SSSW", "MEtDP2SSS", 50, 0, 3.2 );
  MEtDP3SSSW_ = new TH1D ( "MEtDP3SSSW", "MEtDP3SSS", 50, 0, 3.2 );
  UncorrMEtSigSSSW_ = new TH1D ( "UncorrMEtSigSSSW", "UncorrMEtSigSSS", 50, 0, 20 );
  CorrMEtSigSSSW_ = new TH1D ( "CorrMEtSigSSSW", "CorrMEtSigSSS", 50, 0, 20 );
  M3bestSSSW_ = new TH1D ( "M3bestSSSW", "M3bestSSS", 50, 0, 300 );
  MwbestSSSW_ = new TH1D ( "MwbestSSSW", "MwbestSSS", 50, 0, 200 );
  Chi2massSSSW_ = new TH1D ( "Chi2massSSSW", "Chi2massSSS", 50, 0, 50 );
  M45bestSSSW_ = new TH1D ( "M45bestSSSW", "M45best", 50, 0, 300 );
  Chi2extSSSW_ = new TH1D ( "Chi2extSSSW", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSSW_ = new TH2D ( "MEx_SumEtSSSW", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSSW_ = new TH1D ( "DP12SSSW", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSSW_ = new TH1D ( "DPbbSSSW", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSSW_ = new TH1D ( "M_othersSSSW", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSSW_ = new TH1D ( "MbbnohSSSW", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSSSW_ = new TH1D ( "DPbbnohSSSW", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M8SSSW_ = new TH1D ( "M8SSSW", "Mass 8 jets", 50, 0., 2000.) ;
  C8SSSW_ = new TH1D ( "C8SSSW", "Centrality 8 jets", 50, 0., 1.2 );
  N4NJSSSW_ = new TH1D ( "N4NJSSSW", "N of 4HEL tags vs N jets", 20, 0, 20 );
  E4NJSSSW_ = new TH1D ( "E4NJSSSW", "Efficiency of 4HEL tags vs N jets", 20, 0, 20 );
  M45bestallSSSW_ = new TH1D ( "M45bestallSSSW", "M45best", 50, 0, 300 );
  Chi2extallSSSW_ = new TH1D ( "Chi2extallSSSW", "Chi2ext", 50, 0, 50 );
  DPbballSSSW_ = new TH1D ( "DPbballSSSW", "Delta phi bb", 50, 0., 3.2 );

  UncorrMEt_SumEt_ = new TH2D ( "UncorrMEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEt_ = new TH2D ( "CorrMEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEt_ = new TH2D ( "MEt_SumEt", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  UncorrMEt_SumEtC_ = new TH2D ( "UncorrMEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEtC_ = new TH2D ( "CorrMEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEtC_ = new TH2D ( "MEt_SumEtC", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  UncorrMEt_SumEtJ_ = new TH2D ( "UncorrMEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  CorrMEt_SumEtJ_ = new TH2D ( "CorrMEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  MEt_SumEtJ_ = new TH2D ( "MEt_SumEtJ", "MEt vs SumEt", 100, 0, 4000, 100, 0, 1000 );
  
  loose_  = 2.3; // high eff -> 70.49% b / 32.33% c / 8.64% uds / 10.43% g / 9.98% udsg // P.Schilling 23/10/07
  medium_ = 5.3; // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
  tight_  = 4.8; // high pur -> 31.94% b /  2.93% c / 0.10% uds /  0.11% g / 0.10% udsg // P.Schilling 23/10/07

  // Read PTag file
  // --------------
  for ( int i=0; i<1000; i++ ) {
    PHETL[i]=0;
    PHETM[i]=0;
    PHETT[i]=0;
    PHPTL[i]=0;
    PHPTM[i]=0;
    PHPTT[i]=0;
  }
  ifstream PTagFile ("PTagQCD.asc");

  int iet, ieta, in1;
  int j;  
  string line;

  while (PTagFile) {
    getline(PTagFile,line,' ');
    stringstream pippo1(line);
    pippo1 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo2(line); 
    pippo2 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo3(line); 
    pippo3 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo4(line); 
    pippo4 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo5(line); 
    pippo5 >> N1HETL[j];
    getline(PTagFile,line,' ');
    stringstream pippo6(line); 
    pippo6 >> N0HETL[j];
    getline(PTagFile,line,' ');
    stringstream pippo7(line); 
    pippo7 >> PHETL[j];
    getline(PTagFile,line);
    stringstream pippo8(line); 
    pippo8 >> PHETLS[j];
    getline(PTagFile,line,' ');
    stringstream pippo11(line);
    pippo11 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo12(line); 
    pippo12 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo13(line); 
    pippo13 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo14(line); 
    pippo14 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo15(line); 
    pippo15 >> N1HETM[j];
    getline(PTagFile,line,' ');
    stringstream pippo16(line); 
    pippo16 >> N0HETM[j];
    getline(PTagFile,line,' ');
    stringstream pippo17(line); 
    pippo17 >> PHETM[j];
    getline(PTagFile,line);
    stringstream pippo18(line); 
    pippo18 >> PHETMS[j];
    getline(PTagFile,line,' ');
    stringstream pippo21(line);
    pippo21 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo22(line); 
    pippo22 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo23(line); 
    pippo23 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo24(line); 
    pippo24 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo25(line); 
    pippo25 >> N1HETT[j];
    getline(PTagFile,line,' ');
    stringstream pippo26(line); 
    pippo26 >> N0HETT[j];
    getline(PTagFile,line,' ');
    stringstream pippo27(line); 
    pippo27 >> PHETT[j];
    getline(PTagFile,line);
    stringstream pippo28(line); 
    pippo28 >> PHETTS[j];
    getline(PTagFile,line,' ');
    stringstream pippo31(line);
    pippo31 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo32(line); 
    pippo32 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo33(line); 
    pippo33 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo34(line); 
    pippo34 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo35(line); 
    pippo35 >> N1HPTL[j];
    getline(PTagFile,line,' ');
    stringstream pippo36(line); 
    pippo36 >> N0HPTL[j];
    getline(PTagFile,line,' ');
    stringstream pippo37(line); 
    pippo37 >> PHPTL[j];
    getline(PTagFile,line);
    stringstream pippo38(line); 
    pippo38 >> PHPTLS[j];
    getline(PTagFile,line,' ');
    stringstream pippo41(line);
    pippo41 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo42(line); 
    pippo42 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo43(line); 
    pippo43 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo44(line); 
    pippo44 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo45(line); 
    pippo45 >> N1HPTM[j];
    getline(PTagFile,line,' ');
    stringstream pippo46(line); 
    pippo46 >> N0HPTM[j];
    getline(PTagFile,line,' ');
    stringstream pippo47(line); 
    pippo47 >> PHPTM[j];
    getline(PTagFile,line);
    stringstream pippo48(line); 
    pippo48 >> PHPTMS[j];
    getline(PTagFile,line,' ');
    stringstream pippo51(line);
    pippo51 >> iet;
    getline(PTagFile,line,' ');
    stringstream pippo52(line); 
    pippo52 >> ieta;
    getline(PTagFile,line,' ');
    stringstream pippo53(line); 
    pippo53 >> in1;
    getline(PTagFile,line,' ');
    stringstream pippo54(line); 
    pippo54 >> j;
    getline(PTagFile,line,' ');
    stringstream pippo55(line); 
    pippo55 >> N1HPTT[j];
    getline(PTagFile,line,' ');
    stringstream pippo56(line); 
    pippo56 >> N0HPTT[j];
    getline(PTagFile,line,' ');
    stringstream pippo57(line); 
    pippo57 >> PHPTT[j];
    getline(PTagFile,line);
    stringstream pippo58(line); 
    pippo58 >> PHPTTS[j];
  }

  PTagFile.close();

}


TDAna::~TDAna()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  // Draw the histograms
  // -------------------

  OutputFile->Write();

}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void TDAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  double NTagsM = 0;      // Number of medium b-tags 
  double NTagsT = 0;      // Number of tight b-tags

  double UncorrSumEx=0;
  double UncorrSumEy=0;
  double UncorrSumEt=0;  // SumEt with all jets before corrections
  double UncorrHt=0;     // Ht with all jets before corrections
  double UncorrHt0=0;     // Ht with all jets before corrections

  double CorrSumEx=0;
  double CorrSumEy=0;
  double CorrSumEt=0;    // SumEt with all jets after corrections 
  double CorrHt=0;       // Ht with all jets after corrections
  double CorrHt0=0;       // Ht with all jets after corrections

  double GoodSumEt=0;    // SumEt with all selected jets
  double GoodHt=0;       // Ht with all selected jets
  double GoodHt2=0;      // Ht with all selected jets and standard MET
  double GoodHt0=0;       // Ht with all selected jets
  double GoodHt20=0;      // Ht with all selected jets and standard MET

  double dpmin=20;  // DP min MET / selected jets
  double dp1st=20;  // DP MET/leading jet
  double dp2nd=20;  // DP MET/2nd leading jet
  double dp3rd=20;  // DP MET/3rd leading jet

  double m3best=-999; // Massa del tripletto piu' vicino a 172 GeV
  double mwbest=-999;    // Massa del doppietto piu' vicino a 80.4 GeV tra i tre nel tripletto definito sopra
  double chi2=1000;   // Chiquadro costruito con le due masse qui sopra
  double m45best=-999;
  double chi2ext=1000;
  double m45bestall=-999;
  double chi2extall=1000;
  double dp12 = 0;        // Delta phi jet 1 jet 2
  double dpbb = 0;    // delta phi bb
  double dpbball = 0;    // delta phi bb
  double m_others = -999; // mass of jets not part of best triplet and pair
  double mbbnohmax = -999; // mass of b pair with highest mass excluding b-jets from h decay
  double dpbbnohmax = 0; // angle of b pair with highest mass excluding b-jets from h decay

  // Offline cuts
  // ------------
  bool NJetsCut = false;   // Whether event passes NJet cut
  bool MEtSigCut = false;  // Whether event pass MET significance cut

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
  int    JN1[100];
  int    JBin[100];
  bool   JHEM[100];
  double NHSJ = 8; // number of jets from hard subprocess  <---- VERY IMPORTANT PARAMETER
  double NHEM = 0; // number of medium tags

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

      UncorrSumEx+=cal->uncorrEt()*cos(cal->phi());
      UncorrSumEy+=cal->uncorrEt()*sin(cal->phi());
      UncorrSumEt+=cal->uncorrEt();
      UncorrHt0+=cal->uncorrEt();
      CorrSumEx+=cal->et()*cos(cal->phi());
      CorrSumEy+=cal->et()*sin(cal->phi());
      CorrSumEt+=cal->et();
      CorrHt0+=cal->et();

      if ( cal->et() >= 30. && fabs( cal->eta() ) < 2.5 ) {
        goodIc5Jets=goodIc5Jets+1;
	if ( cal->discriminatorHighEff() > loose_ ) NTagsL++;
	if ( cal->discriminatorHighEff() > medium_ ) NTagsM++;
	if ( cal->discriminatorHighEff() > tight_ ) NTagsT++;
	GoodHt0+=cal->et();
	GoodHt20+=cal->et();
	GoodSumEt+=cal->et();
	double dp = 3.1415926-fabs(fabs(cal->phi()-caloMET->phi())-3.1415926);
	if ( dp<dpmin ) dpmin=dp;
	goodIc5JetVec.push_back(*cal);

	// Store information on good jets in arrays
	// ----------------------------------------
	if ( iJ<100 ) {
	  JEt[iJ]=cal->et();
	  Jeta[iJ]=cal->eta();
	  Jphi[iJ]=cal->phi();
	  JHET[iJ]=cal->discriminatorHighEff();
	  JHPT[iJ]=cal->discriminatorHighPur();
	  JN1[iJ]=cal->tkNumS1();
	  JHEM[iJ]=false;
	  if ( JHET[iJ]>medium_ ) {
	    JHEM[iJ]=true;
	    if ( iJ<NHSJ ) NHEM++;   // NNNBBB Tags are counted only if in first 8 jets!
	  }
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
	  int a6 = JN1[k];
	  bool a7 = JHEM[k];
	  JEt[k]=JEt[k-1];
	  Jeta[k]=Jeta[k-1];
	  Jphi[k]=Jphi[k-1];
	  JHET[k]=JHET[k-1];
	  JHPT[k]=JHPT[k-1];
	  JN1[k]=JN1[k-1];
	  JHEM[k]=JHEM[k-1];
	  JEt[k-1]=a1;
	  Jeta[k-1]=a2;
	  Jphi[k-1]=a3;
	  JHET[k-1]=a4;
	  JHPT[k-1]=a5;
	  JN1[k-1]=a6;
	  JHEM[k-1]=a7;
	}
      }
    }
    if ( iJ>1 ) dp12 = 3.1415926-fabs(fabs(Jphi[0]-Jphi[1])-3.1415926);
    if ( iJ>0 ) dp1st= 3.1415926-fabs(fabs(Jphi[0]-caloMET->phi())-3.1415926);
    if ( iJ>1 ) dp2nd= 3.1415926-fabs(fabs(Jphi[1]-caloMET->phi())-3.1415926);
    if ( iJ>2 ) dp3rd= 3.1415926-fabs(fabs(Jphi[2]-caloMET->phi())-3.1415926);

    if ( QCD_ ) { // This is for tag parametrization /////////////////////////////////////

      // Assign tag probability to each jet
      // by finding which bin it belongs to
      // ----------------------------------
      double etbins[9]={30.,50.,70.,90.,110.,150.,200.,300.,10000.};
      double etabins[5]={0.,1.,1.5,2.0,2.5};
      double ns1bins[9]={2.,3.,4.,5.,6.,7.,8.,9.,100.};
      int etbin=-1;
      int etabin=-1;
      int ns1bin=-1;

      for ( int i=0; i<iJ; i++ ) {
	// Find out which bin this jet belongs to
	// --------------------------------------
	for ( int iet=0; iet<8; iet++ ) {
	  if ( JEt[i]>etbins[iet] && JEt[i]<=etbins[iet+1] ) {
	    etbin=iet;
	  }
	}
	for ( int ieta=0; ieta<4; ieta++ ) {
	  if ( fabs(Jeta[i])>=etabins[ieta] && fabs(Jeta[i])<etabins[ieta+1] ) {
	    etabin=ieta;
	  }
	}
	for ( int int1=0; int1<8; int1++ ) {
	  if ( JN1[i]>=ns1bins[int1] && JN1[i]<ns1bins[int1+1] ) {
	    ns1bin=int1;
	  }
	}
	if ( etbin>-1 && etabin>-1 && ns1bin>-1 ) {
	  JBin[i] = etbin*100+etabin*10+ns1bin;
	} else {
	  JBin[i] = 999;
	}
      }
      

      // Assign tags to jets in QCD events according to their probability
      // ----------------------------------------------------------------
      int iJmax=NHSJ;  // Search for tags in the first NHSJ jets only!

      if ( iJ<iJmax ) iJmax=iJ;
      for ( int icomb=0; icomb<(int)pow (2.,iJmax); icomb++ ) { // Combinatorial loop ///////////////
	
	// Reset all variables
	// -------------------
	m3best=-999; // Massa del tripletto piu' vicino a 172 GeV
	mwbest=-999;    // Massa del doppietto piu' vicino a 80.4 GeV tra i tre nel tripletto definito sopra
	chi2=1000;   // Chiquadro costruito con le due masse qui sopra
	m45best=-999; // mass of two b-tags closest to 120 GeV not in triplet
	chi2ext=1000; // chi2 computed by adding contribution from distance to 120 GeV
	dpbb = 0;    // delta phi bb of pair above
	m45bestall=-999;  // Same as above, without bothering to exclude jets in triplet
	chi2extall=1000;  // Same as above
	dpbball = 0;    // delta phi bb as above
	m_others = -999; // mass of jets not part of best triplet and pair
	mbbnohmax = -999; // mass of b pair with highest mass excluding b-jets from h decay
	dpbbnohmax = 0; // angle of b pair with highest mass excluding b-jets from h decay
	
	m8=0;
	set=0;
	c8=0;

	// Ok now start producing probability combinations
	// -----------------------------------------------
	double n=(double)icomb;
	double PTOT=1.;
	double PTOTE2=0.;
	NHEM=0;

	for ( int j=iJmax-1; j>=0; j-- ) {
	  if ( n>=pow(2.,j) ) {
	    JHEM[j]=true;
	    n-=pow(2.,j);
	    PTOT=PTOT*PHETM[JBin[j]];
	    PTOTE2+=pow(PHETMS[JBin[j]]/PHETM[JBin[j]],2);
	    NHEM++;
	  } else {
	    JHEM[j]=false;
	    PTOT=PTOT*(1.-PHETM[JBin[j]]);
	    PTOTE2+=pow(PHETMS[JBin[j]]/(1.-PHETM[JBin[j]]),2);
	  }
	}
	PTOTE2=PTOTE2*pow(PTOT,2);

	if ( NHEM>=2 ) {

	  // Compute best top and W masses and determine simple chi2
	  // -------------------------------------------------------
	  double px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,e1,e2,e3;
	  double px4,px5,py4,py5,pz4,pz5,e4,e5;
	  double px,py,pz;
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
	  double nb;
	  for ( int i=0; i<iJ-2 && i<NHSJ-2; i++ ) {
	    for ( int j=i+1; j<iJ-1 && j<NHSJ-1; j++ ) {
	      for ( int k=j+1; k<iJ && k<NHSJ; k++ ) { // limit to 8 the number of jets in hard subprocess
		
		// Require that the triplet owns at least one b-tag
		// ------------------------------------------------
		nb=0.;
		
		if ( JHEM[i] ) nb++; 
		if ( JHEM[j] ) nb++;
		if ( JHEM[k] ) nb++;
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
		    if ( JHEM[i] || JHEM[j] ) bestWhasb=1;
		    if ( !JHEM[k] ) bestbnob=1;
		  }
		  if ( fabs(m13-80.4)<fabs(m12-80.4) && fabs(m13-80.4)<fabs(m23-80.4) ) { 
		    mw=m13;
		    if ( JHEM[i] || JHEM[k] ) bestWhasb=1;
		    if ( !JHEM[j] ) bestbnob=1;
		  }
		  if ( fabs(m23-80.4)<fabs(m13-80.4) && fabs(m23-80.4)<fabs(m12-80.4) ) {
		    mw=m23;
		    if ( JHEM[j] || JHEM[k] ) bestWhasb=1;
		    if ( !JHEM[i] ) bestbnob=1;
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
		    if ( JHEM[ii] && JHEM[jj] ) {
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
	    // Now redo the same, on all b-jets regardless of best triplet
	    // -----------------------------------------------------------
	    for ( int ii=0; ii<iJ-1 && ii<NHSJ-1; ii++ ) {
	      for ( int jj=ii+1; jj<iJ && jj<NHSJ; jj++ ) { // limit to NHSJ the number of jets from h.s.
		// Demand that both b-jets from H be tagged
		// ----------------------------------------
		if ( JHEM[ii] && JHEM[jj] ) {
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
		    m45bestall=m45;
		    double addchi=0;
		    if ( ii==it1 || ii==it2 || ii==it3 || jj==it1 || jj==it2 || jj==it3 ) addchi=3; 
		    chi2extall=sqrt(chi2*chi2+(m45bestall-120)*(m45bestall-120)/200+addchi*addchi);
		    dpbball=3.145926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
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
	      for ( int kk=0; kk<iJ && kk<NHSJ ; kk++ ) { // limit to NHSJ the number of 
		                                          // jets from hard subprocess
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
		if ( k1!=ih1 && k1!=ih2 && JHEM[k1] ) {
		  for ( int k2=k1+1; k2<iJ && k2<NHSJ; k2++ ) { // limit to 8 the number of jets from h.s.
		    if ( k2!=ih1 && k2!=ih2 && JHEM[k2] ) {
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
	  } // end if it1+it2+it3>0      
	  
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
	
	  ////////////////////////////////////////////////////////////////
	  // Now fill histograms
	  // -------------------
	  // Offline cuts
	  // ------------
	  if ( iJ >= 5 ) NJetsCut = true;
	  if ( caloMET->mEtSig()>3 ) MEtSigCut = true;
	  
	  // Fill kinematical histograms
	  // ---------------------------
	  double UncorrMEt =sqrt(UncorrSumEx*UncorrSumEx+UncorrSumEy*UncorrSumEy);
	  double CorrMEt = sqrt(CorrSumEx*CorrSumEx+CorrSumEy*CorrSumEy);
	  UncorrHt=UncorrHt0+UncorrMEt;
	  CorrHt=CorrHt0+CorrMEt;
	  GoodHt=GoodHt0+CorrMEt;
	  GoodHt2=GoodHt20+caloMET->et();
	  
	  double UncorrMEtSig=0;
	  if ( UncorrSumEt>0 ) UncorrMEtSig = UncorrMEt/sqrt(UncorrSumEt);
	  double CorrMEtSig = 0;
	  if ( CorrSumEt>0 ) CorrMEtSig = CorrMEt/sqrt(CorrSumEt);  
	  
	  // Missing Et significance computed with tuned resolution function
	  // ---------------------------------------------------------------
	  double metsignew=0;
	  if ( caloMET->sumEt()>0 ) metsignew = caloMET->et()/(sqrt(2)*0.8033*pow(caloMET->sumEt(),0.5004));
	  
	  if ( PTOT>0 ) {

	    // Fill regular histograms
	    // -----------------------
	    UncorrMEt_SumEt_->Fill(caloMET->sumEt(),UncorrMEt,PTOT);
	    CorrMEt_SumEt_->Fill(caloMET->sumEt(),CorrMEt,PTOT);
	    MEt_SumEt_->Fill(caloMET->sumEt(),caloMET->et(),PTOT);
	    UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt,PTOT);
	    CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt,PTOT);
	    MEt_SumEtC_->Fill(CorrSumEt,caloMET->et(),PTOT);
	    UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt,PTOT);
	    CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt,PTOT);
	    MEt_SumEtJ_->Fill(GoodSumEt,caloMET->et(),PTOT);
	    
	    // Apply trigger requirement
	    // -------------------------
	    if ( response && NJetsCut && NHEM>=2 ) {
	      NJets_->Fill(goodIc5Jets,PTOT);
	      UncorrHt_->Fill(UncorrHt,PTOT);
	      CorrHt_->Fill(CorrHt,PTOT);
	      GoodHt_->Fill(GoodHt,PTOT);
	      GoodHt2_->Fill(GoodHt2,PTOT);
	      UncorrSumEt_->Fill(UncorrSumEt,PTOT);
	      CorrSumEt_->Fill(CorrSumEt,PTOT);
	      GoodSumEt_->Fill(GoodSumEt,PTOT);
	      MEt_->Fill(caloMET->et(),PTOT);
	      MEtSig_->Fill(caloMET->mEtSig(),PTOT);
	      MEtSigNew_->Fill(metsignew,PTOT);
	      MEtDPM_->Fill(dpmin,PTOT);
	      MEtDP1_->Fill(dp1st,PTOT);
	      MEtDP2_->Fill(dp2nd,PTOT);
	      MEtDP3_->Fill(dp3rd,PTOT);
	      UncorrMEtSig_->Fill(UncorrMEtSig,PTOT);
	      CorrMEtSig_->Fill(CorrMEtSig,PTOT);
	      M3best_->Fill(m3best,PTOT);
	      Mwbest_->Fill(mwbest,PTOT);
	      Chi2mass_->Fill(chi2,PTOT);
	      M45best_->Fill(m45best,PTOT);
	      Chi2ext_->Fill(chi2ext,PTOT);
	      MEx_SumEt_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOT);
	      MEx_SumEt_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOT);
	      DP12_->Fill(dp12,PTOT);
	      DPbb_->Fill(dpbb,PTOT);
	      M_others_->Fill(m_others,PTOT);
	      Mbbnoh_->Fill(mbbnohmax,PTOT);
	      DPbbnoh_->Fill(dpbbnohmax,PTOT);
	      M8_->Fill(m8,PTOT);
	      C8_->Fill(c8,PTOT);
	      M45bestall_->Fill(m45bestall,PTOT);
	      Chi2extall_->Fill(chi2extall,PTOT);
	      DPbball_->Fill(dpbball,PTOT);
	      
	      if ( MEtSigCut ) {
		NJetsS_->Fill(goodIc5Jets,PTOT);
		UncorrHtS_->Fill(UncorrHt,PTOT);
		CorrHtS_->Fill(CorrHt,PTOT);
		GoodHtS_->Fill(GoodHt,PTOT);
		GoodHt2S_->Fill(GoodHt2,PTOT);
		UncorrSumEtS_->Fill(UncorrSumEt,PTOT);
		CorrSumEtS_->Fill(CorrSumEt,PTOT);
		GoodSumEtS_->Fill(GoodSumEt,PTOT);
		MEtS_->Fill(caloMET->et(),PTOT);
		MEtSigS_->Fill(caloMET->mEtSig(),PTOT);
		MEtSigNewS_->Fill(metsignew,PTOT);
		MEtDPMS_->Fill(dpmin,PTOT);
		MEtDP1S_->Fill(dp1st,PTOT);
		MEtDP2S_->Fill(dp2nd,PTOT);
		MEtDP3S_->Fill(dp3rd,PTOT);
		UncorrMEtSigS_->Fill(UncorrMEtSig,PTOT);
		CorrMEtSigS_->Fill(CorrMEtSig,PTOT);
		M3bestS_->Fill(m3best,PTOT);
		MwbestS_->Fill(mwbest,PTOT);
		Chi2massS_->Fill(chi2,PTOT);
		M45bestS_->Fill(m45best,PTOT);
		Chi2extS_->Fill(chi2ext,PTOT);
		MEx_SumEtS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOT);
		MEx_SumEtS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOT);
		DP12S_->Fill(dp12,PTOT);
		DPbbS_->Fill(dpbb,PTOT);
		M_othersS_->Fill(m_others,PTOT);
		MbbnohS_->Fill(mbbnohmax,PTOT);
		DPbbnohS_->Fill(dpbbnohmax,PTOT);
		M8S_->Fill(m8,PTOT);
		C8S_->Fill(c8,PTOT);
		M45bestallS_->Fill(m45bestall,PTOT);
		Chi2extallS_->Fill(chi2extall,PTOT);
		DPbballS_->Fill(dpbball,PTOT);
	      
		if ( NHEM>=3 ) {
		  NJetsSS_->Fill(goodIc5Jets,PTOT);
		  UncorrHtSS_->Fill(UncorrHt,PTOT);
		  CorrHtSS_->Fill(CorrHt,PTOT);
		  GoodHtSS_->Fill(GoodHt,PTOT);
		  GoodHt2SS_->Fill(GoodHt2,PTOT);
		  UncorrSumEtSS_->Fill(UncorrSumEt,PTOT);
		  CorrSumEtSS_->Fill(CorrSumEt,PTOT);
		  GoodSumEtSS_->Fill(GoodSumEt,PTOT);
		  MEtSS_->Fill(caloMET->et(),PTOT);
		  MEtSigSS_->Fill(caloMET->mEtSig(),PTOT);
		  MEtSigNewSS_->Fill(metsignew,PTOT);
		  MEtDPMSS_->Fill(dpmin,PTOT);
		  MEtDP1SS_->Fill(dp1st,PTOT);
		  MEtDP2SS_->Fill(dp2nd,PTOT);
		  MEtDP3SS_->Fill(dp3rd,PTOT);
		  UncorrMEtSigSS_->Fill(UncorrMEtSig,PTOT);
		  CorrMEtSigSS_->Fill(CorrMEtSig,PTOT);
		  M3bestSS_->Fill(m3best,PTOT);
		  MwbestSS_->Fill(mwbest,PTOT);
		  Chi2massSS_->Fill(chi2,PTOT);
		  M45bestSS_->Fill(m45best,PTOT);
		  Chi2extSS_->Fill(chi2ext,PTOT);
		  MEx_SumEtSS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOT);
		  MEx_SumEtSS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOT);
		  DP12SS_->Fill(dp12,PTOT);
		  DPbbSS_->Fill(dpbb,PTOT);
		  M_othersSS_->Fill(m_others,PTOT);      
		  MbbnohSS_->Fill(mbbnohmax,PTOT);
		  DPbbnohSS_->Fill(dpbbnohmax,PTOT);
		  M8SS_->Fill(m8,PTOT);
		  C8SS_->Fill(c8,PTOT);
		  M45bestallSS_->Fill(m45bestall,PTOT);
		  Chi2extallSS_->Fill(chi2extall,PTOT);
		  DPbballSS_->Fill(dpbball,PTOT);
		}
		
		if ( NHEM>=4 ) {
		  NJetsSSS_->Fill(goodIc5Jets,PTOT);
		  UncorrHtSSS_->Fill(UncorrHt,PTOT);
		  CorrHtSSS_->Fill(CorrHt,PTOT);
		  GoodHtSSS_->Fill(GoodHt,PTOT);
		  GoodHt2SSS_->Fill(GoodHt2,PTOT);
		  UncorrSumEtSSS_->Fill(UncorrSumEt,PTOT);
		  CorrSumEtSSS_->Fill(CorrSumEt,PTOT);
		  GoodSumEtSSS_->Fill(GoodSumEt,PTOT);
		  MEtSSS_->Fill(caloMET->et(),PTOT);
		  MEtSigSSS_->Fill(caloMET->mEtSig(),PTOT);
		  MEtSigNewSSS_->Fill(metsignew,PTOT);
		  MEtDPMSSS_->Fill(dpmin,PTOT);
		  MEtDP1SSS_->Fill(dp1st,PTOT);
		  MEtDP2SSS_->Fill(dp2nd,PTOT);
		  MEtDP3SSS_->Fill(dp3rd,PTOT);
		  UncorrMEtSigSSS_->Fill(UncorrMEtSig,PTOT);
		  CorrMEtSigSSS_->Fill(CorrMEtSig,PTOT);
		  M3bestSSS_->Fill(m3best,PTOT);
		  MwbestSSS_->Fill(mwbest,PTOT);
		  Chi2massSSS_->Fill(chi2,PTOT);
		  M45bestSSS_->Fill(m45best,PTOT);
		  Chi2extSSS_->Fill(chi2ext,PTOT);
		  MEx_SumEtSSS_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOT);
		  MEx_SumEtSSS_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOT);
		  DP12SSS_->Fill(dp12,PTOT);
		  DPbbSSS_->Fill(dpbb,PTOT);
		  M_othersSSS_->Fill(m_others,PTOT);      
		  MbbnohSSS_->Fill(mbbnohmax,PTOT);
		  DPbbnohSSS_->Fill(dpbbnohmax,PTOT);
		  M8SSS_->Fill(m8,PTOT);
		  C8SSS_->Fill(c8,PTOT);
		  M45bestallSSS_->Fill(m45bestall,PTOT);
		  Chi2extallSSS_->Fill(chi2extall,PTOT);
		  DPbballSSS_->Fill(dpbball,PTOT);
		}

	      }
	    } 

	    // Fill error histograms now
	    // -------------------------

	    UncorrMEt_SumEt_->Fill(caloMET->sumEt(),UncorrMEt,PTOTE2);
	    CorrMEt_SumEt_->Fill(caloMET->sumEt(),CorrMEt,PTOTE2);
	    MEt_SumEt_->Fill(caloMET->sumEt(),caloMET->et(),PTOTE2);
	    UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt,PTOTE2);
	    CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt,PTOTE2);
	    MEt_SumEtC_->Fill(CorrSumEt,caloMET->et(),PTOTE2);
	    UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt,PTOTE2);
	    CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt,PTOTE2);
	    MEt_SumEtJ_->Fill(GoodSumEt,caloMET->et(),PTOTE2);
	    
	    // Apply trigger requirement
	    // -------------------------
	    if ( response && NJetsCut && NHEM>=2 ) {
	      NJetsW_->Fill(goodIc5Jets,PTOTE2);
	      UncorrHtW_->Fill(UncorrHt,PTOTE2);
	      CorrHtW_->Fill(CorrHt,PTOTE2);
	      GoodHtW_->Fill(GoodHt,PTOTE2);
	      GoodHt2W_->Fill(GoodHt2,PTOTE2);
	      UncorrSumEtW_->Fill(UncorrSumEt,PTOTE2);
	      CorrSumEtW_->Fill(CorrSumEt,PTOTE2);
	      GoodSumEtW_->Fill(GoodSumEt,PTOTE2);
	      MEtW_->Fill(caloMET->et(),PTOTE2);
	      MEtSigW_->Fill(caloMET->mEtSig(),PTOTE2);
	      MEtSigNewW_->Fill(metsignew,PTOTE2);
	      MEtDPMW_->Fill(dpmin,PTOTE2);
	      MEtDP1W_->Fill(dp1st,PTOTE2);
	      MEtDP2W_->Fill(dp2nd,PTOTE2);
	      MEtDP3W_->Fill(dp3rd,PTOTE2);
	      UncorrMEtSigW_->Fill(UncorrMEtSig,PTOTE2);
	      CorrMEtSigW_->Fill(CorrMEtSig,PTOTE2);
	      M3bestW_->Fill(m3best,PTOTE2);
	      MwbestW_->Fill(mwbest,PTOTE2);
	      Chi2massW_->Fill(chi2,PTOTE2);
	      M45bestW_->Fill(m45best,PTOTE2);
	      Chi2extW_->Fill(chi2ext,PTOTE2);
	      MEx_SumEtW_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOTE2);
	      MEx_SumEtW_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOTE2);
	      DP12W_->Fill(dp12,PTOTE2);
	      DPbbW_->Fill(dpbb,PTOTE2);
	      M_othersW_->Fill(m_others,PTOTE2);
	      MbbnohW_->Fill(mbbnohmax,PTOTE2);
	      DPbbnohW_->Fill(dpbbnohmax,PTOTE2);
	      M8W_->Fill(m8,PTOTE2);
	      C8W_->Fill(c8,PTOTE2);
	      M45bestallW_->Fill(m45bestall,PTOTE2);
	      Chi2extallW_->Fill(chi2extall,PTOTE2);
	      DPbballW_->Fill(dpbball,PTOTE2);
	      
	      if ( MEtSigCut ) {
		NJetsSW_->Fill(goodIc5Jets,PTOTE2);
		UncorrHtSW_->Fill(UncorrHt,PTOTE2);
		CorrHtSW_->Fill(CorrHt,PTOTE2);
		GoodHtSW_->Fill(GoodHt,PTOTE2);
		GoodHt2SW_->Fill(GoodHt2,PTOTE2);
		UncorrSumEtSW_->Fill(UncorrSumEt,PTOTE2);
		CorrSumEtSW_->Fill(CorrSumEt,PTOTE2);
		GoodSumEtSW_->Fill(GoodSumEt,PTOTE2);
		MEtSW_->Fill(caloMET->et(),PTOTE2);
		MEtSigSW_->Fill(caloMET->mEtSig(),PTOTE2);
		MEtSigNewSW_->Fill(metsignew,PTOTE2);
		MEtDPMSW_->Fill(dpmin,PTOTE2);
		MEtDP1SW_->Fill(dp1st,PTOTE2);
		MEtDP2SW_->Fill(dp2nd,PTOTE2);
		MEtDP3SW_->Fill(dp3rd,PTOTE2);
		UncorrMEtSigSW_->Fill(UncorrMEtSig,PTOTE2);
		CorrMEtSigSW_->Fill(CorrMEtSig,PTOTE2);
		M3bestSW_->Fill(m3best,PTOTE2);
		MwbestSW_->Fill(mwbest,PTOTE2);
		Chi2massSW_->Fill(chi2,PTOTE2);
		M45bestSW_->Fill(m45best,PTOTE2);
		Chi2extSW_->Fill(chi2ext,PTOTE2);
		MEx_SumEtSW_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOTE2);
		MEx_SumEtSW_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOTE2);
		DP12SW_->Fill(dp12,PTOTE2);
		DPbbSW_->Fill(dpbb,PTOTE2);
		M_othersSW_->Fill(m_others,PTOTE2);
		MbbnohSW_->Fill(mbbnohmax,PTOTE2);
		DPbbnohSW_->Fill(dpbbnohmax,PTOTE2);
		M8SW_->Fill(m8,PTOTE2);
		C8SW_->Fill(c8,PTOTE2);
		M45bestallSW_->Fill(m45bestall,PTOTE2);
		Chi2extallSW_->Fill(chi2extall,PTOTE2);
		DPbballSW_->Fill(dpbball,PTOTE2);
	      
		if ( NHEM>=3 ) {
		  NJetsSSW_->Fill(goodIc5Jets,PTOTE2);
		  UncorrHtSSW_->Fill(UncorrHt,PTOTE2);
		  CorrHtSSW_->Fill(CorrHt,PTOTE2);
		  GoodHtSSW_->Fill(GoodHt,PTOTE2);
		  GoodHt2SSW_->Fill(GoodHt2,PTOTE2);
		  UncorrSumEtSSW_->Fill(UncorrSumEt,PTOTE2);
		  CorrSumEtSSW_->Fill(CorrSumEt,PTOTE2);
		  GoodSumEtSSW_->Fill(GoodSumEt,PTOTE2);
		  MEtSSW_->Fill(caloMET->et(),PTOTE2);
		  MEtSigSSW_->Fill(caloMET->mEtSig(),PTOTE2);
		  MEtSigNewSSW_->Fill(metsignew,PTOTE2);
		  MEtDPMSSW_->Fill(dpmin,PTOTE2);
		  MEtDP1SSW_->Fill(dp1st,PTOTE2);
		  MEtDP2SSW_->Fill(dp2nd,PTOTE2);
		  MEtDP3SSW_->Fill(dp3rd,PTOTE2);
		  UncorrMEtSigSSW_->Fill(UncorrMEtSig,PTOTE2);
		  CorrMEtSigSSW_->Fill(CorrMEtSig,PTOTE2);
		  M3bestSSW_->Fill(m3best,PTOTE2);
		  MwbestSSW_->Fill(mwbest,PTOTE2);
		  Chi2massSSW_->Fill(chi2,PTOTE2);
		  M45bestSSW_->Fill(m45best,PTOTE2);
		  Chi2extSSW_->Fill(chi2ext,PTOTE2);
		  MEx_SumEtSSW_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOTE2);
		  MEx_SumEtSSW_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOTE2);
		  DP12SSW_->Fill(dp12,PTOTE2);
		  DPbbSSW_->Fill(dpbb,PTOTE2);
		  M_othersSSW_->Fill(m_others,PTOTE2);      
		  MbbnohSSW_->Fill(mbbnohmax,PTOTE2);
		  DPbbnohSSW_->Fill(dpbbnohmax,PTOTE2);
		  M8SSW_->Fill(m8,PTOTE2);
		  C8SSW_->Fill(c8,PTOTE2);
		  M45bestallSSW_->Fill(m45bestall,PTOTE2);
		  Chi2extallSSW_->Fill(chi2extall,PTOTE2);
		  DPbballSSW_->Fill(dpbball,PTOTE2);
		}
		
		if ( NHEM>=4 ) {
		  NJetsSSSW_->Fill(goodIc5Jets,PTOTE2);
		  UncorrHtSSSW_->Fill(UncorrHt,PTOTE2);
		  CorrHtSSSW_->Fill(CorrHt,PTOTE2);
		  GoodHtSSSW_->Fill(GoodHt,PTOTE2);
		  GoodHt2SSSW_->Fill(GoodHt2,PTOTE2);
		  UncorrSumEtSSSW_->Fill(UncorrSumEt,PTOTE2);
		  CorrSumEtSSSW_->Fill(CorrSumEt,PTOTE2);
		  GoodSumEtSSSW_->Fill(GoodSumEt,PTOTE2);
		  MEtSSSW_->Fill(caloMET->et(),PTOTE2);
		  MEtSigSSSW_->Fill(caloMET->mEtSig(),PTOTE2);
		  MEtSigNewSSSW_->Fill(metsignew,PTOTE2);
		  MEtDPMSSSW_->Fill(dpmin,PTOTE2);
		  MEtDP1SSSW_->Fill(dp1st,PTOTE2);
		  MEtDP2SSSW_->Fill(dp2nd,PTOTE2);
		  MEtDP3SSSW_->Fill(dp3rd,PTOTE2);
		  UncorrMEtSigSSSW_->Fill(UncorrMEtSig,PTOTE2);
		  CorrMEtSigSSSW_->Fill(CorrMEtSig,PTOTE2);
		  M3bestSSSW_->Fill(m3best,PTOTE2);
		  MwbestSSSW_->Fill(mwbest,PTOTE2);
		  Chi2massSSSW_->Fill(chi2,PTOTE2);
		  M45bestSSSW_->Fill(m45best,PTOTE2);
		  Chi2extSSSW_->Fill(chi2ext,PTOTE2);
		  MEx_SumEtSSSW_->Fill(caloMET->sumEt(),caloMET->et()*cos(caloMET->phi()),PTOTE2);
		  MEx_SumEtSSSW_->Fill(caloMET->sumEt(),caloMET->et()*sin(caloMET->phi()),PTOTE2);
		  DP12SSSW_->Fill(dp12,PTOTE2);
		  DPbbSSSW_->Fill(dpbb,PTOTE2);
		  M_othersSSSW_->Fill(m_others,PTOTE2);      
		  MbbnohSSSW_->Fill(mbbnohmax,PTOTE2);
		  DPbbnohSSSW_->Fill(dpbbnohmax,PTOTE2);
		  M8SSSW_->Fill(m8,PTOTE2);
		  C8SSSW_->Fill(c8,PTOTE2);
		  M45bestallSSSW_->Fill(m45bestall,PTOTE2);
		  Chi2extallSSSW_->Fill(chi2extall,PTOTE2);
		  DPbballSSSW_->Fill(dpbball,PTOTE2);
		}
	      }
	    }

	  } // end if PTOT>0
	}  // end if ntags>=2

      } // end for tag combinatorial
      
    } else { // if not QCD /////////////////////////////////////////////

      // Compute best top and W masses and determine simple chi2
      // -------------------------------------------------------
      double px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,e1,e2,e3;
      double px4,px5,py4,py5,pz4,pz5,e4,e5;
      double px,py,pz;
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

	    // Require that the triplet owns at least one b-tag
	    // ------------------------------------------------
	    double nb=0.;

	    if ( JHEM[i] ) nb++; 
	    if ( JHEM[j] ) nb++;
	    if ( JHEM[k] ) nb++;
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
		if ( JHEM[i] || JHEM[j] ) bestWhasb=1;
		if ( !JHEM[k] ) bestbnob=1;
	      }
	      if ( fabs(m13-80.4)<fabs(m12-80.4) && fabs(m13-80.4)<fabs(m23-80.4) ) { 
		mw=m13;
		if ( JHEM[i] || JHEM[k] ) bestWhasb=1;
		if ( !JHEM[j] ) bestbnob=1;
	      }
	      if ( fabs(m23-80.4)<fabs(m13-80.4) && fabs(m23-80.4)<fabs(m12-80.4) ) {
		mw=m23;
		if ( JHEM[j] || JHEM[k] ) bestWhasb=1;
		if ( !JHEM[i] ) bestbnob=1;
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
		if ( JHEM[ii] && JHEM[jj] ) {
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
	// Now redo the same, on all b-jets regardless of best triplet
	// -----------------------------------------------------------
	for ( int ii=0; ii<iJ-1 && ii<NHSJ-1; ii++ ) {
	  for ( int jj=ii+1; jj<iJ && jj<NHSJ; jj++ ) { // limit to NHSJ the number of jets from h.s.
	    // Demand that both b-jets from H be tagged
	    // ----------------------------------------
	    if ( JHEM[ii] && JHEM[jj] ) {
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
		m45bestall=m45;
		double addchi=0;
		if ( ii==it1 || ii==it2 || ii==it3 || jj==it1 || jj==it2 || jj==it3 ) addchi=3; 
		chi2extall=sqrt(chi2*chi2+(m45bestall-120)*(m45bestall-120)/200+addchi*addchi);
		dpbball=3.145926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
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
	    if ( k1!=ih1 && k1!=ih2 && JHEM[k1] ) {
	      for ( int k2=k1+1; k2<iJ && k2<NHSJ; k2++ ) { // limit to 8 the number of jets from h.s.
		if ( k2!=ih1 && k2!=ih2 && JHEM[k2] ) {
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
      
      ////////////////////////////////////////////////////////////////
      // Now fill histograms
      // -------------------
      // Offline cuts
      // ------------
      if ( iJ >= 5 ) NJetsCut = true;
      if ( caloMET->mEtSig()>3 ) MEtSigCut = true;
      
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
      if ( response && NJetsCut && NHEM>=2 ) {
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
	M45bestall_->Fill(m45bestall);
	Chi2extall_->Fill(chi2extall);
	DPbball_->Fill(dpbball);
	
	if ( MEtSigCut ) {
	  NJetsS_->Fill(goodIc5Jets);
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
	  M45bestallS_->Fill(m45bestall);
	  Chi2extallS_->Fill(chi2extall);
	  DPbballS_->Fill(dpbball);
	  
	  if ( NHEM>=3 ) {
	    NJetsSS_->Fill(goodIc5Jets);
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
	    M45bestallSS_->Fill(m45bestall);
	    Chi2extallSS_->Fill(chi2extall);
	    DPbballSS_->Fill(dpbball);
	  }
	  
	  if ( NHEM>=4 ) {
	    NJetsSSS_->Fill(goodIc5Jets);
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
	    M45bestallSSS_->Fill(m45bestall);
	    Chi2extallSSS_->Fill(chi2extall);
	    DPbballSSS_->Fill(dpbball);

	    // Study number of events with 4HEM tags vs N jets on which to look for them
	    // -------------------------------------------------------------------------
	    if ( NHEM>=4 ) {
	      for ( int i=0; i<iJ; i++ ) {
		int nt=0;
		for ( int j=0; j<=i; j++ ) {
		  if ( JHEM[j] ) nt++;
		}
		N4NJSSS_->Fill((double)i);
		if ( nt>=4 ) E4NJSSS_->Fill((double)i);
	      }
	    }
	  }
	}
      }
    }
  } // end if caloJETS.size()
  else {
    std::cout << "ATTENTION: Jet collection empty" << std::endl;
  }    

}

// ------------ method called once each job just before starting event loop  ------------
// --------------------------------------------------------------------------------------
void TDAna::beginJob(const edm::EventSetup&) {
}


// ------------ method called once each job just after ending the event loop  ------------
// ---------------------------------------------------------------------------------------
void TDAna::endJob() {

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
  M45bestall_->Write();
  Chi2extall_->Write();
  DPbball_->Write();
  
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
  M45bestallS_->Write();
  Chi2extallS_->Write();
  DPbballS_->Write();

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
  M45bestallSS_->Write();
  Chi2extallSS_->Write();
  DPbballSS_->Write();

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
  M45bestallSSS_->Write();
  Chi2extallSSS_->Write();
  DPbballSSS_->Write();
  N4NJSSS_->Write();
  E4NJSSS_->Write();

  // Histograms with errors squared
  // ------------------------------
  NJetsW_->Write();
  UncorrHtW_->Write();
  CorrHtW_->Write();
  GoodHtW_->Write();
  GoodHt2W_->Write();
  UncorrSumEtW_->Write();
  CorrSumEtW_->Write();
  GoodSumEtW_->Write();
  MEtW_->Write();
  MEtSigW_->Write();
  MEtSigNewW_->Write();
  MEtDPMW_->Write();
  MEtDP1W_->Write();
  MEtDP2W_->Write();
  MEtDP3W_->Write();
  UncorrMEtSigW_->Write();
  CorrMEtSigW_->Write();
  M3bestW_->Write();
  MwbestW_->Write();
  Chi2massW_->Write();
  M45bestW_->Write();
  Chi2extW_->Write();
  MEx_SumEtW_->Write();
  DP12W_->Write();
  DPbbW_->Write();
  M_othersW_->Write();
  MbbnohW_->Write();
  DPbbnohW_->Write();
  M8W_->Write();
  C8W_->Write();
  M45bestallW_->Write();
  Chi2extallW_->Write();
  DPbballW_->Write();
  
  // Histograms for events passing trigger and NJet cut
  // --------------------------------------------------
  NJetsSW_->Write();
  UncorrHtSW_->Write();
  CorrHtSW_->Write();
  GoodHtSW_->Write();
  GoodHt2SW_->Write();
  UncorrSumEtSW_->Write();
  CorrSumEtSW_->Write();
  GoodSumEtSW_->Write();
  MEtSW_->Write();
  MEtSigSW_->Write();
  MEtSigNewSW_->Write();
  MEtDPMSW_->Write();
  MEtDP1SW_->Write();
  MEtDP2SW_->Write();
  MEtDP3SW_->Write();
  UncorrMEtSigSW_->Write();
  CorrMEtSigSW_->Write();
  M3bestSW_->Write();
  MwbestSW_->Write();
  Chi2massSW_->Write();
  M45bestSW_->Write();
  Chi2extSW_->Write();
  MEx_SumEtSW_->Write();
  DP12SW_->Write();
  DPbbSW_->Write();
  M_othersSW_->Write();
  MbbnohSW_->Write();
  DPbbnohSW_->Write();
  M8SW_->Write();
  C8SW_->Write();
  M45bestallSW_->Write();
  Chi2extallSW_->Write();
  DPbballSW_->Write();

  // Histograms for events passing trigger, NJet, and METS cut
  // ---------------------------------------------------------
  NJetsSSW_->Write();
  UncorrHtSSW_->Write();
  CorrHtSSW_->Write();
  GoodHtSSW_->Write();
  GoodHt2SSW_->Write();
  UncorrSumEtSSW_->Write();
  CorrSumEtSSW_->Write();
  GoodSumEtSSW_->Write();
  MEtSSW_->Write();
  MEtSigSSW_->Write();
  MEtSigNewSSW_->Write();
  MEtDPMSSW_->Write();
  MEtDP1SSW_->Write();
  MEtDP2SSW_->Write();
  MEtDP3SSW_->Write();
  UncorrMEtSigSSW_->Write();
  CorrMEtSigSSW_->Write();
  M3bestSSW_->Write();
  MwbestSSW_->Write();
  Chi2massSSW_->Write();
  M45bestSSW_->Write();
  Chi2extSSW_->Write();
  MEx_SumEtSSW_->Write();
  DP12SSW_->Write();
  DPbbSSW_->Write();
  M_othersSSW_->Write();
  MbbnohSSW_->Write();
  DPbbnohSSW_->Write();
  M8SSW_->Write();
  C8SSW_->Write();
  M45bestallSSW_->Write();
  Chi2extallSSW_->Write();
  DPbballSSW_->Write();

  // Histograms for events passing trigger, NJet, and METS cut
  // ---------------------------------------------------------
  NJetsSSSW_->Write();
  UncorrHtSSSW_->Write();
  CorrHtSSSW_->Write();
  GoodHtSSSW_->Write();
  GoodHt2SSSW_->Write();
  UncorrSumEtSSSW_->Write();
  CorrSumEtSSSW_->Write();
  GoodSumEtSSSW_->Write();
  MEtSSSW_->Write();
  MEtSigSSSW_->Write();
  MEtSigNewSSSW_->Write();
  MEtDPMSSSW_->Write();
  MEtDP1SSSW_->Write();
  MEtDP2SSSW_->Write();
  MEtDP3SSSW_->Write();
  UncorrMEtSigSSSW_->Write();
  CorrMEtSigSSSW_->Write();
  M3bestSSSW_->Write();
  MwbestSSSW_->Write();
  Chi2massSSSW_->Write();
  M45bestSSSW_->Write();
  Chi2extSSSW_->Write();
  MEx_SumEtSSSW_->Write();
  DP12SSSW_->Write();
  DPbbSSSW_->Write();
  M_othersSSSW_->Write();
  MbbnohSSSW_->Write();
  DPbbnohSSSW_->Write();
  M8SSSW_->Write();
  C8SSSW_->Write();
  M45bestallSSSW_->Write();
  Chi2extallSSSW_->Write();
  DPbballSSSW_->Write();
  N4NJSSSW_->Write();
  E4NJSSSW_->Write();
  
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
DEFINE_FWK_MODULE(TDAna);
