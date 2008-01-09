//////////////////////////////////////////////////////////////////////////////
//
// TDAna.cc
// latest version: 26/12/2007 TD 
//
// To do:
//
// 1) improve likelihood by adding variables included in it
//    (add histograms at SS level, then modify Smooth.C code,
//    run TDAna, run Smooth.C on rootple, run TDAna again with
//    improved functionfileSS.root input)
//    - Variables to add: chi2mass (seems to discriminate better than chi2extall!)  
//    - NJets
//    - GoodHt
//    - M_others
//    - tag information from tag mass
//
// 2) study other kinematical variables and possibly include them:
//    - mass of bb system is not yet very well characterized: need to find 
//      a better way to select the pair
//    - sum of Delta R angles between closest b-tag pairs, E.G. loop on
//      all b-tag jets, compute min(DR) among all, exclude pair, and if
//      at least another pair exists add DR of them; maybe a better possible
//      definition entails computing minimum sum of pair masses. The rationale
//      is to find gluon splitting b-pairs and discriminate from more spread out
//      b-jets as in signal topology
//    - revert logic of M3, M2 selection by first choosing b-pair, then finding
//      best triplet among remaining jets, and look for W boson within those
//    - compute minimum dijet mass in triplet which gives best top mass
//    - compute angle between best triplet momentum and best b-pair momentum
//    - compute mass of those two vectors (just mass of five included jets)
//    - compute average |eta| of first 8 jets
// 
// 3) Compute Q-value of signal vs sum of backgrounds for each variable separately
//    and global likelihood.
//    Then, optimize global likelihood by picking best variables trying to pick
//    those with least correlation (e.g. choose either C6 or C8, not both).
//
// 4) define pseudoexperiments and produce 2-component fits of likelihood, extract
//    distribution of signal significance.
// 
// --------------------------------------------------------------------------------
//
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
#include <iostream>
#include <iomanip>

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
  QCD_( iConfig.getUntrackedParameter<bool> ( "QCD" ) ),
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) )
{

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;
  for ( int i=0; i<10; i++ ) {
    grandtotalh[i]=0;
    grandtotalhpass[i]=0;
    for ( int j=0; j<5; j++ ) {
      total[i][j]=0;
      totalpass[i][j]=0;
    }
  }
  for ( int j=0; j<5; j++ ) {
    grandtotaltt[j]=0;
    grandtotalttpass[j]=0;
  }
  grandgrandtotal=0;
  grandgrandtotalpass=0;

  // Load file with smoothed histograms to use for Likelihood definition
  // -------------------------------------------------------------------
  TFile * FunctionFile = new TFile ("functionfileSS.root");
  FunctionFile->cd();
  HSS_sig[0] = dynamic_cast<TH1D*> ( FunctionFile->Get("C8SS_sigS"));
  HSS_bgr[0] = dynamic_cast<TH1D*> ( FunctionFile->Get("C8SS_bgrS"));
  HSS_sig[1] = dynamic_cast<TH1D*> ( FunctionFile->Get("M45bestallSS_sigS"));
  HSS_bgr[1] = dynamic_cast<TH1D*> ( FunctionFile->Get("M45bestallSS_bgrS"));
  HSS_sig[2] = dynamic_cast<TH1D*> ( FunctionFile->Get("Chi2extallSS_sigS"));
  HSS_bgr[2] = dynamic_cast<TH1D*> ( FunctionFile->Get("Chi2extallSS_bgrS"));
  HSS_sig[3] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSigSS_sigS"));
  HSS_bgr[3] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSigSS_bgrS"));
  HSS_sig[4] = dynamic_cast<TH1D*> ( FunctionFile->Get("M8SS_sigS"));
  HSS_bgr[4] = dynamic_cast<TH1D*> ( FunctionFile->Get("M8SS_bgrS"));
  HSS_sig[5] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSS_sigS"));
  HSS_bgr[5] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSS_bgrS"));
  HSS_sig[6] = dynamic_cast<TH1D*> ( FunctionFile->Get("DP12SS_sigS"));
  HSS_bgr[6] = dynamic_cast<TH1D*> ( FunctionFile->Get("DP12SS_bgrS"));
  HSS_sig[7] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtDP2SS_sigS"));
  HSS_bgr[7] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtDP2SS_bgrS"));
  HSS_sig[8] = dynamic_cast<TH1D*> ( FunctionFile->Get("C6SS_sigS"));
  HSS_bgr[8] = dynamic_cast<TH1D*> ( FunctionFile->Get("C6SS_bgrS"));
  HSS_sig[9] = dynamic_cast<TH1D*> ( FunctionFile->Get("M6SS_sigS"));
  HSS_bgr[9] = dynamic_cast<TH1D*> ( FunctionFile->Get("M6SS_bgrS"));
  HSS_sig[10]= dynamic_cast<TH1D*> ( FunctionFile->Get("CorrSumEtSS_sigS"));
  HSS_bgr[10]= dynamic_cast<TH1D*> ( FunctionFile->Get("CorrSumEtSS_bgrS"));
  //  FunctionFile->Close();

  // File with HED and HPD distributions from QCD jets
  // -------------------------------------------------
  TFile * HEDFile = new TFile ("HED.root");
  HEDFile->cd();
  HEDpdf = dynamic_cast<TH1D*> ( HEDFile->Get("HED"));
  HPDpdf = dynamic_cast<TH1D*> ( HEDFile->Get("HPD"));
  // HEDFile->Close();

  // File for output histograms
  // --------------------------
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");

  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd(); 

  // Histograms
  // ----------
  Nsltt_hjj=0.; // Counter of semileptonic tt decays with h->jj decays
  Drmax_ = new TH2D ( "Drmax", "Drmax for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Drmedall_ = new TH2D ( "Drmedall", "Drmedall for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Drmed07_ = new TH2D ( "Drmed07", "Drmed07 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  N07_ = new TH2D ( "N07", "N07 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  N04_ = new TH2D ( "N04", "N04 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  N02_ = new TH2D ( "N02", "N02 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Nlo_ = new TH2D ( "Nlo", "Nlo for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Detmedall_ = new TH2D ( "Detmedall", "Detmedall for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Detmed07_ = new TH2D ( "Detmed07", "Detmed07 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Perf07_ = new TH2D ( "Perf07", "Perf07 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Perf04_ = new TH2D ( "Perf04", "Perf04 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Perf02_ = new TH2D ( "Perf02", "Perf02 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Det2medall_ = new TH2D ( "Det2medall", "Det2medall for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Det2med07_ = new TH2D ( "Det2med07", "Det2med07 for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6  );
  Hrecfrac_ = new TH2D ( "Hrecfrac", "Hrecfrac for choices of etmin, etamax", 21, 14, 56, 11, 1.4, 3.6 );

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
  M6_ = new TH1D ( "M6", "Mass 6 jets", 50, 0., 2500.); 
  C6_ = new TH1D ( "C6", "Centrality 6 jets", 50, 0., 1.2 );
  M8_ = new TH1D ( "M8", "Mass 8 jets", 50, 0., 2500.);
  C8_ = new TH1D ( "C8", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestall_ = new TH1D ( "M45bestall", "M45best", 50, 0, 300 );
  Chi2extall_ = new TH1D ( "Chi2extall", "Chi2ext", 50, 0, 50 );
  DPbball_ = new TH1D ( "DPbball", "Delta phi bb", 50, 0., 3.2 );
  SumHED4_ = new TH1D ( "SumHED4", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4_ = new TH1D ( "SumHPD4", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6_ = new TH1D ( "SumHED6", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6_ = new TH1D ( "SumHPD6", "Sum of HPD in 6 jets", 50, 0., 50. );
  HED_ = new TH1D ( "HED", "High eff discriminator", 200, 0., 40. );      
  HPD_ = new TH1D ( "HPD", "High pur discriminator", 200, 0., 40. );      

  NJetsN_ = new TH1D ( "NJetsN", "Number of selected jets", 20, 0, 20 );
  UncorrHtN_ = new TH1D ( "UncorrHtN", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtN_ = new TH1D ( "CorrHtN", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtN_ = new TH1D ( "GoodHtN", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2N_ = new TH1D ( "GoodHt2N", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtN_ = new TH1D ( "UncorrSumEtN", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtN_ = new TH1D ( "CorrSumEtN", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtN_ = new TH1D ( "GoodSumEtN", "SumEt with selected jets", 50, 0, 4000 );
  MEtN_ = new TH1D ( "MEtN", "MEt", 50, 0, 500 );
  MEtSigN_ = new TH1D ( "MEtSigN", "MEtSig", 50, 0, 20 );
  MEtSigNewN_ = new TH1D ( "MEtSigNewN", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMN_ = new TH1D ( "MEtDPMN", "MEtDPM", 50, 0, 3.2 );
  MEtDP1N_ = new TH1D ( "MEtDP1N", "MEtDP1", 50, 0, 3.2 );
  MEtDP2N_ = new TH1D ( "MEtDP2N", "MEtDP2", 50, 0, 3.2 );
  MEtDP3N_ = new TH1D ( "MEtDP3N", "MEtDP3", 50, 0, 3.2 );
  UncorrMEtSigN_ = new TH1D ( "UncorrMEtSigN", "MEtSig", 50, 0, 20 );
  CorrMEtSigN_ = new TH1D ( "CorrMEtSigN", "MEtSig", 50, 0, 20 );
  M3bestN_ = new TH1D ( "M3bestN", "M3best", 50, 0, 300 );
  MwbestN_ = new TH1D ( "MwbestN", "Mwbest", 50, 0, 200 );
  Chi2massN_ = new TH1D ( "Chi2massN", "Chi2mass", 50, 0, 50 );
  M45bestN_ = new TH1D ( "M45bestN", "M45best", 50, 0, 300 );
  Chi2extN_ = new TH1D ( "Chi2extN", "Chi2ext", 50, 0, 50 );
  MEx_SumEtN_ = new TH2D ( "MEx_SumEtN", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12N_ = new TH1D ( "DP12N", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbN_ = new TH1D ( "DPbbN", "Delta phi bb", 50, 0., 3.2 );
  M_othersN_ = new TH1D ( "M_othersN", "Mass of other jets", 50, 0., 1500. );
  MbbnohN_ = new TH1D ( "MbbnohN", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohN_ = new TH1D ( "DPbbnohN", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M6N_ = new TH1D ( "M6N", "Mass 6 jets", 50, 0., 2500.);
  C6N_ = new TH1D ( "C6N", "Centrality 6 jets", 50, 0., 1.2 );
  M8N_ = new TH1D ( "M8N", "Mass 8 jets", 50, 0., 2500.);
  C8N_ = new TH1D ( "C8N", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallN_ = new TH1D ( "M45bestallN", "M45best", 50, 0, 300 );
  Chi2extallN_ = new TH1D ( "Chi2extallN", "Chi2ext", 50, 0, 50 );
  DPbballN_ = new TH1D ( "DPbballN", "Delta phi bb", 50, 0., 3.2 );
  SumHED4N_ = new TH1D ( "SumHED4N", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4N_ = new TH1D ( "SumHPD4N", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6N_ = new TH1D ( "SumHED6N", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6N_ = new TH1D ( "SumHPD6N", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDN_ = new TH1D ( "HEDN", "High eff discriminator", 200, 0., 40. );      
  HPDN_ = new TH1D ( "HPDN", "High pur discriminator", 200, 0., 40. );      

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
  M6S_ = new TH1D ( "M6S", "Mass 6 jets", 50, 0., 2500.);
  C6S_ = new TH1D ( "C6S", "Centrality 6 jets", 50, 0., 1.2 );
  M8S_ = new TH1D ( "M8S", "Mass 8 jets", 50, 0., 2500.);
  C8S_ = new TH1D ( "C8S", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallS_ = new TH1D ( "M45bestallS", "M45best", 50, 0, 300 );
  Chi2extallS_ = new TH1D ( "Chi2extallS", "Chi2ext", 50, 0, 50 );
  DPbballS_ = new TH1D ( "DPbballS", "Delta phi bb", 50, 0., 3.2 );
  SumHED4S_ = new TH1D ( "SumHED4S", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4S_ = new TH1D ( "SumHPD4S", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6S_ = new TH1D ( "SumHED6S", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6S_ = new TH1D ( "SumHPD6S", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDS_ = new TH1D ( "HEDS", "High eff discriminator", 200, 0., 40. );       
  HPDS_ = new TH1D ( "HPDS", "High eff discriminator", 200, 0., 40. );       

  NJetsSN_ = new TH1D ( "NJetsSN", "Number of selected jets", 20, 0, 20 );
  UncorrHtSN_ = new TH1D ( "UncorrHtSN", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSN_ = new TH1D ( "CorrHtSN", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSN_ = new TH1D ( "GoodHtSN", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SN_ = new TH1D ( "GoodHt2SN", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSN_ = new TH1D ( "UncorrSumEtSN", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSN_ = new TH1D ( "CorrSumEtSN", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSN_ = new TH1D ( "GoodSumEtSN", "SumEt with selected jets", 50, 0, 4000 );
  MEtSN_ = new TH1D ( "MEtSN", "MEt", 50, 0, 500 );
  MEtSigSN_ = new TH1D ( "MEtSigSN", "MEtSig", 50, 0, 20 );
  MEtSigNewSN_ = new TH1D ( "MEtSigNewSN", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSN_ = new TH1D ( "MEtDPMSN", "MEtDPM", 50, 0, 3.2 );
  MEtDP1SN_ = new TH1D ( "MEtDP1SN", "MEtDP1S", 50, 0, 3.2 );
  MEtDP2SN_ = new TH1D ( "MEtDP2SN", "MEtDP2S", 50, 0, 3.2 );
  MEtDP3SN_ = new TH1D ( "MEtDP3SN", "MEtDP3S", 50, 0, 3.2 );
  UncorrMEtSigSN_ = new TH1D ( "UncorrMEtSigSN", "MEtSig", 50, 0, 20 );
  CorrMEtSigSN_ = new TH1D ( "CorrMEtSigSN", "MEtSig", 50, 0, 20 );
  M3bestSN_ = new TH1D ( "M3bestSN", "M3bestS", 50, 0, 300 );
  MwbestSN_ = new TH1D ( "MwbestSN", "MwbestS", 50, 0, 200 );
  Chi2massSN_ = new TH1D ( "Chi2massSN", "Chi2massS", 50, 0, 50 );
  M45bestSN_ = new TH1D ( "M45bestSN", "M45best", 50, 0, 300 );
  Chi2extSN_ = new TH1D ( "Chi2extSN", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSN_ = new TH2D ( "MEx_SumEtSN", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SN_ = new TH1D ( "DP12SN", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSN_ = new TH1D ( "DPbbSN", "Delta phi bb", 50, 0., 3.2 );
  M_othersSN_ = new TH1D ( "M_othersSN", "Mass of other jets", 50, 0., 1500. );
  MbbnohSN_ = new TH1D ( "MbbnohSN", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSN_ = new TH1D ( "DPbbnohSN", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M6SN_ = new TH1D ( "M6SN", "Mass 6 jets", 50, 0., 2500.);
  C6SN_ = new TH1D ( "C6SN", "Centrality 6 jets", 50, 0., 1.2 );
  M8SN_ = new TH1D ( "M8SN", "Mass 8 jets", 50, 0., 2500.);
  C8SN_ = new TH1D ( "C8SN", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSN_ = new TH1D ( "M45bestallSN", "M45best", 50, 0, 300 );
  Chi2extallSN_ = new TH1D ( "Chi2extallSN", "Chi2ext", 50, 0, 50 );
  DPbballSN_ = new TH1D ( "DPbballSN", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SN_ = new TH1D ( "SumHED4SN", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SN_ = new TH1D ( "SumHPD4SN", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SN_ = new TH1D ( "SumHED6SN", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SN_ = new TH1D ( "SumHPD6SN", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSN_ = new TH1D ( "HEDSN", "High eff discriminator", 200, 0., 40. );     
  HPDSN_ = new TH1D ( "HPDSN", "High pur discriminator", 200, 0., 40. );     

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
  M6SS_ = new TH1D ( "M6SS", "Mass 6 jets", 50, 0., 2500.);
  C6SS_ = new TH1D ( "C6SS", "Centrality 6 jets", 50, 0., 1.2 );
  M8SS_ = new TH1D ( "M8SS", "Mass 8 jets", 50, 0., 2500.);
  C8SS_ = new TH1D ( "C8SS", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSS_ = new TH1D ( "M45bestallSS", "M45best", 50, 0, 300 );
  Chi2extallSS_ = new TH1D ( "Chi2extallSS", "Chi2ext", 50, 0, 50 );
  DPbballSS_ = new TH1D ( "DPbballSS", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SS_ = new TH1D ( "SumHED4SS", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SS_ = new TH1D ( "SumHPD4SS", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SS_ = new TH1D ( "SumHED6SS", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SS_ = new TH1D ( "SumHPD6SS", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSS_ = new TH1D ( "HEDSS", "High eff discriminator", 200, 0., 40. );     
  HPDSS_ = new TH1D ( "HPDSS", "High eff discriminator", 200, 0., 40. );     

  NJetsSSN_ = new TH1D ( "NJetsSSN", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSN_ = new TH1D ( "UncorrHtSSN", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSN_ = new TH1D ( "CorrHtSSN", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSN_ = new TH1D ( "GoodHtSSN", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSN_ = new TH1D ( "GoodHt2SSN", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSN_ = new TH1D ( "UncorrSumEtSSN", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSN_ = new TH1D ( "CorrSumEtSSN", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSN_ = new TH1D ( "GoodSumEtSSN", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSN_ = new TH1D ( "MEtSSN", "MEt", 50, 0, 500 );
  MEtSigSSN_ = new TH1D ( "MEtSigSSN", "MEtSig", 50, 0, 20 );
  MEtSigNewSSN_ = new TH1D ( "MEtSigNewSSN", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSN_ = new TH1D ( "MEtDPMSSN", "MEtDPMSSN", 50, 0, 3.2 );
  MEtDP1SSN_ = new TH1D ( "MEtDP1SSN", "MEtDP1SSN", 50, 0, 3.2 );
  MEtDP2SSN_ = new TH1D ( "MEtDP2SSN", "MEtDP2SSN", 50, 0, 3.2 );
  MEtDP3SSN_ = new TH1D ( "MEtDP3SSN", "MEtDP3SSN", 50, 0, 3.2 );
  UncorrMEtSigSSN_ = new TH1D ( "UncorrMEtSigSSN", "UncorrMEtSigSSN", 50, 0, 20 );
  CorrMEtSigSSN_ = new TH1D ( "CorrMEtSigSSN", "CorrMEtSigSSN", 50, 0, 20 );
  M3bestSSN_ = new TH1D ( "M3bestSSN", "M3bestSSN", 50, 0, 300 );
  MwbestSSN_ = new TH1D ( "MwbestSSN", "MwbestSSN", 50, 0, 200 );
  Chi2massSSN_ = new TH1D ( "Chi2massSSN", "Chi2massSSN", 50, 0, 50 );
  M45bestSSN_ = new TH1D ( "M45bestSSN", "M45best", 50, 0, 300 );
  Chi2extSSN_ = new TH1D ( "Chi2extSSN", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSN_ = new TH2D ( "MEx_SumEtSSN", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSN_ = new TH1D ( "DP12SSN", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSN_ = new TH1D ( "DPbbSSN", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSN_ = new TH1D ( "M_othersSSN", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSN_ = new TH1D ( "MbbnohSSN", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSSN_ = new TH1D ( "DPbbnohSSN", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M6SSN_ = new TH1D ( "M6SSN", "Mass 6 jets", 50, 0., 2500.);
  C6SSN_ = new TH1D ( "C6SSN", "Centrality 6 jets", 50, 0., 1.2 );
  M8SSN_ = new TH1D ( "M8SSN", "Mass 8 jets", 50, 0., 2500.);
  C8SSN_ = new TH1D ( "C8SSN", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSN_ = new TH1D ( "M45bestallSSN", "M45best", 50, 0, 300 );
  Chi2extallSSN_ = new TH1D ( "Chi2extallSSN", "Chi2ext", 50, 0, 50 );
  DPbballSSN_ = new TH1D ( "DPbballSSN", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SSN_ = new TH1D ( "SumHED4SSN", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SSN_ = new TH1D ( "SumHPD4SSN", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SSN_ = new TH1D ( "SumHED6SSN", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SSN_ = new TH1D ( "SumHPD6SSN", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSSN_ = new TH1D ( "HEDSSN", "High eff discriminator", 200, 0., 40. );   
  HPDSSN_ = new TH1D ( "HPDSSN", "High pur discriminator", 200, 0., 40. );   

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
  M6SSS_ = new TH1D ( "M6SSS", "Mass 6 jets", 50, 0., 2500.);
  C6SSS_ = new TH1D ( "C6SSS", "Centrality 6 jets", 50, 0., 1.2 );
  M8SSS_ = new TH1D ( "M8SSS", "Mass 8 jets", 50, 0., 2500.);
  C8SSS_ = new TH1D ( "C8SSS", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSS_ = new TH1D ( "M45bestallSSS", "M45best", 50, 0, 300 );
  Chi2extallSSS_ = new TH1D ( "Chi2extallSSS", "Chi2ext", 50, 0, 50 );
  DPbballSSS_ = new TH1D ( "DPbballSSS", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SSS_ = new TH1D ( "SumHED4SSS", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SSS_ = new TH1D ( "SumHPD4SSS", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SSS_ = new TH1D ( "SumHED6SSS", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SSS_ = new TH1D ( "SumHPD6SSS", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSSS_ = new TH1D ( "HEDSSS_", "High eff discriminator", 200, 0., 40. );  
  HPDSSS_ = new TH1D ( "HPDSSS_", "High eff discriminator", 200, 0., 40. );  

  NJetsSSSN_ = new TH1D ( "NJetsSSSN", "Number of selected jets", 20, 0, 20 );
  UncorrHtSSSN_ = new TH1D ( "UncorrHtSSSN", "Ht with uncorrected jets", 50, 0, 4000 );
  CorrHtSSSN_ = new TH1D ( "CorrHtSSSN", "Ht with corrected jets", 50, 0, 4000 );
  GoodHtSSSN_ = new TH1D ( "GoodHtSSSN", "Ht with corrected good jets", 50, 0, 4000 );
  GoodHt2SSSN_ = new TH1D ( "GoodHt2SSSN", "Ht with corrected good jets and corr MEt", 50, 0, 4000 );
  UncorrSumEtSSSN_ = new TH1D ( "UncorrSumEtSSSN", "SumEt with uncorrected jets", 50, 0, 4000 );
  CorrSumEtSSSN_ = new TH1D ( "CorrSumEtSSSN", "SumEt with corrected jets", 50, 0, 4000 );
  GoodSumEtSSSN_ = new TH1D ( "GoodSumEtSSSN", "SumEt with selected jets", 50, 0, 4000 );
  MEtSSSN_ = new TH1D ( "MEtSSSN", "MEt", 50, 0, 500 );
  MEtSigSSSN_ = new TH1D ( "MEtSigSSSN", "MEtSig", 50, 0, 20 );
  MEtSigNewSSSN_ = new TH1D ( "MEtSigNewSSSN", "MEtSig with tuned resolution", 50, 0., 20. );
  MEtDPMSSSN_ = new TH1D ( "MEtDPMSSSN", "MEtDPMSSSN", 50, 0, 3.2 );
  MEtDP1SSSN_ = new TH1D ( "MEtDP1SSSN", "MEtDP1SSSN", 50, 0, 3.2 );
  MEtDP2SSSN_ = new TH1D ( "MEtDP2SSSN", "MEtDP2SSSN", 50, 0, 3.2 );
  MEtDP3SSSN_ = new TH1D ( "MEtDP3SSSN", "MEtDP3SSSN", 50, 0, 3.2 );
  UncorrMEtSigSSSN_ = new TH1D ( "UncorrMEtSigSSSN", "UncorrMEtSigSSSN", 50, 0, 20 );
  CorrMEtSigSSSN_ = new TH1D ( "CorrMEtSigSSSN", "CorrMEtSigSSSN", 50, 0, 20 );
  M3bestSSSN_ = new TH1D ( "M3bestSSSN", "M3bestSSSN", 50, 0, 300 );
  MwbestSSSN_ = new TH1D ( "MwbestSSSN", "MwbestSSSN", 50, 0, 200 );
  Chi2massSSSN_ = new TH1D ( "Chi2massSSSN", "Chi2massSSSN", 50, 0, 50 );
  M45bestSSSN_ = new TH1D ( "M45bestSSSN", "M45best", 50, 0, 300 );
  Chi2extSSSN_ = new TH1D ( "Chi2extSSSN", "Chi2ext", 50, 0, 50 );
  MEx_SumEtSSSN_ = new TH2D ( "MEx_SumEtSSSN", "MEx vs SumEt", 50, 0., 2000., 50, -200., 200. );
  DP12SSSN_ = new TH1D ( "DP12SSSN", "Delta phi jet 1 jet 2", 50, 0., 3.2 );
  DPbbSSSN_ = new TH1D ( "DPbbSSSN", "Delta phi bb", 50, 0., 3.2 );
  M_othersSSSN_ = new TH1D ( "M_othersSSSN", "Mass of other jets", 50, 0., 1500. );
  MbbnohSSSN_ = new TH1D ( "MbbnohSSSN", "Mass of b pair not from h", 50, 0., 800 );
  DPbbnohSSSN_ = new TH1D ( "DPbbnohSSSN", "Delta phi of b pair not from h", 50, 0., 3.2 );
  M6SSSN_ = new TH1D ( "M6SSSN", "Mass 6 jets", 50, 0., 2500.);
  C6SSSN_ = new TH1D ( "C6SSSN", "Centrality 6 jets", 50, 0., 1.2 );
  M8SSSN_ = new TH1D ( "M8SSSN", "Mass 8 jets", 50, 0., 2500.);
  C8SSSN_ = new TH1D ( "C8SSSN", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSSN_ = new TH1D ( "M45bestallSSSN", "M45best", 50, 0, 300 );
  Chi2extallSSSN_ = new TH1D ( "Chi2extallSSSN", "Chi2ext", 50, 0, 50 );
  DPbballSSSN_ = new TH1D ( "DPbballSSSN", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SSSN_ = new TH1D ( "SumHED4SSSN", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SSSN_ = new TH1D ( "SumHPD4SSSN", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SSSN_ = new TH1D ( "SumHED6SSSN", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SSSN_ = new TH1D ( "SumHPD6SSSN", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSSSN_ = new TH1D ( "HEDSSSN", "High eff discriminator", 200, 0., 40. ); 
  HPDSSSN_ = new TH1D ( "HPDSSSN", "High pur discriminator", 200, 0., 40. ); 

  N4NJSSS_ = new TH1D ( "N4NJSSS", "N of 4HEL tags vs N jets", 
			20, 0, 20 );  // These are filled only for this
  E4NJSSS_ = new TH1D ( "E4NJSSS", "Efficiency of 4HEL tags vs N jets", 
			20, 0, 20 );  // selection and do not require W, N

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
  M6W_ = new TH1D ( "M6W", "Mass 6 jets", 50, 0., 2500.);
  C6W_ = new TH1D ( "C6W", "Centrality 6 jets", 50, 0., 1.2 );
  M8W_ = new TH1D ( "M8W", "Mass 8 jets", 50, 0., 2500.);
  C8W_ = new TH1D ( "C8W", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallW_ = new TH1D ( "M45bestallW", "M45best", 50, 0, 300 );
  Chi2extallW_ = new TH1D ( "Chi2extallW", "Chi2ext", 50, 0, 50 );
  DPbballW_ = new TH1D ( "DPbballW", "Delta phi bb", 50, 0., 3.2 );
  SumHED4W_ = new TH1D ( "SumHED4W", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4W_ = new TH1D ( "SumHPD4W", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6W_ = new TH1D ( "SumHED6W", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6W_ = new TH1D ( "SumHPD6W", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDW_ = new TH1D ( "HEDW", "High eff discriminator", 200, 0., 40. );     
  HPDW_ = new TH1D ( "HPDW", "High pur discriminator", 200, 0., 40. );     

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
  M6SW_ = new TH1D ( "M6SW", "Mass 6 jets", 50, 0., 2500.);
  C6SW_ = new TH1D ( "C6SW", "Centrality 6 jets", 50, 0., 1.2 );
  M8SW_ = new TH1D ( "M8SW", "Mass 8 jets", 50, 0., 2500.);
  C8SW_ = new TH1D ( "C8SW", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSW_ = new TH1D ( "M45bestallSW", "M45best", 50, 0, 300 );
  Chi2extallSW_ = new TH1D ( "Chi2extallSW", "Chi2ext", 50, 0, 50 );
  DPbballSW_ = new TH1D ( "DPbballSW", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SW_ = new TH1D ( "SumHED4SW", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SW_ = new TH1D ( "SumHPD4SW", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SW_ = new TH1D ( "SumHED6SW", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SW_ = new TH1D ( "SumHPD6SW", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSW_ = new TH1D ( "HEDSW", "High eff discriminator", 200, 0., 40. );   
  HPDSW_ = new TH1D ( "HPDSW", "High pur discriminator", 200, 0., 40. );   

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
  M6SSW_ = new TH1D ( "M6SSW", "Mass 6 jets", 50, 0., 2500.);
  C6SSW_ = new TH1D ( "C6SSW", "Centrality 6 jets", 50, 0., 1.2 );
  M8SSW_ = new TH1D ( "M8SSW", "Mass 8 jets", 50, 0., 2500.);
  C8SSW_ = new TH1D ( "C8SSW", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSW_ = new TH1D ( "M45bestallSSW", "M45best", 50, 0, 300 );
  Chi2extallSSW_ = new TH1D ( "Chi2extallSSW", "Chi2ext", 50, 0, 50 );
  DPbballSSW_ = new TH1D ( "DPbballSSW", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SSW_ = new TH1D ( "SumHED4SSW", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SSW_ = new TH1D ( "SumHPD4SSW", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SSW_ = new TH1D ( "SumHED6SSW", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SSW_ = new TH1D ( "SumHPD6SSW", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSSW_ = new TH1D ( "HEDSSW", "High eff discriminator", 200, 0., 40. ); 
  HPDSSW_ = new TH1D ( "HPDSSW", "High pur discriminator", 200, 0., 40. ); 

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
  M6SSSW_ = new TH1D ( "M6SSSW", "Mass 6 jets", 50, 0., 2500.);
  C6SSSW_ = new TH1D ( "C6SSSW", "Centrality 6 jets", 50, 0., 1.2 );
  M8SSSW_ = new TH1D ( "M8SSSW", "Mass 8 jets", 50, 0., 2500.);
  C8SSSW_ = new TH1D ( "C8SSSW", "Centrality 8 jets", 50, 0., 1.2 );
  M45bestallSSSW_ = new TH1D ( "M45bestallSSSW", "M45best", 50, 0, 300 );
  Chi2extallSSSW_ = new TH1D ( "Chi2extallSSSW", "Chi2ext", 50, 0, 50 );
  DPbballSSSW_ = new TH1D ( "DPbballSSSW", "Delta phi bb", 50, 0., 3.2 );
  SumHED4SSSW_ = new TH1D ( "SumHED4SSSW", "Sum of HED in 4 jets", 50, 0., 100. );
  SumHPD4SSSW_ = new TH1D ( "SumHPD4SSSW", "Sum of HPD in 4 jets", 50, 0., 50. );
  SumHED6SSSW_ = new TH1D ( "SumHED6SSSW", "Sum of HED in 6 jets", 50, 0., 100. );
  SumHPD6SSSW_ = new TH1D ( "SumHPD6SSSW", "Sum of HPD in 6 jets", 50, 0., 50. );
  HEDSSSW_ = new TH1D ( "HEDSSSW", "High eff discriminator", 200, 0., 40. ); 
  HPDSSSW_ = new TH1D ( "HPDSSSW", "High pur discriminator", 200, 0., 40. ); 

  N4NJSSSW_ = new TH1D ( "N4NJSSSW", "N of 4HEL tags vs N jets", 20, 0, 20 );
  E4NJSSSW_ = new TH1D ( "E4NJSSSW", "Efficiency of 4HEL tags vs N jets", 20, 0, 20 );

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

  // Relative Likelihood
  // -------------------       
  rel_lik = 0.;

  L_    = new TH1D ("L", "Likelihood 2t passing sel", 50, -10., 10.  );
  LS_   = new TH1D ("LS", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSS_  = new TH1D ("LSS", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSSS_ = new TH1D ("LSSS", "Likelihood 2t passing sel", 50, -10., 10.  );
  LW_    = new TH1D ("LW", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSW_   = new TH1D ("LSW", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSSW_  = new TH1D ("LSSW", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSSSW_ = new TH1D ("LSSSW", "Likelihood 2t passing sel", 50, -10., 10.  );
  LN_    = new TH1D ("LN", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSN_   = new TH1D ("LSN", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSSN_  = new TH1D ("LSSN", "Likelihood 2t passing sel", 50, -10., 10.  );
  LSSSN_ = new TH1D ("LSSSN", "Likelihood 2t passing sel", 50, -10., 10.  );

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

  cout << "Done with constructor" << endl;

}


TDAna::~TDAna()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  // Write out the file
  // ------------------
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
  

  /////////////////////////////////// MC analysis ////////////////////////////

  // Compute events in each final state
  // ----------------------------------
  // Take genParticleCandidates
  // --------------------------
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  // Take the mothers
  // ----------------
  int hdecay=0;
  int tdecay=0;
  int foundw=0;
  double Parton_pt[8]={0.};
  double Parton_eta[8]={0.};
  double Parton_phi[8]={0.};
  double Parton_mass[8]={0.};
  double Parton_dec[8]={0.};
  int iparton=0;

  // NB we know that the block is filled with top and antitop first, then Higgs
  // So we need not worry about messing up W's from top and from H ...
  // --------------------------------------------------------------------------
  for ( MCParticleCollection::const_iterator MCp = MCpartons->begin(); MCp != MCpartons->end(); ++MCp ) {
    // See if this is a daughter of the W
    // ----------------------------------
    if ( fabs(MCp->mPid())==24 && foundw<4 ) { // first two W's in list are always from t, tbar
      if ( MCp->pid()<6 && MCp->pid()>0 ) tdecay+=1; // count particles, not antiparticles, to count once
      if ( fabs(MCp->pid())==11 ) tdecay+=2;
      if ( fabs(MCp->pid())==13 ) tdecay+=4;
      if ( fabs(MCp->pid())==15 ) tdecay+=8;
      foundw++;
      
      // For hadronic W decay from top, store info on partons
      // ----------------------------------------------------
      if ( fabs(MCp->pid())<6 && iparton<8 ) {
	Parton_pt[iparton]=MCp->pt();
	Parton_eta[iparton]=MCp->eta();
	Parton_phi[iparton]=MCp->phi();
	Parton_mass[iparton]=MCp->mass();
	if ( MCp->mPid()==24 )  Parton_dec[iparton] = 1; // 1 for W+ daughters, 2 for b from top
	if ( MCp->mPid()==-24 ) Parton_dec[iparton] = 3; // 3 for W- daughters, 4 for bbar frop antitop
	iparton++;
      }
    }
    // Store information on b partons from t, tbar
    // -------------------------------------------
    if ( fabs(MCp->mPid())==6 && iparton<8 ) {
      if ( MCp->pid()==5 ) {
	Parton_pt[iparton]=MCp->pt();
	Parton_eta[iparton]=MCp->eta();
	Parton_phi[iparton]=MCp->phi();
	Parton_mass[iparton]=MCp->mass();
	Parton_dec[iparton]=2; // 2 for b from top
	iparton++;
      } else if ( MCp->pid()==-5 ) {
	Parton_pt[iparton]=MCp->pt();
	Parton_eta[iparton]=MCp->eta();
	Parton_phi[iparton]=MCp->phi();
	Parton_mass[iparton]=MCp->mass();
	Parton_dec[iparton]=4; // 4 for bbar from antitop
	iparton++;
      }
    }
    if ( MCp->mPid()==25 ) {
      // See if this is a daughter of the Higgs
      // --------------------------------------
      if ( fabs(MCp->pid())== 5 ) hdecay=1;
      if ( fabs(MCp->pid())== 4 ) hdecay=2;
      if ( fabs(MCp->pid())== 15 ) hdecay=3;
      if ( fabs(MCp->pid())== 24 ) hdecay=4;
      
      // Store information on partons from H
      // -----------------------------------
      if ( fabs(MCp->pid())<6 && iparton<8 ) {
	Parton_pt[iparton]=MCp->pt();
	Parton_eta[iparton]=MCp->eta();
	Parton_phi[iparton]=MCp->phi();
	Parton_mass[iparton]=MCp->mass();
	Parton_dec[iparton]=5; // 5 for H daughters
	iparton++;	
      }
    }
  }
  int itdecay = 0;
  if ( tdecay==2 ) itdecay=0; // 6j
  if ( tdecay==3 ) itdecay=1; // en4j
  if ( tdecay==5 ) itdecay=2; // mn4j
  if ( tdecay==9 ) itdecay=3; // tn4j
  if ( tdecay==4 ) itdecay=4; // enen2j
  if ( tdecay==6 ) itdecay=5; // enmn2j
  if ( tdecay==10) itdecay=6; // entn2j
  if ( tdecay==8 ) itdecay=7; // mnmn2j
  if ( tdecay==12) itdecay=8; // mntn2j
  if ( tdecay==16) itdecay=9; // tntn2j
  if ( itdecay==-1 ) {
    itdecay=0;
    cout << "Wrong decay assignment" << endl;
  }
  total[itdecay][hdecay]++;
  grandtotaltt[hdecay]++;
  grandtotalh[itdecay]++;
  grandgrandtotal++;
  if ( response ) {
    totalpass[itdecay][hdecay]++;
    grandtotalttpass[hdecay]++;
    grandtotalhpass[itdecay]++;
    grandgrandtotalpass++;
  }

  // HiVariables
  // -----------
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );
  
  // Preliminary creation of jet array ordered by Et for parton-jet matching study only
  // ----------------------------------------------------------------------------------
  if ( ( itdecay==1 || itdecay==2 || itdecay==3 ) && 
       ( hdecay==1 || hdecay==2 ) && 
       caloJets->size() != 0 && response && iparton==6 ) { // sl decay + h->cc,bb & jets exist
    double JPet[200]={0.};
    double JPeta[200]={0.};
    double JPphi[200]={0.};
    int NJP=0;
    for ( OfflineJetCollection::const_iterator cal = caloJets->begin(); 
	  cal != caloJets->end(); ++cal ) {
      if ( NJP<200 ) {
	JPet[NJP]=cal->et();
	JPeta[NJP]=cal->eta();
	JPphi[NJP]=cal->phi();
	NJP++;
      }
    }
    // Order arrays by Et
    // ------------------
    for ( int times=0; times<NJP; times++ ) {
      for ( int j=0; j<NJP-1; j++ ) {
	int k = NJP-j-1;
	if ( JPet[k]>JPet[k-1] ) {
	  double a1 = JPet[k];
	  double a2 = JPeta[k];
	  double a3 = JPphi[k];
	  JPet[k]=JPet[k-1];
	  JPeta[k]=JPeta[k-1];
	  JPphi[k]=JPphi[k-1];
	  JPet[k-1]=a1;
	  JPeta[k-1]=a2;
	  JPphi[k-1]=a3;
	}
      }
    }
    // Study parton-jet matching
    // -------------------------
    bool ipass[8]; // whether a parton is associated to one of the leading jets
    for ( int ietmin=0; ietmin<21; ietmin++ ) {
      double etmin = 15.+2.*(double)ietmin;
      for ( int ietamax=0; ietamax<11; ietamax++ ) {
	double etamax = 1.5+0.2*(double)ietamax;
	// Variables needed to fill histograms
	// -----------------------------------
	double drmax=0.;
	double drmed07=0.;
	double drmedall=0.;
	double n07=0.;
	double n04=0.;
	double n02=0.;
	double leftover=0.;
	double detmed07=0.;
	double detmedall=0.;
	double det2med07=0.;
	double det2medall=0.;
	double nptotal=0.;
	int ijh[2]={0};
	int njh=0.;
	for ( int ip=0; ip<8; ip++ ) { ipass[ip]=false; }; 
	// Loop on jets
	// ------------
	int considered=0; // Consider at most 8 jets for the association
	for ( int ij=0; ij<NJP; ij++ ) {
	  if ( JPet[ij]>etmin && fabs(JPeta[ij])<etamax && considered<8 ) {
	    considered++;
	    double drmin=1.0;
	    double det=0.;
	    double det2=0.;
	    int tmpind=-1;
	    for ( int ip=0; ip<iparton; ip++ ) {
	      if ( !ipass[ip] ) { // this parton has not been used yet
		double deta = Parton_eta[ip]-JPeta[ij];
		double dphi = 3.1415926-fabs(fabs(Parton_phi[ip]-JPphi[ij])-3.1415926);
		double d_et = fabs(Parton_pt[ip]-JPet[ij])/Parton_pt[ip];
		double dr2  = deta*deta+dphi*dphi; // +d_et*d_et;
		if ( dr2<drmin ) {
		  drmin=dr2;
		  tmpind=ip;
		  det=(-Parton_pt[ip]+JPet[ij])/Parton_pt[ip];
		  det2=pow((Parton_pt[ip]-JPet[ij])/Parton_pt[ip],2);
		}
	      }
	    }
	    drmin = sqrt(drmin);
	    if ( tmpind>=0 ) { 
	      ipass[tmpind]=true;
	      nptotal++;
	      // Compute estimators of association
	      // ---------------------------------
	      if ( drmin>drmax ) drmax=drmin;
	      drmedall+=drmin;
	      detmedall+=det;
	      det2medall+=det2;
	      if ( drmin<0.7 ) { 
		drmed07+=drmin;
		n07++;
		detmed07+=det;
		det2med07+=det2;
	      }
	      if ( drmin<0.4 ) n04++;
	      if ( drmin<0.2 ) n02++;
	      if ( Parton_dec[tmpind]==5 ) {
		// Special treatment for H daughters
		// ---------------------------------
		ijh[njh]=ij;
		njh++;
	      }
	    } else { 
	      // Count here how many jets among the first N (at most equal to the number of partons)
	      // could not be associated reasonably to a parton
	      // -----------------------------------------------------------------------------------
	      if ( considered<iparton ) leftover++;
	    }
	  } // if there is a jet
	} // ij
	if ( ietmin+ietamax==0 ) cout << "Total associated partons = " << nptotal << " N07=" << n07 << endl;
	if ( nptotal>0 ) { // total number of partons matched to one of the leading 8 jets
	  drmedall=drmedall/nptotal;
	  detmedall=detmedall/nptotal;
	  det2medall=det2medall/nptotal;
	}
	if ( n07>0 ) {
	  drmed07=drmed07/n07;
	  detmed07=detmed07/n07;
	  det2med07=det2med07/n07;
	}
	// Mass of two jets from h
	// -----------------------
	double m45=0.;
	if ( njh==2 ) {
	  double px4=JPet[ijh[0]]*cos(JPphi[ijh[0]]);
	  double px5=JPet[ijh[1]]*cos(JPphi[ijh[1]]);
	  double py4=JPet[ijh[0]]*sin(JPphi[ijh[0]]);
	  double py5=JPet[ijh[1]]*sin(JPphi[ijh[1]]);
	  double pz4=JPet[ijh[0]]/tan(2*atan(exp(-JPeta[ijh[0]])));
	  double pz5=JPet[ijh[1]]/tan(2*atan(exp(-JPeta[ijh[1]])));
	  double e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
	  double e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
	  m45=(e4+e5)*(e4+e5)-(px4+px5)*(px4+px5)-(py4+py5)*(py4+py5)-(pz4+pz5)*(pz4+pz5);
	  if ( m45>0 ) m45=sqrt(m45);
	}
	// Fill histograms
	// ---------------
	Drmedall_->Fill(etmin, etamax, drmedall);
	Drmed07_->Fill(etmin, etamax, drmed07);
	Drmax_->Fill(etmin, etamax, drmax);
	N07_->Fill(etmin, etamax, n07);
	N04_->Fill(etmin, etamax, n04);
	N02_->Fill(etmin, etamax, n02);
	Nlo_->Fill(etmin, etamax, leftover);
	Detmed07_->Fill(etmin, etamax, detmed07);
	Detmedall_->Fill(etmin, etamax, detmedall);
	Det2med07_->Fill(etmin, etamax, det2med07);
	Det2medall_->Fill(etmin, etamax, det2medall);
	if ( n07==iparton ) Perf07_->Fill(etmin, etamax);
	if ( n04==iparton ) Perf04_->Fill(etmin, etamax);
	if ( n02==iparton ) Perf02_->Fill(etmin, etamax);
	if ( fabs(m45-120)<20 ) Hrecfrac_->Fill(etmin,etamax);
      } // ietamax
    } // ietmin
    Nsltt_hjj++;
  } // if sl tt decay + h->cc,bb decay
  
  ////////////////////////////////////////////////////////////////////////////////////
  
  // Offline analysis
  // ----------------

  // MEt
  // ---
  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  double sumet = caloMET->sumEt();
  double met = caloMET->et();
  double metsig = caloMET->mEtSig();
  double metphi = caloMET->phi();

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

  double m6=0.;
  double c6=0.;
  double m8=0.;
  double c8=0.;
  double set=0.;

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
  double sumhed4 = 0.; // Sum of high-eff discriminant for 4 highest HET jets
  double sumhpd4 = 0.; // Sum of high-pur discriminant for 4 highest HPT jets
  double sumhed6 = 0.; // Sum of high-eff discriminant for 6 highest HET jets
  double sumhpd6 = 0.; // Sum of high-pur discriminant for 6 highest HPT jets

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

      if ( cal->et() >= 40. && fabs( cal->eta() ) < 2.0 ) {
        goodIc5Jets=goodIc5Jets+1;
	if ( cal->discriminatorHighEff() > loose_ ) NTagsL++;
	if ( cal->discriminatorHighEff() > medium_ ) NTagsM++;
	if ( cal->discriminatorHighEff() > tight_ ) NTagsT++;
	GoodHt0+=cal->et();
	GoodHt20+=cal->et();
	GoodSumEt+=cal->et();
	double dp = 3.1415926-fabs(fabs(cal->phi()-metphi)-3.1415926);
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
    // Compute delta phi angles
    // ------------------------
    if ( iJ>1 ) dp12 = 3.1415926-fabs(fabs(Jphi[0]-Jphi[1])-3.1415926);
    if ( iJ>0 ) dp1st= 3.1415926-fabs(fabs(Jphi[0]-metphi)-3.1415926);
    if ( iJ>1 ) dp2nd= 3.1415926-fabs(fabs(Jphi[1]-metphi)-3.1415926);
    if ( iJ>2 ) dp3rd= 3.1415926-fabs(fabs(Jphi[2]-metphi)-3.1415926);

    double hemax[100]={0.};
    double hpmax[100]={0.};
    
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
	m3best=-999;  // Massa del tripletto piu' vicino a 172 GeV
	mwbest=-999;  // Massa del doppietto piu' vicino a 80.4 GeV 
	              // tra i tre nel tripletto definito sopra
	chi2=1000;    // Chiquadro costruito con le due masse qui sopra
	m45best=-999; // mass of two b-tags closest to 120 GeV not in triplet
	chi2ext=1000; // chi2 computed by adding contribution from distance to 120 GeV
	dpbb = 0;     // delta phi bb of pair above
	m45bestall=-999;  // Same as above, without bothering to exclude jets in triplet
	chi2extall=1000;  // Same as above
	dpbball = 0;      // delta phi bb as above
	m_others = -999;  // mass of jets not part of best triplet and pair
	mbbnohmax = -999; // mass of b pair with highest mass excluding b-jets from h decay
	dpbbnohmax = 0;   // angle of b pair with highest mass excluding b-jets from h decay
	
	m6=0.;
	m8=0.;
	set=0.;
	c8=0.;

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
	    // Extract random HED and HPD values from QCD pdf
	    // ----------------------------------------------
	    double x=0.;
	    for ( ; x<medium_; ) {
	      x=HEDpdf->GetRandom();
	    }   
	    JHET[j]=x;
	    x=0.;
	    for ( ; x<tight_; ) {
	      x=HPDpdf->GetRandom();
	    }   
	    JHPT[j]=x;
	  } else {
	    JHEM[j]=false;
	    PTOT=PTOT*(1.-PHETM[JBin[j]]);
	    PTOTE2+=pow(PHETMS[JBin[j]]/(1.-PHETM[JBin[j]]),2);
	    // Extract random HED and HPD values from QCD pdf
	    // ----------------------------------------------
	    double x=29.9;
	    for ( ; x>=medium_; ) {
	      x=HEDpdf->GetRandom(); 
	    }   
	    JHET[j]=x;
	    x=29.9;
	    for ( ; x>=tight_; ) {
	      x=HPDpdf->GetRandom();
	    }   
	    JHPT[j]=x;
	  }
	}
	PTOTE2=PTOTE2*pow(PTOT,2);

	sumhed4=0.;
	sumhpd4=0.;
	sumhed6=0.;
	sumhpd6=0.;
	// Compute sum of discriminants
	// ----------------------------
	for ( int i=0; i<NHSJ; i++ ) {
	  hemax[i]=JHET[i];
	  hpmax[i]=JHPT[i];
	}
	for ( int times=0; times<NHSJ; times++ ) {
	  for ( int j=NHSJ-1; j>=1; j-- ) {
	    if ( hemax[j]>hemax[j-1] ) {
	      double a1=hemax[j];
	      hemax[j]=hemax[j-1];
	      hemax[j-1]=a1;
	    }
	    if ( hpmax[j]>hpmax[j-1] ) {
	      double a1=hpmax[j];
	      hpmax[j]=hpmax[j-1];
	      hpmax[j-1]=a1;
	    }
	  }
	}
	for ( int i=0; i<4 && i<NHSJ; i++ ) {
	  sumhed4+=hemax[i];
	  sumhpd4+=hpmax[i];
	}
	for ( int i=0; i<6 && i<NHSJ; i++ ) {
	  sumhed6+=hemax[i];
	  sumhpd6+=hpmax[i];
	}

	// if ( NHEM>=2 ) {
          
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
		  if ( fabs(m45-120)<fabs(m45bestall-120) ) {
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
	  // Now compute mass of first 6 jets and their centrality
	  // -----------------------------------------------------
	  spx=0;
	  spy=0;
	  spz=0;
	  se=0;
	  set=0;
	  for ( int k=0; k<iJ && k<6; k++ ) {
	    px=JEt[k]*cos(Jphi[k]);
	    py=JEt[k]*sin(Jphi[k]);
	    pz=JEt[k]/tan(2*atan(exp(-Jeta[k])));
	    spx+=px;
	    spy+=py;
	    spz+=pz;
	    se+=sqrt(px*px+py*py+pz*pz);
	    set+=JEt[k];
	  }
	  m6=se*se-spx*spx-spy*spy-spz*spz;
	  if ( m6>0 ) { 
	    m6=sqrt(m6);
	    c6=set/m6;
	  }
          
	  ////////////////////////////////////////////////////////////////
	  // Now fill histograms
	  // -------------------
	  // Offline cuts
	  // ------------
	  if ( iJ >= 5 ) NJetsCut = true;
	  if ( metsig>3 ) MEtSigCut = true;
	  
	  // Fill kinematical histograms
	  // ---------------------------
	  double UncorrMEt =sqrt(UncorrSumEx*UncorrSumEx+UncorrSumEy*UncorrSumEy);
	  double CorrMEt = sqrt(CorrSumEx*CorrSumEx+CorrSumEy*CorrSumEy);
	  UncorrHt=UncorrHt0+UncorrMEt;
	  CorrHt=CorrHt0+CorrMEt;
	  GoodHt=GoodHt0+CorrMEt;
	  GoodHt2=GoodHt20+met;
	  
	  double UncorrMEtSig=0;
	  if ( UncorrSumEt>0 ) UncorrMEtSig = UncorrMEt/sqrt(UncorrSumEt);
	  double CorrMEtSig = 0;
	  if ( CorrSumEt>0 ) CorrMEtSig = CorrMEt/sqrt(CorrSumEt);  
	  
	  // Missing Et significance computed with tuned resolution function
	  // ---------------------------------------------------------------
	  double metsignew=0;
	  if ( sumet>0 ) 
	    metsignew = met/(sqrt(2)*0.8033*pow(sumet,0.5004));
	  
	  if ( PTOT>0 ) {
	    double weight_N = 1./pow (2.,iJmax);

	    // Fill regular histograms
	    // -----------------------
	    UncorrMEt_SumEt_->Fill(sumet,UncorrMEt,PTOT);
	    CorrMEt_SumEt_->Fill(sumet,CorrMEt,PTOT);
	    MEt_SumEt_->Fill(sumet,met,PTOT);
	    UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt,PTOT);
	    CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt,PTOT);
	    MEt_SumEtC_->Fill(CorrSumEt,met,PTOT);
	    UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt,PTOT);
	    CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt,PTOT);
	    MEt_SumEtJ_->Fill(GoodSumEt,met,PTOT);
	    
	    // Apply trigger requirement
	    // -------------------------
	    if ( response && NJetsCut ) {
	      NJets_->Fill(goodIc5Jets,PTOT);
	      UncorrHt_->Fill(UncorrHt,PTOT);
	      CorrHt_->Fill(CorrHt,PTOT);
	      GoodHt_->Fill(GoodHt,PTOT);
	      GoodHt2_->Fill(GoodHt2,PTOT);
	      UncorrSumEt_->Fill(UncorrSumEt,PTOT);
	      CorrSumEt_->Fill(CorrSumEt,PTOT);
	      GoodSumEt_->Fill(GoodSumEt,PTOT);
	      MEt_->Fill(met,PTOT);
	      MEtSig_->Fill(metsig,PTOT);
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
	      MEx_SumEt_->Fill(sumet,met*cos(metphi),PTOT);
	      MEx_SumEt_->Fill(sumet,met*sin(metphi),PTOT);
	      DP12_->Fill(dp12,PTOT);
	      DPbb_->Fill(dpbb,PTOT);
	      M_others_->Fill(m_others,PTOT);
	      Mbbnoh_->Fill(mbbnohmax,PTOT);
	      DPbbnoh_->Fill(dpbbnohmax,PTOT);
	      M6_->Fill(m6,PTOT);
	      C6_->Fill(c6,PTOT);
	      M8_->Fill(m8,PTOT);
	      C8_->Fill(c8,PTOT);
	      M45bestall_->Fill(m45bestall,PTOT);
	      Chi2extall_->Fill(chi2extall,PTOT);
	      DPbball_->Fill(dpbball,PTOT);
	      SumHED4_->Fill(sumhed4,PTOT);
	      SumHPD4_->Fill(sumhpd4,PTOT);
	      SumHED6_->Fill(sumhed6,PTOT);
	      SumHPD6_->Fill(sumhpd6,PTOT);
	      for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		HED_->Fill(JHET[i],PTOT);
		HPD_->Fill(JHPT[i],PTOT);
	      } 

	      NJetsN_->Fill(goodIc5Jets,weight_N);
	      UncorrHtN_->Fill(UncorrHt,weight_N);
	      CorrHtN_->Fill(CorrHt,weight_N);
	      GoodHtN_->Fill(GoodHt,weight_N);
	      GoodHt2N_->Fill(GoodHt2,weight_N);
	      UncorrSumEtN_->Fill(UncorrSumEt,weight_N);
	      CorrSumEtN_->Fill(CorrSumEt,weight_N);
	      GoodSumEtN_->Fill(GoodSumEt,weight_N);
	      MEtN_->Fill(met,weight_N);
	      MEtSigN_->Fill(metsig,weight_N);
	      MEtSigNewN_->Fill(metsignew,weight_N);
	      MEtDPMN_->Fill(dpmin,weight_N);
	      MEtDP1N_->Fill(dp1st,weight_N);
	      MEtDP2N_->Fill(dp2nd,weight_N);
	      MEtDP3N_->Fill(dp3rd,weight_N);
	      UncorrMEtSigN_->Fill(UncorrMEtSig,weight_N);
	      CorrMEtSigN_->Fill(CorrMEtSig,weight_N);
	      M3bestN_->Fill(m3best,weight_N);
	      MwbestN_->Fill(mwbest,weight_N);
	      Chi2massN_->Fill(chi2,weight_N);
	      M45bestN_->Fill(m45best,weight_N);
	      Chi2extN_->Fill(chi2ext,weight_N);
	      MEx_SumEtN_->Fill(sumet,met*cos(metphi),weight_N);
	      MEx_SumEtN_->Fill(sumet,met*sin(metphi),weight_N);
	      DP12N_->Fill(dp12,weight_N);
	      DPbbN_->Fill(dpbb,weight_N);
	      M_othersN_->Fill(m_others,weight_N);
	      MbbnohN_->Fill(mbbnohmax,weight_N);
	      DPbbnohN_->Fill(dpbbnohmax,weight_N);
	      M6N_->Fill(m6,weight_N);
	      C6N_->Fill(c6,weight_N);
	      M8N_->Fill(m8,weight_N);
	      C8N_->Fill(c8,weight_N);
	      M45bestallN_->Fill(m45bestall,weight_N);
	      Chi2extallN_->Fill(chi2extall,weight_N);
	      DPbballN_->Fill(dpbball,weight_N);
	      SumHED4N_->Fill(sumhed4,weight_N);
	      SumHPD4N_->Fill(sumhpd4,weight_N);
	      SumHED6N_->Fill(sumhed6,weight_N);
	      SumHPD6N_->Fill(sumhpd6,weight_N);
	      for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		HEDN_->Fill(JHET[i],weight_N);
		HPDN_->Fill(JHPT[i],weight_N);
	      } 
	      
	      if ( MEtSigCut && NHEM>=2 ) {
		NJetsS_->Fill(goodIc5Jets,PTOT);
		UncorrHtS_->Fill(UncorrHt,PTOT);
		CorrHtS_->Fill(CorrHt,PTOT);
		GoodHtS_->Fill(GoodHt,PTOT);
		GoodHt2S_->Fill(GoodHt2,PTOT);
		UncorrSumEtS_->Fill(UncorrSumEt,PTOT);
		CorrSumEtS_->Fill(CorrSumEt,PTOT);
		GoodSumEtS_->Fill(GoodSumEt,PTOT);
		MEtS_->Fill(met,PTOT);
		MEtSigS_->Fill(metsig,PTOT);
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
		MEx_SumEtS_->Fill(sumet,met*cos(metphi),PTOT);
		MEx_SumEtS_->Fill(sumet,met*sin(metphi),PTOT);
		DP12S_->Fill(dp12,PTOT);
		DPbbS_->Fill(dpbb,PTOT);
		M_othersS_->Fill(m_others,PTOT);
		MbbnohS_->Fill(mbbnohmax,PTOT);
		DPbbnohS_->Fill(dpbbnohmax,PTOT);
		M6S_->Fill(m6,PTOT);
		C6S_->Fill(c6,PTOT);
		M8S_->Fill(m8,PTOT);
		C8S_->Fill(c8,PTOT);
		M45bestallS_->Fill(m45bestall,PTOT);
		Chi2extallS_->Fill(chi2extall,PTOT);
		DPbballS_->Fill(dpbball,PTOT);
		SumHED4S_->Fill(sumhed4,PTOT);
		SumHPD4S_->Fill(sumhpd4,PTOT);
		SumHED6S_->Fill(sumhed6,PTOT);
		SumHPD6S_->Fill(sumhpd6,PTOT);
		for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		  HEDS_->Fill(JHET[i],PTOT);
		  HPDS_->Fill(JHPT[i],PTOT);
		} 

		NJetsSN_->Fill(goodIc5Jets,weight_N);
		UncorrHtSN_->Fill(UncorrHt,weight_N);
		CorrHtSN_->Fill(CorrHt,weight_N);
		GoodHtSN_->Fill(GoodHt,weight_N);
		GoodHt2SN_->Fill(GoodHt2,weight_N);
		UncorrSumEtSN_->Fill(UncorrSumEt,weight_N);
		CorrSumEtSN_->Fill(CorrSumEt,weight_N);
		GoodSumEtSN_->Fill(GoodSumEt,weight_N);
		MEtSN_->Fill(met,weight_N);
		MEtSigSN_->Fill(metsig,weight_N);
		MEtSigNewSN_->Fill(metsignew,weight_N);
		MEtDPMSN_->Fill(dpmin,weight_N);
		MEtDP1SN_->Fill(dp1st,weight_N);
		MEtDP2SN_->Fill(dp2nd,weight_N);
		MEtDP3SN_->Fill(dp3rd,weight_N);
		UncorrMEtSigSN_->Fill(UncorrMEtSig,weight_N);
		CorrMEtSigSN_->Fill(CorrMEtSig,weight_N);
		M3bestSN_->Fill(m3best,weight_N);
		MwbestSN_->Fill(mwbest,weight_N);
		Chi2massSN_->Fill(chi2,weight_N);
		M45bestSN_->Fill(m45best,weight_N);
		Chi2extSN_->Fill(chi2ext,weight_N);
		MEx_SumEtSN_->Fill(sumet,met*cos(metphi),weight_N);
		MEx_SumEtSN_->Fill(sumet,met*sin(metphi),weight_N);
		DP12SN_->Fill(dp12,weight_N);
		DPbbSN_->Fill(dpbb,weight_N);
		M_othersSN_->Fill(m_others,weight_N);
		MbbnohSN_->Fill(mbbnohmax,weight_N);
		DPbbnohSN_->Fill(dpbbnohmax,weight_N);
		M6SN_->Fill(m6,weight_N);
		C6SN_->Fill(c6,weight_N);
		M8SN_->Fill(m8,weight_N);
		C8SN_->Fill(c8,weight_N);
		M45bestallSN_->Fill(m45bestall,weight_N);
		Chi2extallSN_->Fill(chi2extall,weight_N);
		DPbballSN_->Fill(dpbball,weight_N);
		SumHED4SN_->Fill(sumhed4,weight_N);
		SumHPD4SN_->Fill(sumhpd4,weight_N);
		SumHED6SN_->Fill(sumhed6,weight_N);
		SumHPD6SN_->Fill(sumhpd6,weight_N);
		for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		  HEDSN_->Fill(JHET[i],weight_N);
		  HPDSN_->Fill(JHPT[i],weight_N);
		} 
	      
		if ( NHEM>=3 ) {
		  NJetsSS_->Fill(goodIc5Jets,PTOT);
		  UncorrHtSS_->Fill(UncorrHt,PTOT);
		  CorrHtSS_->Fill(CorrHt,PTOT);
		  GoodHtSS_->Fill(GoodHt,PTOT);
		  GoodHt2SS_->Fill(GoodHt2,PTOT);
		  UncorrSumEtSS_->Fill(UncorrSumEt,PTOT);
		  CorrSumEtSS_->Fill(CorrSumEt,PTOT);
		  GoodSumEtSS_->Fill(GoodSumEt,PTOT);
		  MEtSS_->Fill(met,PTOT);
		  MEtSigSS_->Fill(metsig,PTOT);
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
		  MEx_SumEtSS_->Fill(sumet,met*cos(metphi),PTOT);
		  MEx_SumEtSS_->Fill(sumet,met*sin(metphi),PTOT);
		  DP12SS_->Fill(dp12,PTOT);
		  DPbbSS_->Fill(dpbb,PTOT);
		  M_othersSS_->Fill(m_others,PTOT);      
		  MbbnohSS_->Fill(mbbnohmax,PTOT);
		  DPbbnohSS_->Fill(dpbbnohmax,PTOT);
		  M6SS_->Fill(m6,PTOT);
		  C6SS_->Fill(c6,PTOT);
		  M8SS_->Fill(m8,PTOT);
		  C8SS_->Fill(c8,PTOT);
		  M45bestallSS_->Fill(m45bestall,PTOT);
		  Chi2extallSS_->Fill(chi2extall,PTOT);
		  DPbballSS_->Fill(dpbball,PTOT);
		  SumHED4SS_->Fill(sumhed4,PTOT);
		  SumHPD4SS_->Fill(sumhpd4,PTOT);
		  SumHED6SS_->Fill(sumhed6,PTOT);
		  SumHPD6SS_->Fill(sumhpd6,PTOT);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSS_->Fill(JHET[i],PTOT);
		    HPDSS_->Fill(JHPT[i],PTOT);
		  } 

		  NJetsSSN_->Fill(goodIc5Jets,weight_N);
		  UncorrHtSSN_->Fill(UncorrHt,weight_N);
		  CorrHtSSN_->Fill(CorrHt,weight_N);
		  GoodHtSSN_->Fill(GoodHt,weight_N);
		  GoodHt2SSN_->Fill(GoodHt2,weight_N);
		  UncorrSumEtSSN_->Fill(UncorrSumEt,weight_N);
		  CorrSumEtSSN_->Fill(CorrSumEt,weight_N);
		  GoodSumEtSSN_->Fill(GoodSumEt,weight_N);
		  MEtSSN_->Fill(met,weight_N);
		  MEtSigSSN_->Fill(metsig,weight_N);
		  MEtSigNewSSN_->Fill(metsignew,weight_N);
		  MEtDPMSSN_->Fill(dpmin,weight_N);
		  MEtDP1SSN_->Fill(dp1st,weight_N);
		  MEtDP2SSN_->Fill(dp2nd,weight_N);
		  MEtDP3SSN_->Fill(dp3rd,weight_N);
		  UncorrMEtSigSSN_->Fill(UncorrMEtSig,weight_N);
		  CorrMEtSigSSN_->Fill(CorrMEtSig,weight_N);
		  M3bestSSN_->Fill(m3best,weight_N);
		  MwbestSSN_->Fill(mwbest,weight_N);
		  Chi2massSSN_->Fill(chi2,weight_N);
		  M45bestSSN_->Fill(m45best,weight_N);
		  Chi2extSSN_->Fill(chi2ext,weight_N);
		  MEx_SumEtSSN_->Fill(sumet,met*cos(metphi),weight_N);
		  MEx_SumEtSSN_->Fill(sumet,met*sin(metphi),weight_N);
		  DP12SSN_->Fill(dp12,weight_N);
		  DPbbSSN_->Fill(dpbb,weight_N);
		  M_othersSSN_->Fill(m_others,weight_N);      
		  MbbnohSSN_->Fill(mbbnohmax,weight_N);
		  DPbbnohSSN_->Fill(dpbbnohmax,weight_N);
		  M6SSN_->Fill(m6,weight_N);
		  C6SSN_->Fill(c6,weight_N);
		  M8SSN_->Fill(m8,weight_N);
		  C8SSN_->Fill(c8,weight_N);
		  M45bestallSSN_->Fill(m45bestall,weight_N);
		  Chi2extallSSN_->Fill(chi2extall,weight_N);
		  DPbballSSN_->Fill(dpbball,weight_N);
		  SumHED4SSN_->Fill(sumhed4,weight_N);
		  SumHPD4SSN_->Fill(sumhpd4,weight_N);
		  SumHED6SSN_->Fill(sumhed6,weight_N);
		  SumHPD6SSN_->Fill(sumhpd6,weight_N);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSSN_->Fill(JHET[i],weight_N);
		    HPDSSN_->Fill(JHPT[i],weight_N);
		  } 
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
		  MEtSSS_->Fill(met,PTOT);
		  MEtSigSSS_->Fill(metsig,PTOT);
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
		  MEx_SumEtSSS_->Fill(sumet,met*cos(metphi),PTOT);
		  MEx_SumEtSSS_->Fill(sumet,met*sin(metphi),PTOT);
		  DP12SSS_->Fill(dp12,PTOT);
		  DPbbSSS_->Fill(dpbb,PTOT);
		  M_othersSSS_->Fill(m_others,PTOT);      
		  MbbnohSSS_->Fill(mbbnohmax,PTOT);
		  DPbbnohSSS_->Fill(dpbbnohmax,PTOT);
		  M6SSS_->Fill(m6,PTOT);
		  C6SSS_->Fill(c6,PTOT);
		  M8SSS_->Fill(m8,PTOT);
		  C8SSS_->Fill(c8,PTOT);
		  M45bestallSSS_->Fill(m45bestall,PTOT);
		  Chi2extallSSS_->Fill(chi2extall,PTOT);
		  DPbballSSS_->Fill(dpbball,PTOT);
		  SumHED4SSS_->Fill(sumhed4,PTOT);
		  SumHPD4SSS_->Fill(sumhpd4,PTOT);
		  SumHED6SSS_->Fill(sumhed6,PTOT);
		  SumHPD6SSS_->Fill(sumhpd6,PTOT);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSSS_->Fill(JHET[i],PTOT);
		    HPDSSS_->Fill(JHPT[i],PTOT);
		  } 

		  NJetsSSSN_->Fill(goodIc5Jets,weight_N);
		  UncorrHtSSSN_->Fill(UncorrHt,weight_N);
		  CorrHtSSSN_->Fill(CorrHt,weight_N);
		  GoodHtSSSN_->Fill(GoodHt,weight_N);
		  GoodHt2SSSN_->Fill(GoodHt2,weight_N);
		  UncorrSumEtSSSN_->Fill(UncorrSumEt,weight_N);
		  CorrSumEtSSSN_->Fill(CorrSumEt,weight_N);
		  GoodSumEtSSSN_->Fill(GoodSumEt,weight_N);
		  MEtSSSN_->Fill(met,weight_N);
		  MEtSigSSSN_->Fill(metsig,weight_N);
		  MEtSigNewSSSN_->Fill(metsignew,weight_N);
		  MEtDPMSSSN_->Fill(dpmin,weight_N);
		  MEtDP1SSSN_->Fill(dp1st,weight_N);
		  MEtDP2SSSN_->Fill(dp2nd,weight_N);
		  MEtDP3SSSN_->Fill(dp3rd,weight_N);
		  UncorrMEtSigSSSN_->Fill(UncorrMEtSig,weight_N);
		  CorrMEtSigSSSN_->Fill(CorrMEtSig,weight_N);
		  M3bestSSSN_->Fill(m3best,weight_N);
		  MwbestSSSN_->Fill(mwbest,weight_N);
		  Chi2massSSSN_->Fill(chi2,weight_N);
		  M45bestSSSN_->Fill(m45best,weight_N);
		  Chi2extSSSN_->Fill(chi2ext,weight_N);
		  MEx_SumEtSSSN_->Fill(sumet,met*cos(metphi),weight_N);
		  MEx_SumEtSSSN_->Fill(sumet,met*sin(metphi),weight_N);
		  DP12SSSN_->Fill(dp12,weight_N);
		  DPbbSSSN_->Fill(dpbb,weight_N);
		  M_othersSSSN_->Fill(m_others,weight_N);      
		  MbbnohSSSN_->Fill(mbbnohmax,weight_N);
		  DPbbnohSSSN_->Fill(dpbbnohmax,weight_N);
		  M6SSSN_->Fill(m6,weight_N);
		  C6SSSN_->Fill(c6,weight_N);
		  M8SSSN_->Fill(m8,weight_N);
		  C8SSSN_->Fill(c8,weight_N);
		  M45bestallSSSN_->Fill(m45bestall,weight_N);
		  Chi2extallSSSN_->Fill(chi2extall,weight_N);
		  DPbballSSSN_->Fill(dpbball,weight_N);
		  SumHED4SSSN_->Fill(sumhed4,weight_N);
		  SumHPD4SSSN_->Fill(sumhpd4,weight_N);
		  SumHED6SSSN_->Fill(sumhed6,weight_N);
		  SumHPD6SSSN_->Fill(sumhpd6,weight_N);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSSSN_->Fill(JHET[i],weight_N);
		    HPDSSSN_->Fill(JHPT[i],weight_N);
		  } 
		}

	      }
	    } 

	    // Fill error histograms now
	    // -------------------------

	    UncorrMEt_SumEt_->Fill(sumet,UncorrMEt,PTOTE2);
	    CorrMEt_SumEt_->Fill(sumet,CorrMEt,PTOTE2);
	    MEt_SumEt_->Fill(sumet,met,PTOTE2);
	    UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt,PTOTE2);
	    CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt,PTOTE2);
	    MEt_SumEtC_->Fill(CorrSumEt,met,PTOTE2);
	    UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt,PTOTE2);
	    CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt,PTOTE2);
	    MEt_SumEtJ_->Fill(GoodSumEt,met,PTOTE2);
	    
	    // Apply trigger requirement
	    // -------------------------
	    if ( response && NJetsCut ) {
	      NJetsW_->Fill(goodIc5Jets,PTOTE2);
	      UncorrHtW_->Fill(UncorrHt,PTOTE2);
	      CorrHtW_->Fill(CorrHt,PTOTE2);
	      GoodHtW_->Fill(GoodHt,PTOTE2);
	      GoodHt2W_->Fill(GoodHt2,PTOTE2);
	      UncorrSumEtW_->Fill(UncorrSumEt,PTOTE2);
	      CorrSumEtW_->Fill(CorrSumEt,PTOTE2);
	      GoodSumEtW_->Fill(GoodSumEt,PTOTE2);
	      MEtW_->Fill(met,PTOTE2);
	      MEtSigW_->Fill(metsig,PTOTE2);
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
	      MEx_SumEtW_->Fill(sumet,met*cos(metphi),PTOTE2);
	      MEx_SumEtW_->Fill(sumet,met*sin(metphi),PTOTE2);
	      DP12W_->Fill(dp12,PTOTE2);
	      DPbbW_->Fill(dpbb,PTOTE2);
	      M_othersW_->Fill(m_others,PTOTE2);
	      MbbnohW_->Fill(mbbnohmax,PTOTE2);
	      DPbbnohW_->Fill(dpbbnohmax,PTOTE2);
	      M6W_->Fill(m6,PTOTE2);
	      C6W_->Fill(c6,PTOTE2);
	      M8W_->Fill(m8,PTOTE2);
	      C8W_->Fill(c8,PTOTE2);
	      M45bestallW_->Fill(m45bestall,PTOTE2);
	      Chi2extallW_->Fill(chi2extall,PTOTE2);
	      DPbballW_->Fill(dpbball,PTOTE2);
	      SumHED4W_->Fill(sumhed4,PTOTE2);
	      SumHPD4W_->Fill(sumhpd4,PTOTE2);
	      SumHED6W_->Fill(sumhed6,PTOTE2);
	      SumHPD6W_->Fill(sumhpd6,PTOTE2);
	      for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		HEDW_->Fill(JHET[i],PTOTE2);
		HPDW_->Fill(JHPT[i],PTOTE2);
	      } 
	      
	      if ( MEtSigCut && NHEM>=2 ) {
		NJetsSW_->Fill(goodIc5Jets,PTOTE2);
		UncorrHtSW_->Fill(UncorrHt,PTOTE2);
		CorrHtSW_->Fill(CorrHt,PTOTE2);
		GoodHtSW_->Fill(GoodHt,PTOTE2);
		GoodHt2SW_->Fill(GoodHt2,PTOTE2);
		UncorrSumEtSW_->Fill(UncorrSumEt,PTOTE2);
		CorrSumEtSW_->Fill(CorrSumEt,PTOTE2);
		GoodSumEtSW_->Fill(GoodSumEt,PTOTE2);
		MEtSW_->Fill(met,PTOTE2);
		MEtSigSW_->Fill(metsig,PTOTE2);
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
		MEx_SumEtSW_->Fill(sumet,met*cos(metphi),PTOTE2);
		MEx_SumEtSW_->Fill(sumet,met*sin(metphi),PTOTE2);
		DP12SW_->Fill(dp12,PTOTE2);
		DPbbSW_->Fill(dpbb,PTOTE2);
		M_othersSW_->Fill(m_others,PTOTE2);
		MbbnohSW_->Fill(mbbnohmax,PTOTE2);
		DPbbnohSW_->Fill(dpbbnohmax,PTOTE2);
		M6SW_->Fill(m6,PTOTE2);
		C6SW_->Fill(c6,PTOTE2);
		M8SW_->Fill(m8,PTOTE2);
		C8SW_->Fill(c8,PTOTE2);
		M45bestallSW_->Fill(m45bestall,PTOTE2);
		Chi2extallSW_->Fill(chi2extall,PTOTE2);
		DPbballSW_->Fill(dpbball,PTOTE2);
		SumHED4SW_->Fill(sumhed4,PTOTE2);
		SumHPD4SW_->Fill(sumhpd4,PTOTE2);
		SumHED6SW_->Fill(sumhed6,PTOTE2);
		SumHPD6SW_->Fill(sumhpd6,PTOTE2);
		for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		  HEDSW_->Fill(JHET[i],PTOTE2);
		  HPDSW_->Fill(JHPT[i],PTOTE2);
		} 
	      
		if ( NHEM>=3 ) {
		  NJetsSSW_->Fill(goodIc5Jets,PTOTE2);
		  UncorrHtSSW_->Fill(UncorrHt,PTOTE2);
		  CorrHtSSW_->Fill(CorrHt,PTOTE2);
		  GoodHtSSW_->Fill(GoodHt,PTOTE2);
		  GoodHt2SSW_->Fill(GoodHt2,PTOTE2);
		  UncorrSumEtSSW_->Fill(UncorrSumEt,PTOTE2);
		  CorrSumEtSSW_->Fill(CorrSumEt,PTOTE2);
		  GoodSumEtSSW_->Fill(GoodSumEt,PTOTE2);
		  MEtSSW_->Fill(met,PTOTE2);
		  MEtSigSSW_->Fill(metsig,PTOTE2);
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
		  MEx_SumEtSSW_->Fill(sumet,met*cos(metphi),PTOTE2);
		  MEx_SumEtSSW_->Fill(sumet,met*sin(metphi),PTOTE2);
		  DP12SSW_->Fill(dp12,PTOTE2);
		  DPbbSSW_->Fill(dpbb,PTOTE2);
		  M_othersSSW_->Fill(m_others,PTOTE2);      
		  MbbnohSSW_->Fill(mbbnohmax,PTOTE2);
		  DPbbnohSSW_->Fill(dpbbnohmax,PTOTE2);
		  M6SSW_->Fill(m6,PTOTE2);
		  C6SSW_->Fill(c6,PTOTE2);
		  M8SSW_->Fill(m8,PTOTE2);
		  C8SSW_->Fill(c8,PTOTE2);
		  M45bestallSSW_->Fill(m45bestall,PTOTE2);
		  Chi2extallSSW_->Fill(chi2extall,PTOTE2);
		  DPbballSSW_->Fill(dpbball,PTOTE2);
		  SumHED4SSW_->Fill(sumhed4,PTOTE2);
		  SumHPD4SSW_->Fill(sumhpd4,PTOTE2);
		  SumHED6SSW_->Fill(sumhed6,PTOTE2);
		  SumHPD6SSW_->Fill(sumhpd6,PTOTE2);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSSW_->Fill(JHET[i],PTOTE2);
		    HPDSSW_->Fill(JHPT[i],PTOTE2);
		  } 
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
		  MEtSSSW_->Fill(met,PTOTE2);
		  MEtSigSSSW_->Fill(metsig,PTOTE2);
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
		  MEx_SumEtSSSW_->Fill(sumet,met*cos(metphi),PTOTE2);
		  MEx_SumEtSSSW_->Fill(sumet,met*sin(metphi),PTOTE2);
		  DP12SSSW_->Fill(dp12,PTOTE2);
		  DPbbSSSW_->Fill(dpbb,PTOTE2);
		  M_othersSSSW_->Fill(m_others,PTOTE2);      
		  MbbnohSSSW_->Fill(mbbnohmax,PTOTE2);
		  DPbbnohSSSW_->Fill(dpbbnohmax,PTOTE2);
		  M6SSSW_->Fill(m6,PTOTE2);
		  C6SSSW_->Fill(c6,PTOTE2);
		  M8SSSW_->Fill(m8,PTOTE2);
		  C8SSSW_->Fill(c8,PTOTE2);
		  M45bestallSSSW_->Fill(m45bestall,PTOTE2);
		  Chi2extallSSSW_->Fill(chi2extall,PTOTE2);
		  DPbballSSSW_->Fill(dpbball,PTOTE2);
		  SumHED4SSSW_->Fill(sumhed4,PTOTE2);
		  SumHPD4SSSW_->Fill(sumhpd4,PTOTE2);
		  SumHED6SSSW_->Fill(sumhed6,PTOTE2);
		  SumHPD6SSSW_->Fill(sumhpd6,PTOTE2);
		  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
		    HEDSSSW_->Fill(JHET[i],PTOTE2);
		    HPDSSSW_->Fill(JHPT[i],PTOTE2);
		  } 
		}
	      }
	    }

	    // Compute Likelihood
	    // NB need to work on templates to make sure
	    // they are well-defined and not null over the
	    // support segment -> check Smooth.C
	    // --------------------------------------------
	    double r=1.0;
	    double fs=0.;
	    double fb=0.;
	    rel_lik = 0.;
	    int bx[11];
	    bx[0]= (int)((c8/1.2)*50)+1;
	    bx[1]= (int)((m45bestall/300.)*50)+1;
	    bx[2]= (int)(chi2extall)+1;
	    bx[3]= (int)((metsig/20.)*50)+1;
	    bx[4]= (int)((m8/2500.)*50.)+1;
	    bx[5]= (int)((met/500.)*50)+1;
	    bx[6]= (int)((dp12/3.2)*50.)+1;
	    bx[7]= (int)((dp2nd/3.2)*50)+1;
	    bx[8]= (int)((c6/1.2)*50)+1;
	    bx[9]= (int)((m6/2500.)*50.)+1;
	    bx[10]= (int)((sumet/4000.)*50)+1;
	    for ( int ivar=0; ivar<11; ivar++ ) {
	      if ( bx[ivar]>=1 && bx[ivar]<=50 ) {
		fs = HSS_sig[ivar]->GetBinContent(bx[ivar]);
		fb = HSS_bgr[ivar]->GetBinContent(bx[ivar]);
		if (fb>0 && fs>0 ) {
		  r *= fs/fb;
		} else {
		  if ( fb==0 ) r *= 10.;  // Need to improve this kludge
		  if ( fs==0 ) r *= 0.1; 
		}
	      }
	    }
	    if ( r>=0 ) {
	      rel_lik += log(r);
	    } else {
              rel_lik=-9.99;	    
	    }
	    if ( rel_lik<-10. ) rel_lik = -9.99;
	    if ( rel_lik>=10. ) rel_lik = 9.99;
	    if ( response && NJetsCut ) {
	      L_->Fill(rel_lik,PTOT);
	      LW_->Fill(rel_lik,PTOTE2);
	      LN_->Fill(rel_lik,weight_N);
	      if ( MEtSigCut && NHEM>=2 ) {
		LS_->Fill(rel_lik,PTOT);
		LSW_->Fill(rel_lik,PTOTE2);
		LSN_->Fill(rel_lik,weight_N);
		if ( NHEM>=3 ) {
		  LSS_->Fill(rel_lik,PTOT);
		  LSSW_->Fill(rel_lik,PTOTE2);
		  LSSN_->Fill(rel_lik,weight_N);
		}
		if ( NHEM>=4 ) {
		  LSSS_->Fill(rel_lik,PTOT);
		  LSSSW_->Fill(rel_lik,PTOTE2);
		  LSSSN_->Fill(rel_lik,weight_N);
		}
	      }
	    }

	  } // end if PTOT>0
	  //	}  // end if ntags>=2

      } // end for tag combinatorial
      
    } else { // if not QCD /////////////////////////////////////////////

      sumhed4=0.;
      sumhpd4=0.;
      sumhed6=0.;
      sumhpd6=0.;
      // Compute sum of discriminants
      // ----------------------------
      for ( int i=0; i<NHSJ; i++ ) {
	hemax[i]=JHET[i];
	hpmax[i]=JHPT[i];
      }
      for ( int times=0; times<NHSJ; times++ ) {
	for ( int j=NHSJ-1; j>=1; j-- ) {
	  if ( hemax[j]>hemax[j-1] ) {
	    double a1=hemax[j];
	    hemax[j]=hemax[j-1];
	    hemax[j-1]=a1;
	  }
	  if ( hpmax[j]>hpmax[j-1] ) {
	    double a1=hpmax[j];
	    hpmax[j]=hpmax[j-1];
	    hpmax[j-1]=a1;
	  }
	}
      }
      for ( int i=0; i<4 && i<NHSJ; i++ ) {
	sumhed4+=hemax[i];
	sumhpd4+=hpmax[i];
      }
      for ( int i=0; i<6 && i<NHSJ; i++ ) {
	sumhed6+=hemax[i];
	sumhpd6+=hpmax[i];
      }

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
	      if ( fabs(m45-120)<fabs(m45bestall-120) ) {
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
      // Now compute mass of first 6 jets and their centrality
      // -----------------------------------------------------
      spx=0;
      spy=0;
      spz=0;
      se=0;
      set=0;
      for ( int k=0; k<iJ && k<6; k++ ) {
	px=JEt[k]*cos(Jphi[k]);
	py=JEt[k]*sin(Jphi[k]);
	pz=JEt[k]/tan(2*atan(exp(-Jeta[k])));
	spx+=px;
	spy+=py;
	spz+=pz;
	se+=sqrt(px*px+py*py+pz*pz);
	set+=JEt[k];
      }
      m6=se*se-spx*spx-spy*spy-spz*spz;
      if ( m6>0 ) { 
	m6=sqrt(m6);
	c6=set/m6;
      }

      ////////////////////////////////////////////////////////////////
      // Now fill histograms
      // -------------------
      // Offline cuts
      // ------------
      if ( iJ >= 5 ) NJetsCut = true;
      if ( metsig>3 ) MEtSigCut = true;
      
      // Fill kinematical histograms
      // ---------------------------
      double UncorrMEt =sqrt(UncorrSumEx*UncorrSumEx+UncorrSumEy*UncorrSumEy);
      double CorrMEt = sqrt(CorrSumEx*CorrSumEx+CorrSumEy*CorrSumEy);
      UncorrHt=UncorrHt0+UncorrMEt;
      CorrHt=CorrHt0+CorrMEt;
      GoodHt=GoodHt0+CorrMEt;
      GoodHt2=GoodHt20+met;
      
      double UncorrMEtSig=0;
      if ( UncorrSumEt>0 ) UncorrMEtSig = UncorrMEt/sqrt(UncorrSumEt);
      double CorrMEtSig = 0;
      if ( CorrSumEt>0 ) CorrMEtSig = CorrMEt/sqrt(CorrSumEt);  
      
      // Missing Et significance computed with tuned resolution function
      // ---------------------------------------------------------------
      double metsignew=0;
      if ( sumet>0 ) metsignew = met/(sqrt(2)*0.8033*pow(sumet,0.5004));
      
      UncorrMEt_SumEt_->Fill(sumet,UncorrMEt);
      CorrMEt_SumEt_->Fill(sumet,CorrMEt);
      MEt_SumEt_->Fill(sumet,met);
      UncorrMEt_SumEtC_->Fill(CorrSumEt,UncorrMEt);
      CorrMEt_SumEtC_->Fill(CorrSumEt,CorrMEt);
      MEt_SumEtC_->Fill(CorrSumEt,met);
      UncorrMEt_SumEtJ_->Fill(GoodSumEt,UncorrMEt);
      CorrMEt_SumEtJ_->Fill(GoodSumEt,CorrMEt);
      MEt_SumEtJ_->Fill(GoodSumEt,met);
      
      // Apply trigger requirement
      // -------------------------
      if ( response && NJetsCut ) {
	NJets_->Fill(goodIc5Jets);
	UncorrHt_->Fill(UncorrHt);
	CorrHt_->Fill(CorrHt);
	GoodHt_->Fill(GoodHt);
	GoodHt2_->Fill(GoodHt2);
	UncorrSumEt_->Fill(UncorrSumEt);
	CorrSumEt_->Fill(CorrSumEt);
	GoodSumEt_->Fill(GoodSumEt);
	MEt_->Fill(met);
	MEtSig_->Fill(metsig);
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
	MEx_SumEt_->Fill(sumet,met*cos(metphi));
	MEx_SumEt_->Fill(sumet,met*sin(metphi));
	DP12_->Fill(dp12);
	DPbb_->Fill(dpbb);
	M_others_->Fill(m_others);
	Mbbnoh_->Fill(mbbnohmax);
	DPbbnoh_->Fill(dpbbnohmax);
	M6_->Fill(m6);
	C6_->Fill(c6);
	M8_->Fill(m8);
	C8_->Fill(c8);
	M45bestall_->Fill(m45bestall);
	Chi2extall_->Fill(chi2extall);
	DPbball_->Fill(dpbball);
	SumHED4_->Fill(sumhed4);
	SumHPD4_->Fill(sumhpd4);
	SumHED6_->Fill(sumhed6);
	SumHPD6_->Fill(sumhpd6);
	for ( int i=0; i<iJ && i<NHSJ; i++ ) {
	  HED_->Fill(JHET[i]);
	  HPD_->Fill(JHPT[i]);
	} 
	
	if ( MEtSigCut && NHEM>=2 ) {
	  NJetsS_->Fill(goodIc5Jets);
	  UncorrHtS_->Fill(UncorrHt);
	  CorrHtS_->Fill(CorrHt);
	  GoodHtS_->Fill(GoodHt);
	  GoodHt2S_->Fill(GoodHt2);
	  UncorrSumEtS_->Fill(UncorrSumEt);
	  CorrSumEtS_->Fill(CorrSumEt);
	  GoodSumEtS_->Fill(GoodSumEt);
	  MEtS_->Fill(met);
	  MEtSigS_->Fill(metsig);
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
	  MEx_SumEtS_->Fill(sumet,met*cos(metphi));
	  MEx_SumEtS_->Fill(sumet,met*sin(metphi));
	  DP12S_->Fill(dp12);
	  DPbbS_->Fill(dpbb);
	  M_othersS_->Fill(m_others);
	  MbbnohS_->Fill(mbbnohmax);
	  DPbbnohS_->Fill(dpbbnohmax);
	  M6S_->Fill(m6);
	  C6S_->Fill(c6);
	  M8S_->Fill(m8);
	  C8S_->Fill(c8);
	  M45bestallS_->Fill(m45bestall);
	  Chi2extallS_->Fill(chi2extall);
	  DPbballS_->Fill(dpbball);
	  SumHED4S_->Fill(sumhed4);
	  SumHPD4S_->Fill(sumhpd4);
	  SumHED6S_->Fill(sumhed6);
	  SumHPD6S_->Fill(sumhpd6);
	  for ( int i=0; i<iJ && i<NHSJ; i++ ) {
	    HEDS_->Fill(JHET[i]);
	    HPDS_->Fill(JHPT[i]);
	  } 
	  
	  if ( NHEM>=3 ) {
	    NJetsSS_->Fill(goodIc5Jets);
	    UncorrHtSS_->Fill(UncorrHt);
	    CorrHtSS_->Fill(CorrHt);
	    GoodHtSS_->Fill(GoodHt);
	    GoodHt2SS_->Fill(GoodHt2);
	    UncorrSumEtSS_->Fill(UncorrSumEt);
	    CorrSumEtSS_->Fill(CorrSumEt);
	    GoodSumEtSS_->Fill(GoodSumEt);
	    MEtSS_->Fill(met);
	    MEtSigSS_->Fill(metsig);
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
	    MEx_SumEtSS_->Fill(sumet,met*cos(metphi));
	    MEx_SumEtSS_->Fill(sumet,met*sin(metphi));
	    DP12SS_->Fill(dp12);
	    DPbbSS_->Fill(dpbb);
	    M_othersSS_->Fill(m_others);      
	    MbbnohSS_->Fill(mbbnohmax);
	    DPbbnohSS_->Fill(dpbbnohmax);
	    M6SS_->Fill(m6);
	    C6SS_->Fill(c6);
	    M8SS_->Fill(m8);
	    C8SS_->Fill(c8);
	    M45bestallSS_->Fill(m45bestall);
	    Chi2extallSS_->Fill(chi2extall);
	    DPbballSS_->Fill(dpbball);
	    SumHED4SS_->Fill(sumhed4);
	    SumHPD4SS_->Fill(sumhpd4);
	    SumHED6SS_->Fill(sumhed6);
	    SumHPD6SS_->Fill(sumhpd6);
	    for ( int i=0; i<iJ && i<NHSJ; i++ ) {
	      HEDSS_->Fill(JHET[i]);
	      HPDSS_->Fill(JHPT[i]);
	    } 
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
	    MEtSSS_->Fill(met);
	    MEtSigSSS_->Fill(metsig);
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
	    MEx_SumEtSSS_->Fill(sumet,met*cos(metphi));
	    MEx_SumEtSSS_->Fill(sumet,met*sin(metphi));
	    DP12SSS_->Fill(dp12);
	    DPbbSSS_->Fill(dpbb);
	    M_othersSSS_->Fill(m_others);      
	    MbbnohSSS_->Fill(mbbnohmax);
	    DPbbnohSSS_->Fill(dpbbnohmax);
	    M6SSS_->Fill(m6);
	    C6SSS_->Fill(c6);
	    M8SSS_->Fill(m8);
	    C8SSS_->Fill(c8);
	    M45bestallSSS_->Fill(m45bestall);
	    Chi2extallSSS_->Fill(chi2extall);
	    DPbballSSS_->Fill(dpbball);
	    SumHED4SSS_->Fill(sumhed4);
	    SumHPD4SSS_->Fill(sumhpd4);
	    SumHED6SSS_->Fill(sumhed6);
	    SumHPD6SSS_->Fill(sumhpd6);
	    for ( int i=0; i<iJ && i<NHSJ; i++ ) {
	      HEDSSS_->Fill(JHET[i]);
	      HPDSSS_->Fill(JHPT[i]);
	    } 

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

      // Compute Likelihood
      // ------------------
      double r=1.0;
      double fs=0.;
      double fb=0.;
      rel_lik = 0.;
      int bx[11];
      bx[0]= (int)((c8/1.2)*50)+1;
      bx[1]= (int)((m45bestall/300.)*50)+1;
      bx[2]= (int)(chi2extall)+1;
      bx[3]= (int)((metsig/20.)*50)+1;
      bx[4]= (int)((m8/2000)*50.)+1;
      bx[5]= (int)((met/500.)*50)+1;
      bx[6]= (int)((dp12/3.2)*50.)+1;
      bx[7]= (int)((dp2nd/3.2)*50)+1;
      bx[8]= (int)((c6/1.2)*50)+1;
      bx[9]= (int)((m6/2500.)*50.)+1;
      bx[10]= (int)((sumet/4000.)*50)+1;
      for ( int ivar=0; ivar<11; ivar++ ) {
	if ( bx[ivar]>=1 && bx[ivar]<=50 ) {
	  fs = HSS_sig[ivar]->GetBinContent(bx[ivar]);
	  fb = HSS_bgr[ivar]->GetBinContent(bx[ivar]);
	  if (fb>0 && fs>0 ) {
	    r *= fs/fb;
	  } else {
	    if ( fb==0 ) r *= 10.;  // Need to improve this kludge
	    if ( fs==0 ) r *= 0.1; 
	  }
	}
      }
      if ( r>=0 ) {
	rel_lik += log(r);
      } else {
	rel_lik=-9.99;	    
      }
      if ( rel_lik<-10. ) rel_lik = -9.99;
      if ( rel_lik>=10. ) rel_lik = 9.99;
      if ( response && NJetsCut ) {
	L_->Fill(rel_lik);
	if ( MEtSigCut && NHEM>=2 ) {
	  LS_->Fill(rel_lik);
	  if ( NHEM>=3 ) {
	    LSS_->Fill(rel_lik);
	  }
	  if ( NHEM>=4 ) {
	    LSSS_->Fill(rel_lik);
	  }
	}
      }

    }
  } // end if caloJETS.size()
  else {
    std::cout << "ATTENTION: Jet collection empty" << std::endl;
  }    


}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void TDAna::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void TDAna::endJob() {

  // Write stats on decays
  // ---------------------
  double f0[5]={0.};
  double sf0[5]={0.};
  double f0t=0.;
  double sf0t=0.;
  TString Hname[5] = { "  H->others ", "    H->bb   ", "    H->cc   ", 
		       " H->tau tau ", "    H->WW   " };
  TString Tname[10] = { "tt->jjbjjb  ", "tt->evbjjb  ", "tt->mvbjjb  ", "tt->tvbjjb  ", 
			"tt->evbenb  ", "tt->evbmvb  ", "tt->evbtvb  ", "tt->mnbmnb  ", 
			"tt->mvbtvb  ", "tt->tnbtnb  " };
  // File with HEP table
  // -------------------
  ofstream decayfile("tthdecays.asc");
  decayfile << endl;
  decayfile << "            " << Hname[1] << "  " << Hname[2] << "  " << Hname[3] << "  " 
	    << Hname[4] << "  " << Hname[0] << "   Total   " << endl;
  decayfile << "--------------------------------------------------" 
	    << "-------------------------------------" << endl;
  for ( int itdecay=0; itdecay<10; itdecay++ ) {
    cout << endl;
    cout << Tname[itdecay] << " : " << endl;
    cout << "--------------" << endl;
    if ( grandtotalh[itdecay]>0 ) {
      f0t = grandtotalhpass[itdecay]/grandtotalh[itdecay];
      sf0t = sqrt(f0t*(1-f0t)/grandtotalh[itdecay]);
    }
    cout << "Total       " << setprecision(6) << grandtotalhpass[itdecay] << "/" 
	 << grandtotalh[itdecay] << " = " 
	 << "(" <<setprecision(5) << f0t*100. << "+-" 
	 << setprecision(5) << sf0t*100. << ") %" << endl;
    for ( int hdecay=0; hdecay<5; hdecay++ ) {
      if ( total[itdecay][hdecay]>0 ) {
	f0[hdecay] = totalpass[itdecay][hdecay]/total[itdecay][hdecay];
	sf0[hdecay] = sqrt(f0[hdecay]*(1-f0[hdecay])/total[itdecay][hdecay]);
      }
      cout << Hname[hdecay] << setprecision(6) << totalpass[itdecay][hdecay] << "/" 
	   << total[itdecay][hdecay] << " = " << "(" << setprecision(5) << f0[hdecay]*100. 
	   << "+-" << setprecision(5) << sf0[hdecay]*100. << ") %" << endl;
    } 
    
    decayfile << Tname[itdecay] << setprecision(6) << setw(7) << totalpass[itdecay][1] << "/" 
	      << setprecision(6) << total[itdecay][1] << " "  << setprecision(6) 
	      << setw(8) << totalpass[itdecay][2] << "/" 
	      << setprecision(6) << total[itdecay][2] << " "  << setprecision(6) 
	      << setw(8) << totalpass[itdecay][3] << "/" 
	      << setprecision(6) << total[itdecay][3] << " "  << setprecision(6) 
	      << setw(8) << totalpass[itdecay][4] << "/" 
	      << setprecision(6) << total[itdecay][4] << " "  << setprecision(6) 
	      << setw(8) << totalpass[itdecay][0] << "/" 
	      << setprecision(6) << total[itdecay][0] << " "  << setprecision(6) 
	      << setw(8) << grandtotalhpass[itdecay] << "/" 
	      << setprecision(6) << grandtotalh[itdecay] << endl;
    decayfile << "            "  
	      << setprecision(5) << setw(5) << f0[1]*100. << "+-" << setprecision(4) << sf0[1]*100. << " " 
	      << setprecision(5) << setw(5) << f0[2]*100. << "+-" << setprecision(4) << sf0[2]*100. << " " 
	      << setprecision(5) << setw(5) << f0[3]*100. << "+-" << setprecision(4) << sf0[3]*100. << " " 
	      << setprecision(5) << setw(5) << f0[4]*100. << "+-" << setprecision(4) << sf0[4]*100. << " " 
	      << setprecision(5) << setw(5) << f0[0]*100. << "+-" << setprecision(4) << sf0[0]*100. << " "
	      << setprecision(5) << setw(5) << f0t*100.   << "+-" << setprecision(4) << sf0t*100. << endl;
    decayfile << endl;
  }

  cout << endl;
  cout << "Totals : " << endl;
  cout << "-------- " << endl;
  f0t=0.;
  sf0t=0.;
  for ( int hdecay=0; hdecay<5; hdecay++ ) {
    if ( grandtotaltt[hdecay]>0 ) {
      f0[hdecay] = grandtotalttpass[hdecay]/grandtotaltt[hdecay];
      sf0[hdecay] = sqrt(f0[hdecay]*(1-f0[hdecay])/grandtotaltt[hdecay]);
    }
    cout << Hname[hdecay] << setprecision(6) << grandtotalttpass[hdecay] << "/" 
	 << grandtotaltt[hdecay] << " = " << "(" <<setprecision(5) << f0[hdecay]*100. 
	 << "+-" << setprecision(5) << sf0[hdecay]*100. << ") %" << endl;
  }
  cout << endl;
  cout << "Grandtotal : " << endl;
  cout << "------------ " << endl;
  f0t=0.;
  sf0t=0.;
  if ( grandgrandtotal>0 ) {
    f0t = grandgrandtotalpass/grandgrandtotal;
    sf0t = sqrt(f0t*(1-f0t)/grandgrandtotal);
  }
  cout << setprecision(6) << grandgrandtotalpass << "/" << grandgrandtotal << " = " 
       << "(" <<setprecision(5) << f0t*100. 
       << "+-" << setprecision(5) << sf0t*100. << ") %" << endl;

  decayfile << "Totals:     " 
	    << setprecision(6) << setw(7) << grandtotalttpass[1] << "/" 
	    <<  setprecision(5) << grandtotaltt[1] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[2] << "/" 
	    <<  setprecision(5) << grandtotaltt[2] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[3] << "/" 
	    <<  setprecision(5) << grandtotaltt[3] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[4] << "/" 
	    <<  setprecision(5) << grandtotaltt[4] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[0] << "/" 
	    <<  setprecision(5) << grandtotaltt[0] << " " 
	    << setprecision(6) << setw(8) << grandgrandtotalpass << "/" 
	    <<  setprecision(5) << grandgrandtotal << endl;
  decayfile << "            " 
	    << setprecision(5) << setw(5) << f0[1]*100. << "+-" 
	    << setprecision(4) << sf0[1]*100.  << " " 
	    << setprecision(5) << setw(5) << f0[2]*100.  << "+-" 
	    << setprecision(4) << sf0[2]*100.  << " " 
	    << setprecision(5) << setw(5) << f0[3]*100.  << "+-" 
	    << setprecision(4) << sf0[3]*100.  << " " 
	    << setprecision(5) << setw(5) << f0[4]*100.  << "+-" 
	    << setprecision(4) << sf0[4]*100.  << " " 
	    << setprecision(5) << setw(5) << f0[0]*100.  << "+-" 
	    << setprecision(4) << sf0[0]*100.  << " " 
	    << setprecision(5) << setw(5) << f0t*100.    << "+-" 
	    << setprecision(4) << sf0t*100.  << endl;

  decayfile.close();

  // Histograms for the study of jet/parton association
  // --------------------------------------------------
  Drmax_->Scale(1./Nsltt_hjj);
  Drmedall_->Scale(1./Nsltt_hjj);
  Drmed07_->Scale(1./Nsltt_hjj);
  N07_->Scale(1./Nsltt_hjj);
  N04_->Scale(1./Nsltt_hjj);
  N02_->Scale(1./Nsltt_hjj);
  Nlo_->Scale(1./Nsltt_hjj);
  Detmed07_->Scale(1./Nsltt_hjj);
  Detmedall_->Scale(1./Nsltt_hjj);
  Perf07_->Scale(1./Nsltt_hjj);
  Perf04_->Scale(1./Nsltt_hjj);
  Perf02_->Scale(1./Nsltt_hjj);
  Det2med07_->Scale(1./Nsltt_hjj);
  Det2medall_->Scale(1./Nsltt_hjj);
  Hrecfrac_->Scale(1./Nsltt_hjj);
  Drmax_->Write();
  Drmedall_->Write();
  Drmed07_->Write();
  N07_->Write();
  N04_->Write();
  N02_->Write();
  Nlo_->Write();
  Detmedall_->Write();
  Detmed07_->Write();
  Perf07_->Write();
  Perf04_->Write();
  Perf02_->Write();
  Det2med07_->Write();
  Det2medall_->Write();
  Hrecfrac_->Write();

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
  M6_->Write();
  C6_->Write();
  M8_->Write();
  C8_->Write();
  M45bestall_->Write();
  Chi2extall_->Write();
  DPbball_->Write();
  SumHED4_->Write();
  SumHPD4_->Write();
  SumHED6_->Write();
  SumHPD6_->Write();
  HED_->Write();
  HPD_->Write();

  NJetsN_->Write();
  UncorrHtN_->Write();
  CorrHtN_->Write();
  GoodHtN_->Write();
  GoodHt2N_->Write();
  UncorrSumEtN_->Write();
  CorrSumEtN_->Write();
  GoodSumEtN_->Write();
  MEtN_->Write();
  MEtSigN_->Write();
  MEtSigNewN_->Write();
  MEtDPMN_->Write();
  MEtDP1N_->Write();
  MEtDP2N_->Write();
  MEtDP3N_->Write();
  UncorrMEtSigN_->Write();
  CorrMEtSigN_->Write();
  M3bestN_->Write();
  MwbestN_->Write();
  Chi2massN_->Write();
  M45bestN_->Write();
  Chi2extN_->Write();
  MEx_SumEtN_->Write();
  DP12N_->Write();
  DPbbN_->Write();
  M_othersN_->Write();
  MbbnohN_->Write();
  DPbbnohN_->Write();
  M6N_->Write();
  C6N_->Write();
  M8N_->Write();
  C8N_->Write();
  M45bestallN_->Write();
  Chi2extallN_->Write();
  DPbballN_->Write();
  SumHED4N_->Write();
  SumHPD4N_->Write();
  SumHED6N_->Write();
  SumHPD6N_->Write();
  HEDN_->Write();
  HPDN_->Write();
  
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
  M6S_->Write();
  C6S_->Write();
  M8S_->Write();
  C8S_->Write();
  M45bestallS_->Write();
  Chi2extallS_->Write();
  DPbballS_->Write();
  SumHED4S_->Write();
  SumHPD4S_->Write();
  SumHED6S_->Write();
  SumHPD6S_->Write();
  HEDS_->Write();
  HPDS_->Write();

  NJetsSN_->Write();
  UncorrHtSN_->Write();
  CorrHtSN_->Write();
  GoodHtSN_->Write();
  GoodHt2SN_->Write();
  UncorrSumEtSN_->Write();
  CorrSumEtSN_->Write();
  GoodSumEtSN_->Write();
  MEtSN_->Write();
  MEtSigSN_->Write();
  MEtSigNewSN_->Write();
  MEtDPMSN_->Write();
  MEtDP1SN_->Write();
  MEtDP2SN_->Write();
  MEtDP3SN_->Write();
  UncorrMEtSigSN_->Write();
  CorrMEtSigSN_->Write();
  M3bestSN_->Write();
  MwbestSN_->Write();
  Chi2massSN_->Write();
  M45bestSN_->Write();
  Chi2extSN_->Write();
  MEx_SumEtSN_->Write();
  DP12SN_->Write();
  DPbbSN_->Write();
  M_othersSN_->Write();
  MbbnohSN_->Write();
  DPbbnohSN_->Write();
  M6SN_->Write();
  C6SN_->Write();
  M8SN_->Write();
  C8SN_->Write();
  M45bestallSN_->Write();
  Chi2extallSN_->Write();
  DPbballSN_->Write();
  SumHED4SN_->Write();
  SumHPD4SN_->Write();
  SumHED6SN_->Write();
  SumHPD6SN_->Write();
  HEDSN_->Write();
  HPDSN_->Write();

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
  M6SS_->Write();
  C6SS_->Write();
  M8SS_->Write();
  C8SS_->Write();
  M45bestallSS_->Write();
  Chi2extallSS_->Write();
  DPbballSS_->Write();
  SumHED4SS_->Write();
  SumHPD4SS_->Write();
  SumHED6SS_->Write();
  SumHPD6SS_->Write();
  HEDSS_->Write();
  HPDSS_->Write();

  NJetsSSN_->Write();
  UncorrHtSSN_->Write();
  CorrHtSSN_->Write();
  GoodHtSSN_->Write();
  GoodHt2SSN_->Write();
  UncorrSumEtSSN_->Write();
  CorrSumEtSSN_->Write();
  GoodSumEtSSN_->Write();
  MEtSSN_->Write();
  MEtSigSSN_->Write();
  MEtSigNewSSN_->Write();
  MEtDPMSSN_->Write();
  MEtDP1SSN_->Write();
  MEtDP2SSN_->Write();
  MEtDP3SSN_->Write();
  UncorrMEtSigSSN_->Write();
  CorrMEtSigSSN_->Write();
  M3bestSSN_->Write();
  MwbestSSN_->Write();
  Chi2massSSN_->Write();
  M45bestSSN_->Write();
  Chi2extSSN_->Write();
  MEx_SumEtSSN_->Write();
  DP12SSN_->Write();
  DPbbSSN_->Write();
  M_othersSSN_->Write();
  MbbnohSSN_->Write();
  DPbbnohSSN_->Write();
  M6SSN_->Write();
  C6SSN_->Write();
  M8SSN_->Write();
  C8SSN_->Write();
  M45bestallSSN_->Write();
  Chi2extallSSN_->Write();
  DPbballSSN_->Write();
  SumHED4SSN_->Write();
  SumHPD4SSN_->Write();
  SumHED6SSN_->Write();
  SumHPD6SSN_->Write();
  HEDSSN_->Write();
  HPDSSN_->Write();

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
  M6SSS_->Write();
  C6SSS_->Write();
  M8SSS_->Write();
  C8SSS_->Write();
  M45bestallSSS_->Write();
  Chi2extallSSS_->Write();
  DPbballSSS_->Write();
  SumHED4SSS_->Write();
  SumHPD4SSS_->Write();
  SumHED6SSS_->Write();
  SumHPD6SSS_->Write();
  HEDSSS_->Write();
  HPDSSS_->Write();

  NJetsSSSN_->Write();
  UncorrHtSSSN_->Write();
  CorrHtSSSN_->Write();
  GoodHtSSSN_->Write();
  GoodHt2SSSN_->Write();
  UncorrSumEtSSSN_->Write();
  CorrSumEtSSSN_->Write();
  GoodSumEtSSSN_->Write();
  MEtSSSN_->Write();
  MEtSigSSSN_->Write();
  MEtSigNewSSSN_->Write();
  MEtDPMSSSN_->Write();
  MEtDP1SSSN_->Write();
  MEtDP2SSSN_->Write();
  MEtDP3SSSN_->Write();
  UncorrMEtSigSSSN_->Write();
  CorrMEtSigSSSN_->Write();
  M3bestSSSN_->Write();
  MwbestSSSN_->Write();
  Chi2massSSSN_->Write();
  M45bestSSSN_->Write();
  Chi2extSSSN_->Write();
  MEx_SumEtSSSN_->Write();
  DP12SSSN_->Write();
  DPbbSSSN_->Write();
  M_othersSSSN_->Write();
  MbbnohSSSN_->Write();
  DPbbnohSSSN_->Write();
  M6SSSN_->Write();
  C6SSSN_->Write();
  M8SSSN_->Write();
  C8SSSN_->Write();
  M45bestallSSSN_->Write();
  Chi2extallSSSN_->Write();
  DPbballSSSN_->Write();
  SumHED4SSSN_->Write();
  SumHPD4SSSN_->Write();
  SumHED6SSSN_->Write();
  SumHPD6SSSN_->Write();
  HEDSSSN_->Write();
  HPDSSSN_->Write();

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
  M6W_->Write();
  C6W_->Write();
  M8W_->Write();
  C8W_->Write();
  M45bestallW_->Write();
  Chi2extallW_->Write();
  DPbballW_->Write();
  SumHED4W_->Write();
  SumHPD4W_->Write();
  SumHED6W_->Write();
  SumHPD6W_->Write();
  HEDW_->Write();
  HPDW_->Write();
  
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
  M6SW_->Write();
  C6SW_->Write();
  M8SW_->Write();
  C8SW_->Write();
  M45bestallSW_->Write();
  Chi2extallSW_->Write();
  DPbballSW_->Write();
  SumHED4SW_->Write();
  SumHPD4SW_->Write();
  SumHED6SW_->Write();
  SumHPD6SW_->Write();
  HEDSW_->Write();
  HPDSW_->Write();

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
  M6SSW_->Write();
  C6SSW_->Write();
  M8SSW_->Write();
  C8SSW_->Write();
  M45bestallSSW_->Write();
  Chi2extallSSW_->Write();
  DPbballSSW_->Write();
  SumHED4SSW_->Write();
  SumHPD4SSW_->Write();
  SumHED6SSW_->Write();
  SumHPD6SSW_->Write();
  HEDSSW_->Write();
  HPDSSW_->Write();

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
  M6SSSW_->Write();
  C6SSSW_->Write();
  M8SSSW_->Write();
  C8SSSW_->Write();
  M45bestallSSSW_->Write();
  Chi2extallSSSW_->Write();
  DPbballSSSW_->Write();
  SumHED4SSSW_->Write();
  SumHPD4SSSW_->Write();
  SumHED6SSSW_->Write();
  SumHPD6SSSW_->Write();
  HEDSSSW_->Write();
  HPDSSSW_->Write();

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

  // Likelihood
  // ----------
  L_->Write();
  LS_->Write();
  LSS_->Write();
  LSSS_->Write();
  LW_->Write();
  LSW_->Write();
  LSSW_->Write();
  LSSSW_->Write();
  LN_->Write();
  LSN_->Write();
  LSSN_->Write();
  LSSSN_->Write();

}

// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(TDAna);
