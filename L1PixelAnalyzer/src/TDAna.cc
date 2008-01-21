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
//    x compute minimum dijet mass in triplet which gives best top mass
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
// 5) Also, NB: corrected jet energies appear overestimated by a good 15%. Because of 
//    this, the top mass window in Trecfrac is moved to 200 GeV, and the fits likewise
//    (H window is set at 123). See MTbest, MHbest. MWbest indicates the W peaks at 94.
//    These values are used for the reference points in the extraction of chisquares.
//
// 6) Note: the tag rate matrix should, if possible, be redone to account for the changes,
//    and if so, some smart way of taking into account correlations could be envisioned.
//    (changes: etamax 2.5->3.0 and etmin 30->25).
// 
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
  for ( int i=0; i<1000; i++ ) {
    H[i]=0.;
    Hnot[i]=0.;
  }
  for ( int i=0; i<1000; i++ ) {
    T[i]=0.;
    Tnot[i]=0.;
  }

  // Parameters of Et correction functions
  // -------------------------------------
//   upar[2] =  1.512588;
//   upar[1] = -0.05629097;
//   upar[0] = -0.0002950658;
//   tpar[2] =  9.017700;
//   tpar[1] = -0.08954207;
//   tpar[0] = -0.0003376759;
  // Pol6 parameters:
  // ----------------
  // b-tags
  // ------
  //  1  p0           6.55150e-15   1.09289e-16   9.16348e-20   3.26447e+08
  //  2  p1          -2.07528e-11   1.45156e-13   7.84133e-17  -1.94822e+05
  //  3  p2           2.42585e-08   1.04566e-10   6.14276e-14   1.60783e+02
  //  4  p3          -1.28471e-05   6.17959e-08   3.91532e-11   7.25911e-02
  //  5  p4           2.65480e-03   2.13231e-05   1.54571e-08  -1.39745e-03
  //  6  p5          -3.66889e-01   3.42189e-03   3.02745e-06  -8.47267e-06
  //  7  p6           1.69886e+01   1.77644e-01   3.25881e-04  -5.99603e-08
  tpar[0] = 6.55150e-15;
  tpar[1] = -2.07528e-11;
  tpar[2] = 2.42585e-08;
  tpar[3] = -1.28471e-05;
  tpar[4] = 2.65480e-03;
  tpar[5] = -3.66889e-01;
  tpar[6] = 1.69886e+01;
  // non-b-tags
  // ----------
  //  1  p0           4.23909e-15   5.48031e-17   4.06693e-20   2.68032e+13
  //  2  p1          -1.78658e-11   4.85489e-14   3.15999e-17  -1.05256e+11
  //  3  p2           2.52304e-08   3.86624e-11   2.44879e-14   1.49660e+08
  //  4  p3          -1.49399e-05   2.84783e-08   1.88330e-11  -9.48611e+04
  //  5  p4           3.37074e-03   1.53051e-05   1.35516e-08   3.03554e+01
  //  6  p5          -3.93742e-01   2.88642e-03   4.43593e-06  -1.91399e-02
  //  7  p6           1.05092e+01   1.42165e-01   3.74549e-04   5.74116e-05
  upar[0] = 4.23909e-15;
  upar[1] = -1.78658e-11;
  upar[2] = 2.52304e-08;
  upar[3] = -1.49399e-05;
  upar[4] = 3.37074e-03;
  upar[5] = -3.93742e-01;
  upar[6] = 1.05092e+01;

  // Load file with smoothed histograms to use for Likelihood definition
  // -------------------------------------------------------------------
  TFile * FunctionFile = new TFile ("functionfileSS.root");
  FunctionFile->cd();
  HSS_sig[0] = dynamic_cast<TH1D*> ( FunctionFile->Get("C8SS_sigS"));
  HSS_bgr[0] = dynamic_cast<TH1D*> ( FunctionFile->Get("C8SS_bgrS"));
  HSS_sig[1] = dynamic_cast<TH1D*> ( FunctionFile->Get("M8SS_sigS"));
  HSS_bgr[1] = dynamic_cast<TH1D*> ( FunctionFile->Get("M8SS_bgrS"));
  HSS_sig[2] = dynamic_cast<TH1D*> ( FunctionFile->Get("C6SS_sigS"));
  HSS_bgr[2] = dynamic_cast<TH1D*> ( FunctionFile->Get("C6SS_bgrS"));
  HSS_sig[3] = dynamic_cast<TH1D*> ( FunctionFile->Get("GoodHtSS_sigS"));
  HSS_bgr[3] = dynamic_cast<TH1D*> ( FunctionFile->Get("GoodHtSS_bgrS"));
  HSS_sig[4] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSigSS_sigS"));
  HSS_bgr[4] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtSigSS_bgrS"));
  HSS_sig[5] = dynamic_cast<TH1D*> ( FunctionFile->Get("ThdetaSS_sigS"));
  HSS_bgr[5] = dynamic_cast<TH1D*> ( FunctionFile->Get("ThdetaSS_bgrS"));
  HSS_sig[6] = dynamic_cast<TH1D*> ( FunctionFile->Get("HbestcombSS_sigS"));
  HSS_bgr[6] = dynamic_cast<TH1D*> ( FunctionFile->Get("HbestcombSS_bgrS"));
  HSS_sig[7] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtDP2SS_sigS"));
  HSS_bgr[7] = dynamic_cast<TH1D*> ( FunctionFile->Get("MEtDP2SS_bgrS"));
  HSS_sig[8] = dynamic_cast<TH1D*> ( FunctionFile->Get("M5SS_sigS"));
  HSS_bgr[8] = dynamic_cast<TH1D*> ( FunctionFile->Get("M5SS_bgrS"));
  HSS_sig[9] = dynamic_cast<TH1D*> ( FunctionFile->Get("M_othersSS_sigS"));
  HSS_bgr[9] = dynamic_cast<TH1D*> ( FunctionFile->Get("M_othersSS_bgrS"));
  HSS_sig[10]= dynamic_cast<TH1D*> ( FunctionFile->Get("Et6SS_sigS"));
  HSS_bgr[10]= dynamic_cast<TH1D*> ( FunctionFile->Get("Et6SS_bgrS"));
  //  FunctionFile->Close();

  // File with HED and HPD distributions from QCD jets
  // -------------------------------------------------
  TFile * HEDFile = new TFile ("HEDn.root");
  HEDFile->cd();
  HEDpdf[0] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_0"));
  HPDpdf[0] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_0"));
  HEDpdf[1] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_1"));
  HPDpdf[1] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_1"));
  HEDpdf[2] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_2"));
  HPDpdf[2] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_2"));
  HEDpdf[3] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_3"));
  HPDpdf[3] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_3"));
  HEDpdf[4] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_4"));
  HPDpdf[4] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_4"));
  HEDpdf[5] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_5"));
  HPDpdf[5] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_5"));
  HEDpdf[6] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_6"));
  HPDpdf[6] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_6"));
  HEDpdf[7] = dynamic_cast<TH1D*> ( HEDFile->Get("HED_7"));
  HPDpdf[7] = dynamic_cast<TH1D*> ( HEDFile->Get("HPD_7"));
 // HEDFile->Close();

  // File with Total tag mass distributions (S1) from QCD jets
  // ---------------------------------------------------------
  TFile * TagMassFileS1= new TFile ("TTMnS1.root");
  TagMassFileS1->cd();
  MTS1pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_0"));
  MNS1pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_0"));
  MTS1pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_1"));
  MNS1pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_1"));
  MTS1pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_2"));
  MNS1pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_2"));
  MTS1pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_3"));
  MNS1pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_3"));
  MTS1pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_4"));
  MNS1pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_4"));
  MTS1pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_5"));
  MNS1pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_5"));
  MTS1pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_6"));
  MNS1pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_6"));
  MTS1pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MTS1_7"));
  MNS1pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS1->Get("MNS1_7"));
  // TagMassFileS1->Close();
  // File with Total tag mass distributions (S1) from QCD jets
  // ---------------------------------------------------------
  TFile * TagMassFileS2 = new TFile ("TTMnS2.root");
  TagMassFileS2->cd();
  MTS2pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_0"));
  MNS2pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_0"));
  MTS2pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_1"));
  MNS2pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_1"));
  MTS2pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_2"));
  MNS2pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_2"));
  MTS2pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_3"));
  MNS2pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_3"));
  MTS2pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_4"));
  MNS2pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_4"));
  MTS2pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_5"));
  MNS2pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_5"));
  MTS2pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_6"));
  MNS2pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_6"));
  MTS2pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MTS2_7"));
  MNS2pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS2->Get("MNS2_7"));
  // TagMassFileS2->Close();
  // File with Total tag mass distributions (S1) from QCD jets
  // ---------------------------------------------------------
  TFile * TagMassFileS3 = new TFile ("TTMnS3.root");
  TagMassFileS3->cd();
  MTS3pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_0"));
  MNS3pdf[0] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_0"));
  MTS3pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_1"));
  MNS3pdf[1] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_1"));
  MTS3pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_2"));
  MNS3pdf[2] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_2"));
  MTS3pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_3"));
  MNS3pdf[3] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_3"));
  MTS3pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_4"));
  MNS3pdf[4] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_4"));
  MTS3pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_5"));
  MNS3pdf[5] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_5"));
  MTS3pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_6"));
  MNS3pdf[6] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_6"));
  MTS3pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MTS3_7"));
  MNS3pdf[7] =  dynamic_cast<TH1D*> ( TagMassFileS3->Get("MNS3_7"));
  // TagMassFileS3->Close();

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

  // Histograms of delta Et vs Et for jet energy corrections
  // -------------------------------------------------------
  DEtb_prof_ = new TProfile ( "DEtb_prof", "Et mismeasurement vs Et - b-tagged jets", 100, 0., 1000., 
			      -1000., 1000. );
  DEtq_prof_ = new TProfile ( "DEtq_prof", "Et mismeasurement vs Et - non b-tagged jets", 100, 0., 1000., 
			      -1000., 1000. );
  DEtcb_prof_ = new TProfile ( "DEtcb_prof", "Et mismeasurement vs Etc - b-tagged jets", 100, 0., 1000., 
			      -1000., 1000. );
  DEtcq_prof_ = new TProfile ( "DEtcq_prof", "Et mismeasurement vs Etc - non b-tagged jets", 100, 0., 1000., 
			      -1000., 1000. );

  // Histograms for jet-parton matching
  // ----------------------------------
  Drmax_ = new TH2D ( "Drmax", "Drmax for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Drmedall_ = new TH2D ( "Drmedall", "Drmedall for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Drmed07_ = new TH2D ( "Drmed07", "Drmed07 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  N07_ = new TH2D ( "N07", "N07 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  N04_ = new TH2D ( "N04", "N04 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  N02_ = new TH2D ( "N02", "N02 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Nlo_ = new TH2D ( "Nlo", "Nlo for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Detmedall_ = new TH2D ( "Detmedall", "Detmedall for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Detmed07_ = new TH2D ( "Detmed07", "Detmed07 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Perf07_ = new TH2D ( "Perf07", "Perf07 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Perf04_ = new TH2D ( "Perf04", "Perf04 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Perf02_ = new TH2D ( "Perf02", "Perf02 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Det2medall_ = new TH2D ( "Det2medall", "Det2medall for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Det2med07_ = new TH2D ( "Det2med07", "Det2med07 for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7  );
  Hrecfrac_ = new TH2D ( "Hrecfrac", "Hrecfrac for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7 );
  Trecfrac_ = new TH2D ( "Trecfrac", "Trecfrac for choices of etmin, etamax", 21, 14, 56, 11, 1.5, 3.7 );

  // Histograms for MC kinematics reconstruction
  // -------------------------------------------
  MHbest_ = new TH1D ( "MHbest", "Best H mass from jets ass to H partons", 50, 0, 500 );
  MTbest_ = new TH1D ( "MTbest", "Best T mass from jets ass to T partons", 50, 0, 800 );
  MWbest_ = new TH1D ( "MWbest", "Best W mass from jets ass to W partons", 50, 0, 200 );
  HBJ_etrank_ = new TH1D ( "HBJ_etrank", "Rank in Et of b-jets from H decay", 8, -0.5, 7.5 );
  Hpt_ = new TH1D ( "Hpt", "Pt of best H mass", 50, 0, 500 );
  Heta_ = new TH1D ( "Heta", "Eta of H from jets", 50, 0., 4.);
  Hdr_ = new TH1D ( "Hdr", "Delta R of b-jets from H", 50, 0., 5. );
  MHnot_ = new TH1D ( "MHnot", "H mass from jets NOT ass to H partons", 50, 0, 500 ); 
  Hnotpt_ = new TH1D ( "Hnotpt", "Pt of b jets NOT from H", 50, 0, 500 );
  Hnoteta_ = new TH1D ( "Hnoteta", "Eta of H from NOT H jets", 50, 0., 4.);
  Hnotdr_ = new TH1D ( "Hnotdr", "Delta R of b-jets NOT from H", 50, 0., 5. );
  MTnotbest_ = new TH1D ( "MTnotbest", "Best T mass from jets NOT ass to T partons", 50, 0, 800 );
  Tpt_ = new TH1D ( "Tpt", "Pt of jet triplet ass to T partons", 50, 0, 600 );
  Teta_ = new TH1D ( "Teta", "Eta of jet triplet ass to T partons", 50, 0., 4. );
  THdeta_ = new TH1D ( "THdeta", "Eta difference of T and H", 50, -5., 5. );
  THdphi_ = new TH1D ( "THdphi", "Phi difference of T and H", 50, 0., 3.15 );
  THproj_ = new TH1D ( "THproj", "Momentum of T projected on H direction", 50, -1000., 1000. );
  Tnotpt_ = new TH1D ( "Tnotpt", "Pt of jet triplet NOT ass to T partons", 50, 0, 600 );
  Tnoteta_ = new TH1D ( "Tnoteta", "Eta of jet triplet NOT ass to T partons", 50, 0., 4. );
  THnotdeta_ = new TH1D ( "THnotdeta", "Eta difference of NOT T and H", 50, -5., 5. );
  THnotdphi_ = new TH1D ( "THnotdphi", "Phi difference of NOT T and H", 50, 0., 3.15 );
  THnotproj_ = new TH1D ( "THnotproj", "Momentum of NOT T projected on H direction", 50, -1000., 1000. );

  // Histograms used to obtain shapes of SumHED and SumHPD
  // for tag matrix, by summing QCD shapes (see HEDn.C)
  // -----------------------------------------------------
  HED1_ = new TH1D ( "HED1", "HED 1", 200, 0., 40. );
  HPD1_ = new TH1D ( "HPD1", "HPD 1", 200, 0., 40. );
  HED2_ = new TH1D ( "HED2", "HED 2", 200, 0., 40. );
  HPD2_ = new TH1D ( "HPD2", "HPD 2", 200, 0., 40. );
  HED3_ = new TH1D ( "HED3", "HED 3", 200, 0., 40. );
  HPD3_ = new TH1D ( "HPD3", "HPD 3", 200, 0., 40. );
  HED4_ = new TH1D ( "HED4", "HED 4", 200, 0., 40. );
  HPD4_ = new TH1D ( "HPD4", "HPD 4", 200, 0., 40. );
  HED5_ = new TH1D ( "HED5", "HED 5", 200, 0., 40. );
  HPD5_ = new TH1D ( "HPD5", "HPD 5", 200, 0., 40. );
  HED6_ = new TH1D ( "HED6", "HED 6", 200, 0., 40. );
  HPD6_ = new TH1D ( "HPD6", "HPD 6", 200, 0., 40. );
  HED7_ = new TH1D ( "HED7", "HED 7", 200, 0., 40. );
  HPD7_ = new TH1D ( "HPD7", "HPD 7", 200, 0., 40. );
  HED8_ = new TH1D ( "HED8", "HED 8", 200, 0., 40. );
  HPD8_ = new TH1D ( "HPD8", "HPD 8", 200, 0., 40. );

  // Histograms for different levels of selection
  // --------------------------------------------
  NJets_ = new TH1D ( "NJets", "Number of selected jets", 50, 0, 50 );
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
  M3best_ = new TH1D ( "M3best", "M3best", 50, 0., 800. );
  Mwbest_ = new TH1D ( "Mwbest", "Mwbest", 50, 0., 500. );
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
  Et6_ = new TH1D ( "Et6", "Et of sixth jet", 50, 0., 150. );
  Mwmin_ = new TH1D ( "Mwmin", "Minimum dijet mass in best triplet", 50, 0., 200. );
  Hbestcomb_ = new TH1D ( "Hbestcomb", "Best H mass combination", 50, 0., 400. );
  Drpairbestall_ = new TH1D ( "Drpairbestall", "Dr of best bb pair", 50, 0., 4. );
  M3a_ = new TH1D ( "M3a", "Highest Et triplet besides best H combo", 50, 0., 500. );
  Mwa_ = new TH1D ( "Mwa", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  Scprod_ = new TH1D ( "Scprod", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  Thdeta_ = new TH1D ( "Thdeta", "Delta eta top-higgs", 50, -5., 5. );
  M5_ = new TH1D ( "M5", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1_ = new TH1D ( "TTMS1", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2_ = new TH1D ( "TTMS2", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3_ = new TH1D ( "TTMS3", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsN_ = new TH1D ( "NJetsN", "Number of selected jets", 50, 0, 50 );
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
  M3bestN_ = new TH1D ( "M3bestN", "M3best", 50, 0., 800. );
  MwbestN_ = new TH1D ( "MwbestN", "Mwbest", 50, 0., 500. );
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
  Et6N_ = new TH1D ( "Et6N", "Et of sixth jet", 50, 0., 150. );
  MwminN_ = new TH1D ( "MwminN", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombN_ = new TH1D ( "HbestcombN", "Best H mass combination", 50, 0., 400. );
  DrpairbestallN_ = new TH1D ( "DrpairbestallN", "Dr of best bb pair", 50, 0., 4. );
  M3aN_ = new TH1D ( "M3aN", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaN_ = new TH1D ( "MwaN", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodN_ = new TH1D ( "ScprodN", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaN_ = new TH1D ( "ThdetaN", "Delta eta top-higgs", 50, -5., 5. );
  M5N_ = new TH1D ( "M5N", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1N_ = new TH1D ( "TTMS1N", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2N_ = new TH1D ( "TTMS2N", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3N_ = new TH1D ( "TTMS3N", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsS_ = new TH1D ( "NJetsS", "Number of selected jets", 50, 0, 50 );
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
  M3bestS_ = new TH1D ( "M3bestS", "M3bestS", 50, 0., 800. );
  MwbestS_ = new TH1D ( "MwbestS", "MwbestS", 50, 0., 500. );
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
  Et6S_ = new TH1D ( "Et6S", "Et of sixth jet", 50, 0., 150. );
  MwminS_ = new TH1D ( "MwminS", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombS_ = new TH1D ( "HbestcombS", "Best H mass combination", 50, 0., 400. );
  DrpairbestallS_ = new TH1D ( "DrpairbestallS", "Dr of best bb pair", 50, 0., 4. );
  M3aS_ = new TH1D ( "M3aS", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaS_ = new TH1D ( "MwaS", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodS_ = new TH1D ( "ScprodS", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaS_ = new TH1D ( "ThdetaS", "Delta eta top-higgs", 50, -5., 5. );
  M5S_ = new TH1D ( "M5S", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1S_ = new TH1D ( "TTMS1S", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2S_ = new TH1D ( "TTMS2S", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3S_ = new TH1D ( "TTMS3S", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSN_ = new TH1D ( "NJetsSN", "Number of selected jets", 50, 0, 50 );
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
  M3bestSN_ = new TH1D ( "M3bestSN", "M3bestS", 50, 0., 800. );
  MwbestSN_ = new TH1D ( "MwbestSN", "MwbestS", 50, 0., 500. );
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
  Et6SN_ = new TH1D ( "Et6SN", "Et of sixth jet", 50, 0., 150. );
  MwminSN_ = new TH1D ( "MwminSN", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSN_ = new TH1D ( "HbestcombSN", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSN_ = new TH1D ( "DrpairbestallSN", "Dr of best bb pair", 50, 0., 4. );
  M3aSN_ = new TH1D ( "M3aSN", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSN_ = new TH1D ( "MwaSN", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSN_ = new TH1D ( "ScprodSN", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSN_ = new TH1D ( "ThdetaSN", "Delta eta top-higgs", 50, -5., 5. );
  M5SN_ = new TH1D ( "M5SN", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SN_ = new TH1D ( "TTMS1SN", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SN_ = new TH1D ( "TTMS2SN", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SN_ = new TH1D ( "TTMS3SN", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSS_ = new TH1D ( "NJetsSS", "Number of selected jets", 50, 0, 50 );
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
  M3bestSS_ = new TH1D ( "M3bestSS", "M3bestSS", 50, 0., 800. );
  MwbestSS_ = new TH1D ( "MwbestSS", "MwbestSS", 50, 0., 500. );
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
  Et6SS_ = new TH1D ( "Et6SS", "Et of sixth jet", 50, 0., 150. );
  MwminSS_ = new TH1D ( "MwminSS", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSS_ = new TH1D ( "HbestcombSS", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSS_ = new TH1D ( "DrpairbestallSS", "Dr of best bb pair", 50, 0., 4. );
  M3aSS_ = new TH1D ( "M3aSS", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSS_ = new TH1D ( "MwaSS", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSS_ = new TH1D ( "ScprodSS", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSS_ = new TH1D ( "ThdetaSS", "Delta eta top-higgs", 50, -5., 5. );
  M5SS_ = new TH1D ( "M5SS", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SS_ = new TH1D ( "TTMS1SS", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SS_ = new TH1D ( "TTMS2SS", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SS_ = new TH1D ( "TTMS3SS", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSSN_ = new TH1D ( "NJetsSSN", "Number of selected jets", 50, 0, 50 );
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
  M3bestSSN_ = new TH1D ( "M3bestSSN", "M3bestSSN", 50, 0., 800. );
  MwbestSSN_ = new TH1D ( "MwbestSSN", "MwbestSSN", 50, 0., 500. );
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
  Et6SSN_ = new TH1D ( "Et6SSN", "Et of sixth jet", 50, 0., 150. );
  MwminSSN_ = new TH1D ( "MwminSSN", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSSN_ = new TH1D ( "HbestcombSSN", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSSN_ = new TH1D ( "DrpairbestallSSN", "Dr of best bb pair", 50, 0., 4. );
  M3aSSN_ = new TH1D ( "M3aSSN", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSSN_ = new TH1D ( "MwaSSN", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSSN_ = new TH1D ( "ScprodSSN", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSSN_ = new TH1D ( "ThdetaSSN", "Delta eta top-higgs", 50, -5., 5. );
  M5SSN_ = new TH1D ( "M5SSN", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SSN_ = new TH1D ( "TTMS1SSN", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SSN_ = new TH1D ( "TTMS2SSN", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SSN_ = new TH1D ( "TTMS3SSN", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSSS_ = new TH1D ( "NJetsSSS", "Number of selected jets", 50, 0, 50 );
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
  M3bestSSS_ = new TH1D ( "M3bestSSS", "M3bestSSS", 50, 0., 800. );
  MwbestSSS_ = new TH1D ( "MwbestSSS", "MwbestSSS", 50, 0., 500. );
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
  Et6SSS_ = new TH1D ( "Et6SSS", "Et of sixth jet", 50, 0., 150. );
  MwminSSS_ = new TH1D ( "MwminSSS", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSSS_ = new TH1D ( "HbestcombSSS", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSSS_ = new TH1D ( "DrpairbestallSSS", "Dr of best bb pair", 50, 0., 4. );
  M3aSSS_ = new TH1D ( "M3aSSS", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSSS_ = new TH1D ( "MwaSSS", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSSS_ = new TH1D ( "ScprodSSS", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSSS_ = new TH1D ( "ThdetaSSS", "Delta eta top-higgs", 50, -5., 5. );
  M5SSS_ = new TH1D ( "M5SSS", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SSS_ = new TH1D ( "TTMS1SSS", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SSS_ = new TH1D ( "TTMS2SSS", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SSS_ = new TH1D ( "TTMS3SSS", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSSSN_ = new TH1D ( "NJetsSSSN", "Number of selected jets", 50, 0, 50 );
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
  M3bestSSSN_ = new TH1D ( "M3bestSSSN", "M3bestSSSN", 50, 0., 800. );
  MwbestSSSN_ = new TH1D ( "MwbestSSSN", "MwbestSSSN", 50, 0., 500. );
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
  Et6SSSN_ = new TH1D ( "Et6SSSN", "Et of sixth jet", 50, 0., 150. );
  MwminSSSN_ = new TH1D ( "MwminSSSN", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSSSN_ = new TH1D ( "HbestcombSSSN", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSSSN_ = new TH1D ( "DrpairbestallSSSN", "Dr of best bb pair", 50, 0., 4. );
  M3aSSSN_ = new TH1D ( "M3aSSSN", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSSSN_ = new TH1D ( "MwaSSSN", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSSSN_ = new TH1D ( "ScprodSSSN", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSSSN_ = new TH1D ( "ThdetaSSSN", "Delta eta top-higgs", 50, -5., 5. );
  M5SSSN_ = new TH1D ( "M5SSSN", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SSSN_ = new TH1D ( "TTMS1SSSN", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SSSN_ = new TH1D ( "TTMS2SSSN", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SSSN_ = new TH1D ( "TTMS3SSSN", "Total tag mass with S3 tracks", 50, 0., 80. );

  // Histograms to check where are the tags in the Et-ordered jet list
  // -----------------------------------------------------------------
  N4NJSSS_ = new TH1D ( "N4NJSSS", "N of 4HEL tags vs N jets", 
			20, 0, 20 );  // These are filled only for this
  E4NJSSS_ = new TH1D ( "E4NJSSS", "Efficiency of 4HEL tags vs N jets", 
			20, 0, 20 );  // selection and do not require W, N

  NJetsW_ = new TH1D ( "NJetsW", "Number of selected jets", 50, 0, 50 );
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
  M3bestW_ = new TH1D ( "M3bestW", "M3best", 50, 0., 800. );
  MwbestW_ = new TH1D ( "MwbestW", "Mwbest", 50, 0., 500. );
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
  Et6W_ = new TH1D ( "Et6W", "Et of sixth jet", 50, 0., 150. );
  MwminW_ = new TH1D ( "MwminW", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombW_ = new TH1D ( "HbestcombW", "Best H mass combination", 50, 0., 400. );
  DrpairbestallW_ = new TH1D ( "DrpairbestallW", "Dr of best bb pair", 50, 0., 4. );
  M3aW_ = new TH1D ( "M3aW", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaW_ = new TH1D ( "MwaW", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodW_ = new TH1D ( "ScprodW", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaW_ = new TH1D ( "ThdetaW", "Delta eta top-higgs", 50, -5., 5. );
  M5W_ = new TH1D ( "M5W", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1W_ = new TH1D ( "TTMS1W", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2W_ = new TH1D ( "TTMS2W", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3W_ = new TH1D ( "TTMS3W", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSW_ = new TH1D ( "NJetsSW", "Number of selected jets", 50, 0, 50 );
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
  M3bestSW_ = new TH1D ( "M3bestSW", "M3bestS", 50, 0., 800. );
  MwbestSW_ = new TH1D ( "MwbestSW", "MwbestS", 50, 0., 500. );
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
  Et6SW_ = new TH1D ( "Et6SW", "Et of sixth jet", 50, 0., 150. );
  MwminSW_ = new TH1D ( "MwminSW", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSW_ = new TH1D ( "HbestcombSW", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSW_ = new TH1D ( "DrpairbestallSW", "Dr of best bb pair", 50, 0., 4. );
  M3aSW_ = new TH1D ( "M3aSW", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSW_ = new TH1D ( "MwaSW", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSW_ = new TH1D ( "ScprodSW", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSW_ = new TH1D ( "ThdetaSW", "Delta eta top-higgs", 50, -5., 5. );
  M5SW_ = new TH1D ( "M5SW", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SW_ = new TH1D ( "TTMS1SW", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SW_ = new TH1D ( "TTMS2SW", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SW_ = new TH1D ( "TTMS3SW", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSSW_ = new TH1D ( "NJetsSSW", "Number of selected jets", 50, 0, 50 );
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
  M3bestSSW_ = new TH1D ( "M3bestSSW", "M3bestSS", 50, 0., 800. );
  MwbestSSW_ = new TH1D ( "MwbestSSW", "MwbestSS", 50, 0., 500. );
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
  Et6SSW_ = new TH1D ( "Et6SSW", "Et of sixth jet", 50, 0., 150. );
  MwminSSW_ = new TH1D ( "MwminSSW", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSSW_ = new TH1D ( "HbestcombSSW", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSSW_ = new TH1D ( "DrpairbestallSSW", "Dr of best bb pair", 50, 0., 4. );
  M3aSSW_ = new TH1D ( "M3aSSW", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSSW_ = new TH1D ( "MwaSSW", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSSW_ = new TH1D ( "ScprodSSW", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSSW_ = new TH1D ( "ThdetaSSW", "Delta eta top-higgs", 50, -5., 5. );
  M5SSW_ = new TH1D ( "M5SSW", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SSW_ = new TH1D ( "TTMS1SSW", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SSW_ = new TH1D ( "TTMS2SSW", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SSW_ = new TH1D ( "TTMS3SSW", "Total tag mass with S3 tracks", 50, 0., 80. );

  NJetsSSSW_ = new TH1D ( "NJetsSSSW", "Number of selected jets", 50, 0, 50 );
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
  M3bestSSSW_ = new TH1D ( "M3bestSSSW", "M3bestSSS", 50, 0., 800. );
  MwbestSSSW_ = new TH1D ( "MwbestSSSW", "MwbestSSS", 50, 0., 500. );
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
  Et6SSSW_ = new TH1D ( "Et6SSSW", "Et of sixth jet", 50, 0., 150. );
  MwminSSSW_ = new TH1D ( "MwminSSSW", "Minimum dijet mass in best triplet", 50, 0., 200. );
  HbestcombSSSW_ = new TH1D ( "HbestcombSSSW", "Best H mass combination", 50, 0., 400. );
  DrpairbestallSSSW_ = new TH1D ( "DrpairbestallSSSW", "Dr of best bb pair", 50, 0., 4. );
  M3aSSSW_ = new TH1D ( "M3aSSSW", "Highest Et triplet besides best H combo", 50, 0., 500. );
  MwaSSSW_ = new TH1D ( "MwaSSSW", "Doublet not tagged in highest Et triplet", 50, 0., 300. );
  ScprodSSSW_ = new TH1D ( "ScprodSSSW", "Scalar product of top vector and H versor", 50, -1200., 1200. );
  ThdetaSSSW_ = new TH1D ( "ThdetaSSSW", "Delta eta top-higgs", 50, -5., 5. );
  M5SSSW_ = new TH1D ( "M5SSSW", "Mass of five jets from t and h", 50, 0., 2000. );
  TTMS1SSSW_ = new TH1D ( "TTMS1SSSW", "Total tag mass with S1 tracks", 50, 0., 80. );
  TTMS2SSSW_ = new TH1D ( "TTMS2SSSW", "Total tag mass with S2 tracks", 50, 0., 80. );
  TTMS3SSSW_ = new TH1D ( "TTMS3SSSW", "Total tag mass with S3 tracks", 50, 0., 80. );

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
  
  // Definitions for b-tagging
  // -------------------------
  loose_  = 2.3; 
  // high eff -> 70.49% b / 32.33% c / 8.64% uds / 10.43% g / 9.98% udsg // P.Schilling 23/10/07
  medium_ = 5.3; 
  // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
  tight_  = 4.8; 
  // high pur -> 31.94% b /  2.93% c / 0.10% uds /  0.11% g / 0.10% udsg // P.Schilling 23/10/07

  // Relative Likelihood
  // -------------------       
  rel_lik = 0.;

  // Likelihood histograms
  // ---------------------
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

  // Read H decay file
  // -----------------
  for ( int i=0; i<1000; i++ ) {
    Hread[i]=0.;
    Hnotread[i]=0.;
  }
  int counter=0;
  ifstream Hreadfile ("Hread.asc");
  string line;
  while (Hreadfile) {
    getline(Hreadfile,line,' ');
    stringstream pippo1(line);
    pippo1 >> Hread[counter];
    getline(Hreadfile,line);
    stringstream pippo2(line);
    pippo1 >> Hnotread[counter];
    counter++;
  }
  cout << "Read " << counter << " lines from H file." << endl;
  double toth=0.;
  double tothnot=0.;
  for ( int i=0; i<1000; i++ ) {
    toth+=Hread[i];
    tothnot+=Hnotread[i];
  }
  for ( int i=0; i<1000; i++ ) {
    Hread[i]=Hread[i]/toth;
    Hnotread[i]=Hnotread[i]/tothnot;
  }

  // Read T decay file
  // -----------------
  for ( int i=0; i<1000; i++ ) {
    Tread[i]=0.;
    Tnotread[i]=0.;
  }
  counter=0;
  ifstream Treadfile ("Tread.asc");
  while (Treadfile) {
    getline(Treadfile,line,' ');
    stringstream pippo1(line);
    pippo1 >> Tread[counter];
    getline(Treadfile,line);
    stringstream pippo2(line);
    pippo1 >> Tnotread[counter];
    counter++;
  }
  cout << "Read " << counter << " lines from T file." << endl;
  double tott=0.;
  double tottnot=0.;
  for ( int i=0; i<1000; i++ ) {
    tott+=Tread[i];
    tottnot+=Tnotread[i];
  }
  for ( int i=0; i<1000; i++ ) {
    Tread[i]=Tread[i]/tott;
    Tnotread[i]=Tnotread[i]/tottnot;
  }

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

  // End of initializations
  // ----------------------
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

  // General definitions
  // -------------------
  double mhref=123.; 
  double mtref=200.;   // Massa top ricostruita, >172 per via dello shift in JEcorrs.
  double mwref=94.;

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
    double JPet[100]={0.};
    double JPeta[100]={0.};
    double JPphi[100]={0.};
    double JPetc[100]={0.};
    bool JPtag[100]={false};
    int NJP=0;
    for ( OfflineJetCollection::const_iterator cal = caloJets->begin(); 
	  cal != caloJets->end(); ++cal ) {
      if ( NJP<100 ) {
	JPet[NJP]=cal->et();
	JPeta[NJP]=cal->eta();
	JPphi[NJP]=cal->phi();
	if ( cal->discriminatorHighEff()>loose_ ) JPtag[NJP]=true;
	// Correct jet Et
	// --------------
	if ( JPtag[NJP] ) {
	  JPetc[NJP]=JPet[NJP]+
	    tpar[0]*pow(JPet[NJP],6)+
	    tpar[1]*pow(JPet[NJP],5)+
	    tpar[2]*pow(JPet[NJP],4)+
	    tpar[3]*pow(JPet[NJP],3)+
	    tpar[4]*pow(JPet[NJP],2)+
	    tpar[5]*JPet[NJP]+
	    tpar[6];
	  if ( JPetc[NJP]<0. ) JPetc[NJP]= 0.;
	} else {
	  JPetc[NJP]=JPet[NJP]+
	    upar[0]*pow(JPet[NJP],6)+
	    upar[1]*pow(JPet[NJP],5)+
	    upar[2]*pow(JPet[NJP],4)+
	    upar[3]*pow(JPet[NJP],3)+
	    upar[4]*pow(JPet[NJP],2)+
	    upar[5]*JPet[NJP]+
	    upar[6];
	  if ( JPetc[NJP]<0. ) JPetc[NJP]= 0.;
	}
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
	  bool   a4 = JPtag[k];
	  double a5 = JPetc[k];
	  JPet[k]=JPet[k-1];
	  JPeta[k]=JPeta[k-1];
	  JPphi[k]=JPphi[k-1];
	  JPtag[k]=JPtag[k-1];
	  JPetc[k]=JPetc[k-1];
	  JPet[k-1]=a1;
	  JPeta[k-1]=a2;
	  JPphi[k-1]=a3;
	  JPtag[k-1]=a4;
	  JPetc[k-1]=a5;
	}
      }
    }

    // Construct a pt res function
    // ---------------------------
    bool ijass[20]={false};
    int nass=0;
    for ( int ip=0; ip<iparton; ip++ ) {
      double drmin=0.04;
      int tmpind=-1;
      for ( int ij=0; ij<NJP; ij++ ) {
	if ( JPet[ij]>25. && fabs(JPeta[ij])<3.0 && !ijass[ij] ) {
	  double deta = Parton_eta[ip]-JPeta[ij];
	  double dphi = 3.1415926-fabs(fabs(Parton_phi[ip]-JPphi[ij])-3.1415926);
	  double dr2  = deta*deta+dphi*dphi;
	  if ( dr2<drmin ) {
	    drmin=dr2;
	    tmpind=ij;
	  }
	}
      }
      if ( tmpind>-1 && nass<20 ) {
	ijass[nass]=true;
	nass++;
	if ( JPtag[tmpind] ) {
	  DEtb_prof_->Fill(JPet[tmpind],Parton_pt[ip]-JPet[tmpind]);
	  DEtcb_prof_->Fill(JPetc[tmpind],Parton_pt[ip]-JPetc[tmpind]);
	} else {
	  DEtq_prof_->Fill(JPet[tmpind],Parton_pt[ip]-JPet[tmpind]);
	  DEtcq_prof_->Fill(JPetc[tmpind],Parton_pt[ip]-JPetc[tmpind]);
	}
      }
    }

    // Understand if it is the t or the tbar the one going to 3 jets
    // -------------------------------------------------------------
    bool thad=false;
    for ( int ip=0; ip<iparton; ip++ ) {
      if ( Parton_dec[ip]==1 ) thad=true;
    }
    // Study parton-jet matching
    // -------------------------
    bool ipass[8]; // whether a parton is associated to one of the leading jets
    for ( int ietmin=0; ietmin<21; ietmin++ ) {
      double etmin = 15.+2.*(double)ietmin;
      for ( int ietamax=0; ietamax<11; ietamax++ ) {
	double etamax = 1.6+0.2*(double)ietamax;
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
	int njh=0;
	int ijt[3]={0};
	int njt=0;
	int ijw[2]={0};
	int njw=0;
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
		double d_et = fabs(Parton_pt[ip]-JPetc[ij])/Parton_pt[ip];
		double dr2  = deta*deta+dphi*dphi; // +d_et*d_et;
		if ( dr2<drmin ) {
		  drmin=dr2;
		  tmpind=ip;
		  det=(-Parton_pt[ip]+JPetc[ij])/Parton_pt[ip];
		  det2=pow((Parton_pt[ip]-JPetc[ij])/Parton_pt[ip],2);
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
	      if ( (  thad && ( Parton_dec[tmpind]==1 || Parton_dec[tmpind]==2 ) ) ||
		   ( !thad && ( Parton_dec[tmpind]==3 || Parton_dec[tmpind]==4 ) ) ) {
		     // Special treatment for H daughters
		     // ---------------------------------
		ijt[njt]=ij;
		njt++;
		if ( Parton_dec[tmpind]==1 || Parton_dec[tmpind]==3 ) {
		  // W decay daughters
		  // -----------------
		  ijw[njw]=ij;
		  njw++;
		}
	      }
	    } else { 
	      // Count here how many jets among the first N (at most equal to the number of partons)
	      // could not be associated reasonably to a parton
	      // -----------------------------------------------------------------------------------
	      if ( considered<iparton ) leftover++;
	    }
	  } // if there is a jet
	} // ij
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
	// Mass of two jets from h when tagged
	// -----------------------------------
	double m45=0.;
	double etapair=0.;
	double phipair=0.;
	double pxh=0.;
	double pyh=0.;
	double pzh=0.;
	double ph2=0.;
	if ( njh==2 && njt==3 && JPtag[ijh[0]] && JPtag[ijh[1]] ) {
	  double px4=JPetc[ijh[0]]*cos(JPphi[ijh[0]]);
	  double px5=JPetc[ijh[1]]*cos(JPphi[ijh[1]]);
	  double py4=JPetc[ijh[0]]*sin(JPphi[ijh[0]]);
	  double py5=JPetc[ijh[1]]*sin(JPphi[ijh[1]]);
	  double pz4=JPetc[ijh[0]]/tan(2*atan(exp(-JPeta[ijh[0]])));
	  double pz5=JPetc[ijh[1]]/tan(2*atan(exp(-JPeta[ijh[1]])));
	  double e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
	  double e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
	  pxh=px4+px5;
	  pyh=py4+py5;
	  pzh=pz4+pz5;
	  ph2=pxh*pxh+pyh*pyh+pzh*pzh; 
	  m45=(e4+e5)*(e4+e5)-ph2;
	  if ( m45>0 ) m45=sqrt(m45);
	  if ( ietmin==5 && ietamax==7 ) { 
	    double ptpair = sqrt(pxh*pxh+pyh*pyh);
	    double deta = JPeta[ijh[0]]-JPeta[ijh[1]];
	    double dphi = 3.1415926-fabs(fabs(JPphi[ijh[0]]-JPphi[ijh[1]])-3.1415926);
	    double drpair = sqrt(deta*deta+dphi*dphi);
	    double ptotpair = sqrt(pxh*pxh+pyh*pyh+pzh*pzh);
	    double etapair=0.;
	    if ( ptotpair-pzh!=0 ) etapair = 0.5*log((ptotpair+pzh)/(ptotpair-pzh));
	    phipair = atan2(pyh,pxh);
	    MHbest_->Fill(m45);
	    Hpt_->Fill(ptpair);
	    Heta_->Fill(fabs(etapair));
	    Hdr_->Fill(drpair);
	    // Fill matrix with probability of pair
	    // ------------------------------------
	    int ipt=(int)((ptpair/500.)*10);
	    int ieta=(int)((fabs(etapair)/4.)*10);
	    int idr=(int)((drpair/4.)*10);
	    if ( ipt<0 ) ipt=0;
	    if ( ipt>9 ) ipt=9;
	    if ( ieta<0 ) ieta=0;
	    if ( ieta>9) ieta=9;
	    if ( idr<0 ) idr=0;
	    if ( idr>9 ) idr=9;
	    H[idr*100+ieta*10+ipt]++;
	  }
	}

	// Only for et,eta cuts chosen: study which is most likely b-tag pair
	// for events where both h products have been b-tagged
	// ------------------------------------------------------------------
	if ( njh==2 && njt==3 && JPtag[ijh[0]] && JPtag[ijh[1]] && ietmin==5 && ietamax==7 ) {
	  considered=0; // Consider at most 8 jets for the association
	  // First get b-tag indices:
	  // ------------------------
	  int indtag=0;
	  int ijtag[100]={0.};
	  for ( int ij=0; ij<NJP; ij++ ) {
	    if ( JPet[ij]>etmin && fabs(JPeta[ij])<etamax && considered<8 ) {
	      considered++;
	      if ( JPtag[ij] ) {
		ijtag[indtag]=ij;
		indtag++;
	      }
	    }
	  }
	  // Now find rank in Et
	  // -------------------
	  if ( indtag>=4 ) {
	    for ( int ind=0; ind<indtag; ind++ ) {
	      if ( ijh[0]==ijtag[ind] || ijh[1]==ijtag[ind] ) {
		HBJ_etrank_->Fill((double)ind);
	      }
	    }
	  }
	  // Fill pt, eta, dr distribution for wrong tagged pairs
	  // ----------------------------------------------------
	  for ( int ind1=0; ind1<indtag-1; ind1++ ) {
	    for ( int ind2=ind1+1; ind2<indtag; ind2++ ) {
	      if ( ( ijh[0]!=ind1 && ijh[0]!=ind2 ) ||  
		   ( ijh[1]!=ind1 && ijh[1]!=ind2 ) ) {
		double px4=JPetc[ind1]*cos(JPphi[ind1]);
		double px5=JPetc[ind2]*cos(JPphi[ind2]);
		double py4=JPetc[ind1]*sin(JPphi[ind1]);
		double py5=JPetc[ind2]*sin(JPphi[ind2]);
		double pz4=JPetc[ind1]/tan(2*atan(exp(-JPeta[ind1])));
		double pz5=JPetc[ind2]/tan(2*atan(exp(-JPeta[ind2])));
		double e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
		double e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
		double m45not=(e4+e5)*(e4+e5)-(px4+px5)*(px4+px5)-(py4+py5)*(py4+py5)-(pz4+pz5)*(pz4+pz5);
		if ( m45not>0 ) m45not=sqrt(m45not);
		double ptpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5));
		double deta = JPeta[ind1]-JPeta[ind2];
		double dphi = 3.1415926-fabs(fabs(JPphi[ind1]-JPphi[ind2])-3.1415926);
		double drpair = sqrt(deta*deta+dphi*dphi);
		double ptotpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5)+(pz4+pz5)*(pz4+pz5));
		double etanotpair=0.;
		if ( ptotpair-pz4-pz5!=0 ) etanotpair = 0.5*log((ptotpair+pz4+pz5)/(ptotpair-pz4-pz5));
		double phinotpair=atan2(py4+py5,px4+px5);
		MHnot_->Fill(m45not);
		Hnotpt_->Fill(ptpair);
		Hnoteta_->Fill(fabs(etanotpair));
		Hnotdr_->Fill(drpair);
		// Fill matrix with probability of pair
		// ------------------------------------
		int ipt=(int)((ptpair/500.)*10);
		int ieta=(int)((fabs(etanotpair)/4.)*10);
		int idr=(int)((drpair/4.)*10);
		if ( ipt<0 ) ipt=0;
		if ( ipt>9 ) ipt=9;
		if ( ieta<0 ) ieta=0;
		if ( ieta>9) ieta=9;
		if ( idr<0 ) idr=0;
		if ( idr>9 ) idr=9;
		Hnot[idr*100+ieta*10+ipt]++;
	      }
	    }
	  }

	} // ietmin=5, ietamax=7, two tags from hbb jets
	
	// Mass of three jets from t
	// -------------------------
	double m123=0.;
 	if ( ietmin==5 && ietamax==7 && njt==3 && njh==2 ) {
	  if (  JPet[ijt[0]]>etmin && fabs(JPeta[ijt[0]])<etamax  &&
		JPet[ijt[1]]>etmin && fabs(JPeta[ijt[1]])<etamax  && 
		JPet[ijt[2]]>etmin && fabs(JPeta[ijt[2]])<etamax ) {
	    double px1=JPetc[ijt[0]]*cos(JPphi[ijt[0]]);
	    double px2=JPetc[ijt[1]]*cos(JPphi[ijt[1]]);
	    double px3=JPetc[ijt[2]]*cos(JPphi[ijt[2]]);
	    double py1=JPetc[ijt[0]]*sin(JPphi[ijt[0]]);
	    double py2=JPetc[ijt[1]]*sin(JPphi[ijt[1]]);
	    double py3=JPetc[ijt[2]]*sin(JPphi[ijt[2]]);
	    double pz1=JPetc[ijt[0]]/tan(2*atan(exp(-JPeta[ijt[0]])));
	    double pz2=JPetc[ijt[1]]/tan(2*atan(exp(-JPeta[ijt[1]])));
	    double pz3=JPetc[ijt[2]]/tan(2*atan(exp(-JPeta[ijt[2]])));
	    double e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
	    double e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
	    double e3 =sqrt(px3*px3+py3*py3+pz3*pz3);
	    double spx=px1+px2+px3;
	    double spy=py1+py2+py3;
	    double spz=pz1+pz2+pz3;
	    m123=pow(e1+e2+e3,2)-pow(spx,2)-pow(spy,2)-pow(spz,2);
	    if ( m123>0 ) m123=sqrt(m123);
	    
	    double pt3 = sqrt(spx*spx+spy*spy);
	    double scprodth=0.;          // projection of top momentum in h direction
	    if ( ph2>0 ) scprodth = (spx*pxh+spy*pyh+spz*pzh)/sqrt(ph2);
	    double ptot3 = sqrt(spx*spx+spy*spy+spz*spz);
	    double eta3=0.;
	    if ( ptot3-spz!=0 ) eta3 = 0.5*log((ptot3+spz)/(ptot3-spz));
	    double phi3 = atan2(spy,spx);
	    double thdphi = 3.141592-fabs(fabs(phi3-phipair)-3.141592);
	    double thdeta = eta3-etapair; // angle between momenta of t and h
	    MTbest_->Fill(m123);
	    Tpt_->Fill(pt3);
	    Teta_->Fill(fabs(eta3));
	    THdeta_->Fill(thdeta);
	    if ( ph2>0 ) THproj_->Fill(scprodth);
	    THdphi_->Fill(thdphi);
	    // Fill matrix with probability of pair
	    // ------------------------------------
 	    int ipt=(int)((pt3/600.)*10);
 	    int idp=(int)((thdphi/3.15)*10);
 	    int iet=(int)((fabs(eta3)/4.)*10);
 	    if ( ipt<0 ) ipt=0;
 	    if ( ipt>9 ) ipt=9;
 	    if ( idp<0 ) idp=0;
 	    if ( idp>9 ) idp=9;
 	    if ( iet<0 ) iet=0;
 	    if ( iet>9 ) iet=9;
 	    T[iet*100+idp*10+ipt]++;
	  }
	  // Mass of two jets from W
	  // -----------------------
	  double m12=0.;
	  if ( njw==2 ) {
	    double px1=JPetc[ijw[0]]*cos(JPphi[ijw[0]]);
	    double px2=JPetc[ijw[1]]*cos(JPphi[ijw[1]]);
	    double py1=JPetc[ijw[0]]*sin(JPphi[ijw[0]]);
	    double py2=JPetc[ijw[1]]*sin(JPphi[ijw[1]]);
	    double pz1=JPetc[ijw[0]]/tan(2*atan(exp(-JPeta[ijw[0]])));
	    double pz2=JPetc[ijw[1]]/tan(2*atan(exp(-JPeta[ijw[1]])));
	    double e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
	    double e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
	    m12=(e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2);
	    if ( m12>0 ) m12=sqrt(m12);
	    if ( ietmin==5 && ietamax==7 ) MWbest_->Fill(m12);
	  }
	  
	  // Mass and kinematics of three jets NOT from t
	  // --------------------------------------------
	  if (  ietmin==5 && ietamax==7 && njt==3 && njh==2 ) {
	    int considered1=0; // Consider at most 8 jets for the association
	    for ( int i1=0; i1<NJP-2; i1++ ) {
	      if ( JPet[i1]>etmin && fabs(JPeta[i1])<etamax && considered1<8 ) {
		considered1++;
		int considered2=0;
		for ( int i2=0; i2<NJP-1; i2++ ) {
		  if ( JPet[i2]>etmin && fabs(JPeta[i2])<etamax && considered2<8 ) {
		    considered2++;
		    int considered3=0;
		    for ( int i3=0; i3<NJP; i3++ ) {
		      if ( JPet[i3]>etmin && fabs(JPeta[i3])<etamax && considered3<8 ) {
			considered3++;
			if ( ( ijt[0]!=i1 && ijt[0]!=i2 && ijt[0]!=i3 ) ||  
			     ( ijt[1]!=i1 && ijt[1]!=i2 && ijt[1]!=i3 ) ||
			     ( ijt[2]!=i1 && ijt[2]!=i2 && ijt[2]!=i3 ) ) {
			  double px1=JPetc[i1]*cos(JPphi[i1]);
			  double px2=JPetc[i2]*cos(JPphi[i2]);
			  double px3=JPetc[i3]*cos(JPphi[i3]);
			  double py1=JPetc[i1]*sin(JPphi[i1]);
			  double py2=JPetc[i2]*sin(JPphi[i2]);
			  double py3=JPetc[i3]*sin(JPphi[i3]);
			  double pz1=JPetc[i1]/tan(2*atan(exp(-JPeta[i1])));
			  double pz2=JPetc[i2]/tan(2*atan(exp(-JPeta[i2])));
			  double pz3=JPetc[i3]/tan(2*atan(exp(-JPeta[i3])));
			  double e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
			  double e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
			  double e3 =sqrt(px3*px3+py3*py3+pz3*pz3);
			  double spx=px1+px2+px3;
			  double spy=py1+py2+py3;
			  double spz=pz1+pz2+pz3;
			  double m123not=pow(e1+e2+e3,2)-pow(spx,2)-pow(spy,2)-pow(spz,2);
			  if ( m123not>0 ) m123not=sqrt(m123not);
			  
			  double pt3 = sqrt(spx*spx+spy*spy);
			  double scprodth=0.;          // projection of top momentum in h direction
			  if ( ph2>0 ) scprodth = (spx*pxh+spy*pyh+spz*pzh)/sqrt(ph2);
			  double ptot3 = sqrt(spx*spx+spy*spy+spz*spz);
			  double eta3=0.;
			  if ( ptot3-spz!=0 ) eta3 = 0.5*log((ptot3+spz)/(ptot3-spz));
			  double phi3 = atan2(spy,spx);
			  double thdphi = 3.141592-fabs(fabs(phi3-phipair)-3.141592);
			  double thdeta = eta3-etapair; // angle between momenta of t and h
			  MTnotbest_->Fill(m123not);
			  Tnotpt_->Fill(pt3);
			  Tnoteta_->Fill(fabs(eta3));
			  THnotdeta_->Fill(thdeta);
			  if ( ph2>0 ) THnotproj_->Fill(scprodth);
			  THnotdphi_->Fill(thdphi);
			  // Fill matrix with probability of pair
			  // ------------------------------------
			  int ipt=(int)((pt3/600.)*10);
			  int idp=(int)((thdphi/3.15)*10);
			  int iet=(int)((fabs(eta3)/4.)*10);
			  if ( ipt<0 ) ipt=0;
			  if ( ipt>9 ) ipt=9;
			  if ( idp<0 ) idp=0;
			  if ( idp>9 ) idp=9;
			  if ( iet<0 ) iet=0;
			  if ( iet>9 ) iet=9;
			  Tnot[iet*100+idp*10+ipt]++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
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
	if ( fabs(m45-mhref)<25 ) Hrecfrac_->Fill(etmin,etamax);
	if ( fabs(m123-mtref)<40 ) Trecfrac_->Fill(etmin,etamax);
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

  double etbounds[9]={25.,50.,70.,90.,110.,150.,200.,300.,10000.};
  double etabounds[5]={0.,1.,1.5,2.0,3.0};
  double ns1bounds[9]={2.,3.,4.,5.,6.,7.,8.,9.,100.};

  double sumet = caloMET->sumEt();
  double met = caloMET->et();
  double metsig = caloMET->mEtSig();
  double metphi = caloMET->phi();

  // Count IC5 jets with Et>=30 GeV and |eta| < 3.0
  // ----------------------------------------------
  double goodIc5Jets = 0.; // Number of selected jets
  double NTagsL = 0.;      // Number of loose b-tags
  double NTagsM = 0.;      // Number of medium b-tags 
  double NTagsT = 0.;      // Number of tight b-tags

  double UncorrSumEx=0;
  double UncorrSumEy=0.;
  double UncorrSumEt=0.;  // SumEt with all jets before corrections
  double UncorrHt=0.;     // Ht with all jets before corrections
  double UncorrHt0=0.;     // Ht with all jets before corrections

  double CorrSumEx=0.;
  double CorrSumEy=0.;
  double CorrSumEt=0.;    // SumEt with all jets after corrections 
  double CorrHt=0.;       // Ht with all jets after corrections
  double CorrHt0=0.;       // Ht with all jets after corrections

  double GoodSumEt=0.;    // SumEt with all selected jets
  double GoodHt=0.;       // Ht with all selected jets
  double GoodHt2=0.;      // Ht with all selected jets and standard MET
  double GoodHt0=0.;       // Ht with all selected jets
  double GoodHt20=0.;      // Ht with all selected jets and standard MET

  double dpmin=20.;  // DP min MET / selected jets
  double dp1st=20.;  // DP MET/leading jet
  double dp2nd=20.;  // DP MET/2nd leading jet
  double dp3rd=20.;  // DP MET/3rd leading jet

  double m3best=-999.; // Massa del tripletto piu' vicino a mtref (200 GeV, see above)
  double mwbest=-999.;    // Massa del doppietto piu' vicino a 80.4 GeV tra i tre nel tripletto definito sopra
  double chi2=1000.;   // Chiquadro costruito con le due masse qui sopra
  double m45best=-999.;
  double mwmin=0.;
  double chi2ext=1000.;
  double m45bestall=-999.;
  double drpairbestall=3.99;
  double hbestcomb=0.;
  double chi2extall=1000.;
  double dp12 = 0.;        // Delta phi jet 1 jet 2
  double dpbb = 0.;    // delta phi bb
  double dpbball = 0.;    // delta phi bb
  double m_others = -999.; // mass of jets not part of best triplet and pair
  double mbbnohmax = -999.; // mass of b pair with highest mass excluding b-jets from h decay
  double dpbbnohmax = 0.; // angle of b pair with highest mass excluding b-jets from h decay
  double scprodthbest=0.; // scalar product between top and higgs
  double thdetabest=0.; // delta eta between top and higgs
  double m5=0.;         // mass of five jets assigned to top and higgs
  double ttms1=0.;      // total track mass of jets, S1 significance
  double ttms2=0.;
  double ttms3=0.;

  // Offline cuts
  // ------------
  bool NJetsCut = false;   // Whether event passes NJet cut
  bool MEtSigCut = false;  // Whether event pass MET significance cut

  double m6=0.;
  double c6=0.;
  double m8=0.;
  double c8=0.;
  double set=0.;
  double Et6=0.;

  // good jets array
  // ---------------
  int iJ=0;
  double JEt[100]={0.};
  double JEtc[100]={0.};
  double Jeta[100]={0.};
  double Jphi[100]={0.};
  double JHET[100]={0.};
  double JHPT[100]={0.};
  double JMT[100]={0.};
  int    JN1[100]={0};
  int    JBin[100]={0};
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

      if ( cal->et() > 25. && fabs( cal->eta() ) < 3.0 ) {
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
	  JMT[iJ]=cal->tagTkMassS1();
	  if ( JHET[iJ]>medium_ ) {
	    JHEM[iJ]=true;
	  }
	  // Correct jet Et
	  // --------------
	  if ( JHEM[iJ] ) {
	    JEtc[iJ]=JEt[iJ]+
	      tpar[0]*pow(JEt[iJ],6)+
	      tpar[1]*pow(JEt[iJ],5)+
	      tpar[2]*pow(JEt[iJ],4)+
	      tpar[3]*pow(JEt[iJ],3)+
	      tpar[4]*pow(JEt[iJ],2)+
	      tpar[5]*JEt[iJ]+
	      tpar[6];
	    if ( JEtc[iJ]<0. ) JEtc[iJ] = 0.;
	  } else {
	    JEtc[iJ]=JEt[iJ]+
	      upar[0]*pow(JEt[iJ],6)+
	      upar[1]*pow(JEt[iJ],5)+
	      upar[2]*pow(JEt[iJ],4)+
	      upar[3]*pow(JEt[iJ],3)+
	      upar[4]*pow(JEt[iJ],2)+
	      upar[5]*JEt[iJ]+
	      upar[6];
	    if ( JEtc[iJ]<0. ) JEtc[iJ] = 0.;
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
	  double a8 = JEtc[k];
	  double a9 = JMT[k];
	  JEt[k]=JEt[k-1];
	  Jeta[k]=Jeta[k-1];
	  Jphi[k]=Jphi[k-1];
	  JHET[k]=JHET[k-1];
	  JHPT[k]=JHPT[k-1];
	  JN1[k]=JN1[k-1];
	  JHEM[k]=JHEM[k-1];
	  JEtc[k]=JEtc[k-1];
	  JMT[k] = JMT[k-1];
	  JEt[k-1]=a1;
	  Jeta[k-1]=a2;
	  Jphi[k-1]=a3;
	  JHET[k-1]=a4;
	  JHPT[k-1]=a5;
	  JN1[k-1]=a6;
	  JHEM[k-1]=a7;
	  JEtc[k-1]=a8;
	  JMT[k-1]=a9;
	}
      }
    }
    for ( int ij=0; ij<iJ && ij<NHSJ; ij++ ) { 
      if ( JHEM[ij] ) NHEM++;   // NNNBBB Tags are counted only if in first 8 jets!
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

      // Assign tags to jets in QCD events according to their probability
      // ----------------------------------------------------------------
      int iJmax=NHSJ;  // Search for tags in the first NHSJ jets only!
      if ( iJ<iJmax ) iJmax=iJ;

      // Assign tag probability to each jet
      // by finding which bin it belongs to
      // ----------------------------------
      int etbin[8]={-1};
      int etabin[8]={-1};
      int ns1bin[8]={-1};
      for ( int i=0; i<iJmax; i++ ) {
	// Find out which bin this jet belongs to
	// --------------------------------------
	for ( int iet=0; iet<8; iet++ ) {
	  if ( JEt[i]>etbounds[iet] && JEt[i]<=etbounds[iet+1] ) {
	    etbin[i]=iet;
	  }
	}
	for ( int ieta=0; ieta<4; ieta++ ) {
	  if ( fabs(Jeta[i])>=etabounds[ieta] && fabs(Jeta[i])<etabounds[ieta+1] ) {
	    etabin[i]=ieta;
	  }
	}
	for ( int int1=0; int1<8; int1++ ) {
	  if ( JN1[i]>=ns1bounds[int1] && JN1[i]<ns1bounds[int1+1] ) {
	    ns1bin[i]=int1;
	  }
	}
	if ( etbin[i]>-1 && etabin[i]>-1 && ns1bin[i]>-1 ) {
	  JBin[i] = etbin[i]*100+etabin[i]*10+ns1bin[i];
	} else {
	  JBin[i] = 999;
	}
      }

      for ( int icomb=0; icomb<(int)pow(2.,iJmax); icomb++ ) { // Combinatorial loop ///////////////
	// Reset all variables
	// -------------------
	m3best=-999;  // Massa del tripletto piu' vicino a 172 GeV
	mwbest=-999;  // Massa del doppietto piu' vicino a 80.4 GeV 
	              // tra i tre nel tripletto definito sopra
	chi2=1000;    // Chiquadro costruito con le due masse qui sopra
	m45best=-999; // mass of two b-tags closest to 120 GeV not in triplet
        mwmin=0.;
	chi2ext=1000; // chi2 computed by adding contribution from distance to 120 GeV
	dp12 = 0.;
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
	c6=0.;
	c8=0.;
	Et6=0.;

	// Ok now start producing probability combinations
	// -----------------------------------------------
	double n=(double)icomb;
	double PTOT=1.;
	double PTOTE2=0.;
	NHEM=0;
	for ( int j=iJmax-1; j>=0; j-- ) {
	  if ( etbin[j]>-1 && etabin[j]>-1 && ns1bin[j]>-1 ) {
	    if ( n>=pow(2.,j) ) {
	      JHEM[j]=true;
	      n-=pow(2.,j);
	      PTOT=PTOT*PHETM[JBin[j]];
	      if ( PHETM[JBin[j]]>0 ) PTOTE2+=pow(PHETMS[JBin[j]]/PHETM[JBin[j]],2);
	      NHEM++;
	      // Extract random TM values from QCD pdf
	      // -------------------------------------
	      JMT[j]=MTS1pdf[ns1bin[j]]->GetRandom();
	      // Extract random HED and HPD values from QCD pdf
	      // ----------------------------------------------
	      double x=0.;
	      for ( ; x<medium_; ) {
		x=HEDpdf[ns1bin[j]]->GetRandom();
	      }   
	      JHET[j]=x;
	      if ( ns1bin[j]==0 ) {
		JHPT[j]=0.; // NNBB for 2-track jets, HPD=0
	      } else {
		x=0.;
		for ( ; x<tight_; ) { 
		  x=HPDpdf[ns1bin[j]]->GetRandom();
		}   
		JHPT[j]=x;
	      }
	    } else {
	      JHEM[j]=false;
	      PTOT=PTOT*(1.-PHETM[JBin[j]]);
	      if ( PHETM[JBin[j]]!=1. ) PTOTE2+=pow(PHETMS[JBin[j]]/(1.-PHETM[JBin[j]]),2);
	      // Extract random TM values from QCD pdf
	      // -------------------------------------
	      JMT[j]=MNS1pdf[ns1bin[j]]->GetRandom();
	      // Extract random HED and HPD values from QCD pdf
	      // ----------------------------------------------
	      double x=29.9;
	      for ( ; x>=medium_; ) {
		x=HEDpdf[ns1bin[j]]->GetRandom(); 
	      }   
	      JHET[j]=x;
	      JHPT[j]=x;
	      if ( ns1bin[j]==0 ) {
		JHPT[j]=0.;  // NNBB for 2-track jets, HPD=0
	      } else {
	      x=29.9;
	      for ( ; x>=tight_; ) {
		x=HPDpdf[ns1bin[j]]->GetRandom();
	      }   
	      }
	    }
	  } else {
	    if ( n>=pow(2.,j) ) n-=pow(2.,j);
	    JHEM[j]=false;
	    JMT[j]=0.;
	    JHET[j]=0.;
	    JHPT[j]=0.;
	  }
	}
	PTOTE2=PTOTE2*pow(PTOT,2);

	sumhed4=0.;
	sumhpd4=0.;
	sumhed6=0.;
	sumhpd6=0.;
	ttms1=0.;

	// Compute sum of discriminants
	// ----------------------------
	for ( int i=0; i<iJmax; i++ ) {
	  hemax[i]=JHET[i];
	  hpmax[i]=JHPT[i];
	}
	for ( int times=0; times<iJmax; times++ ) {
	  for ( int j=iJmax-1; j>=1; j-- ) {
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
	for ( int i=0; i<4 && i<iJmax; i++ ) {
	  sumhed4+=hemax[i];
	  sumhpd4+=hpmax[i];
	}
	for ( int i=0; i<6 && i<iJmax; i++ ) {
	  sumhed6+=hemax[i];
	  sumhpd6+=hpmax[i];
	}

	// Compute sum of tag masses
	// -------------------------
	for ( int i=0; i<iJmax; i++ ) {
	  ttms1+=JMT[i];
	}
          

	// ////////////////////////////////////////////////////////////////////////
	// Top and H kinematics
	// --------------------

	// Find best b-jet doublet
	// -----------------------
	double px4,px5,py4,py5,pz4,pz5,e4,e5;
	double m45=0.;
	double maxr=-1.;
	int ih1a=-1;
	int ih2a=-1;
	double ph2=0.;
	double phih=0.;
	double etah=0.;
	double pxh=0.;
	double pyh=0.;
	double pzh=0.;
	for ( int ii=0; ii<iJ-1 && ii<NHSJ-1; ii++ ) {
	  for ( int jj=ii+1; jj<iJ && jj<NHSJ; jj++ ) { // limit to NHSJ the number of jets from h.s.
	    
	    // Demand that both b-jets from H be tagged
	    // ----------------------------------------
	    if ( JHEM[ii] && JHEM[jj] ) {
	      px4=JEtc[ii]*cos(Jphi[ii]);
	      px5=JEtc[jj]*cos(Jphi[jj]);
	      py4=JEtc[ii]*sin(Jphi[ii]);
	      py5=JEtc[jj]*sin(Jphi[jj]);
	      pz4=JEtc[ii]/tan(2*atan(exp(-Jeta[ii])));
	      pz5=JEtc[jj]/tan(2*atan(exp(-Jeta[jj])));
	      e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
	      e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
	      m45=(e4+e5)*(e4+e5)-(px4+px5)*(px4+px5)-(py4+py5)*(py4+py5)-(pz4+pz5)*(pz4+pz5);
	      if ( m45>0 ) m45=sqrt(m45);
	      // Now define best bb combination based on H matrices
	      // --------------------------------------------------
	      double ptpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5));
	      double deta = Jeta[ii]-Jeta[jj];
	      double dphi = 3.1415926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
	      double drpair = sqrt(deta*deta+dphi*dphi);
	      double ptotpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5)+(pz4+pz5)*(pz4+pz5));
	      double etapair=0.;
	      if ( ptotpair-pz4-pz5!=0 ) etapair = 0.5*log((ptotpair+pz4+pz5)/(ptotpair-pz4-pz5));
	      double phipair = atan2(py4+py5,px4+px5);
	      int ipt=(int)((ptpair/500.)*10);
	      int ieta=(int)((fabs(etapair)/4.)*10);
	      int idr=(int)((drpair/4.)*10);
	      if ( ipt<0 ) ipt=0;
	      if ( ipt>9 ) ipt=9;
	      if ( ieta<0 ) ieta=0;
	      if ( ieta>9) ieta=9;
	      if ( idr<0 ) idr=0;
	      if ( idr>9 ) idr=9;
	      double r = Hread[100*idr+10*ieta+ipt];
	      if ( Hnotread[100*idr+10*ieta+ipt]>0 ) r=r/Hnotread[100*idr+10*ieta+ipt];
	      if ( r>maxr ) {
		maxr = r;
		hbestcomb=m45;
		ih1a=ii;
		ih2a=jj;
		ph2=ptotpair*ptotpair;
		phih=phipair;
		etah=etapair;
		pxh=px4+px5;
		pyh=py4+py5;
		pzh=pz4+pz5;
		dpbball=3.145926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
		drpairbestall = drpair;
	      }
	    }
	  }
	}
	double px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,e1,e2,e3;
	double px,py,pz;
	double spx,spy,spz,se;
	double m3a=0.;
	double mwa=0.;
	int it1a=-1;
	int it2a=-1;
	int it3a=-1;
	double tbestcomb=0.;
	double wbestcomb=0.;
	maxr=-1.;
	// First triplet with one b-jet
	// ----------------------------
	for ( int i=0; i<iJ-2 && i<NHSJ-2; i++ ) {
	  for ( int j=i+1; j<iJ-1 && j<NHSJ-1; j++ ) {
	    for ( int k=j+1; k<iJ && k<NHSJ; k++ ) { // limit to 8 the number of jets 
	                                                       // in hard subprocess
	      // Require that the triplet owns at least one b-tag
	      // ------------------------------------------------
	      int ib=-1;	      
	      if ( JHEM[i] ) ib=i;
	      if ( JHEM[j] ) ib=j;
	      if ( JHEM[k] ) ib=k;
	      if ( (JHEM[i]||JHEM[j]||JHEM[k]) && i!=ih1a && i!=ih2a && j!=ih1a && j!=ih2a 
		   && k!=ih1a && k!=ih2a ) {
		px1=JEtc[i]*cos(Jphi[i]);
		px2=JEtc[j]*cos(Jphi[j]);
		px3=JEtc[k]*cos(Jphi[k]);
		py1=JEtc[i]*sin(Jphi[i]);
		py2=JEtc[j]*sin(Jphi[j]);
		py3=JEtc[k]*sin(Jphi[k]);
		pz1=JEtc[i]/tan(2*atan(exp(-Jeta[i])));
		pz2=JEtc[j]/tan(2*atan(exp(-Jeta[j])));
		pz3=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
		e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
		e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
		e3 =sqrt(px3*px3+py3*py3+pz3*pz3);
		if ( m3a==0 ) {
		  m3a=(e1+e2+e3)*(e1+e2+e3)-(px1+px2+px3)*(px1+px2+px3)-
		    (py1+py2+py3)*(py1+py2+py3)-(pz1+pz2+pz3)*(pz1+pz2+pz3);
		  if ( m3a>0 ) m3a=sqrt(m3a);
		}
		double mw=0.;
		if ( ib==k ) {
		  mw=(e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2);
		} else if ( ib==j ) {
		  mw=(e1+e3)*(e1+e3)-(px1+px3)*(px1+px3)-(py1+py3)*(py1+py3)-(pz1+pz3)*(pz1+pz3);
		} else if ( ib==i ) {
		  mw=(e3+e2)*(e3+e2)-(px3+px2)*(px3+px2)-(py3+py2)*(py3+py2)-(pz3+pz2)*(pz3+pz2);
		}
		if ( mw>0 ) mw=sqrt(mw);
		if ( mwa==0 ) mwa==mw;
		double spx=px1+px2+px3;
		double spy=py1+py2+py3;
		double spz=pz1+pz2+pz3;
		double m123=pow(e1+e2+e3,2)-pow(spx,2)-pow(spy,2)-pow(spz,2);
		if ( m123>0 ) m123=sqrt(m123);		
		double pt3 = sqrt(spx*spx+spy*spy);
		double scprodth=0.;          // projection of top momentum in h direction
		if ( ph2>0 ) scprodth = (spx*pxh+spy*pyh+spz*pzh)/sqrt(ph2);
		double ptot3 = sqrt(spx*spx+spy*spy+spz*spz);
		double eta3=0.;
		if ( ptot3-spz!=0 ) eta3 = 0.5*log((ptot3+spz)/(ptot3-spz));
		double phi3 = atan2(spy,spx);
		double thdphi = 3.141592-fabs(fabs(phi3-phih)-3.141592);
		double thdeta = eta3-etah; // angle between momenta of t and h		
		int ipt=(int)((pt3/600.)*10);
		int idp=(int)((thdphi/3.15)*10);
		int iet=(int)((fabs(eta3)/4.)*10);
		if ( ipt<0 ) ipt=0;
		if ( ipt>9 ) ipt=9;
		if ( idp<0 ) idp=0;
		if ( idp>9 ) idp=9;
		if ( iet<0 ) iet=0;
		if ( iet>9 ) iet=9;
		double r = Tread[100*iet+10*idp+ipt];
		if ( Tnotread[100*iet+10*idp+ipt]>0 ) r=r/Tnotread[100*iet+10*idp+ipt];
		if ( r>maxr ) {
		  maxr = r;
		  tbestcomb=m123;
		  wbestcomb=mw;
		  scprodthbest=scprodth;
		  thdetabest=thdeta;
		  it1a=i;
		  it2a=j;
		  it3a=k;
		}

	      }
	    }
	  }
	}
	  
	// Compute mass of jets not selected in triplet of top and doublet of H decay
	// and mass of 5 jets chosen to be H and T decay products
	// --------------------------------------------------------------------------
	if ( ih1a+ih2a>0 && it1a+it2a+it3a>0 ) {
	  double spxo=0.;
	  double spyo=0.;
	  double spzo=0.;
	  double seo=0.;
	  double spx5=0.;
	  double spy5=0.;
	  double spz5=0.;
	  double se5=0.;
	  for ( int kk=0; kk<iJ && kk<NHSJ ; kk++ ) { // limit to NHSJ the number of 
	    // jets from hard subprocess
	    if ( kk!= ih1a && kk!=ih2a && kk!=it1a && kk!=it2a && kk!=it3a ) { 
	      px=JEtc[kk]*cos(Jphi[kk]);
	      py=JEtc[kk]*sin(Jphi[kk]);
	      pz=JEtc[kk]/tan(2*atan(exp(-Jeta[kk])));
	      spxo+=px;
	      spyo+=py;
	      spzo+=pz;
	      seo+=sqrt(px*px+py*py+pz*pz);
	    } else {
	      px=JEtc[kk]*cos(Jphi[kk]);
	      py=JEtc[kk]*sin(Jphi[kk]);
	      pz=JEtc[kk]/tan(2*atan(exp(-Jeta[kk])));
	      spx5+=px;
	      spy5+=py;
	      spz5+=pz;
	      se5+=sqrt(px*px+py*py+pz*pz);
	    }
	  }
	  m_others=seo*seo-spxo*spxo-spyo*spyo-spzo*spzo;
	  if ( m_others>0 ) m_others=sqrt(m_others);
	  m5=se5*se5-spx5*spx5-spy5*spy5-spz5*spz5;
	  if ( m5>0 ) m5=sqrt(m5);
	  
	  // Now compute mass and angle of two b-jets not chosen as H decay
	  // --------------------------------------------------------------
	  spx=0;
	  spy=0;
	  spz=0;
	  se=0;
	  double mbbnoh;
	  for ( int k1=0; k1<iJ-1 && k1<NHSJ-1; k1++ ) {
	    if ( k1!=ih1a && k1!=ih2a && JHEM[k1] ) {
	      for ( int k2=k1+1; k2<iJ && k2<NHSJ; k2++ ) { // limit to 8 the number of jets from h.s.
		if ( k2!=ih1a && k2!=ih2a && JHEM[k2] ) {
		  px1=JEtc[k1]*cos(Jphi[k1]);
		  py1=JEtc[k1]*sin(Jphi[k1]);
		  pz1=JEtc[k1]/tan(2*atan(exp(-Jeta[k1])));
		  px2=JEtc[k2]*cos(Jphi[k2]);
		  py2=JEtc[k2]*sin(Jphi[k2]);
		  pz2=JEtc[k2]/tan(2*atan(exp(-Jeta[k2])));
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

	  // Now compute global mass chisquared
	  // ----------------------------------
	  chi2 = pow((tbestcomb-mtref)/40.,2)+pow((wbestcomb-mwref)/20.,2)+pow((hbestcomb-mhref)/25.,2);


	}      // all 5 jets from h and t defined
	
	// Now compute mass of first 8 jets and their centrality
        // -----------------------------------------------------
	spx=0;
	spy=0;
	spz=0;
	se=0;
	set=0;
	for ( int k=0; k<iJ && k<NHSJ; k++ ) {
	  px=JEtc[k]*cos(Jphi[k]);
	  py=JEtc[k]*sin(Jphi[k]);
	  pz=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
	  spx+=px;
	  spy+=py;
	  spz+=pz;
	  se+=sqrt(px*px+py*py+pz*pz);
	  set+=JEtc[k];
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
	  px=JEtc[k]*cos(Jphi[k]);
	  py=JEtc[k]*sin(Jphi[k]);
	  pz=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
	  spx+=px;
	  spy+=py;
	  spz+=pz;
	  se+=sqrt(px*px+py*py+pz*pz);
	  set+=JEtc[k];
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
	
	if ( iJ>5 ) Et6=JEtc[5];
	
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
	    M3best_->Fill(tbestcomb,PTOT);
	    Mwbest_->Fill(wbestcomb,PTOT);
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
	    Et6_->Fill(Et6,PTOT);
	    Mwmin_->Fill(mwmin,PTOT);
	    Hbestcomb_->Fill(hbestcomb,PTOT);
	    Drpairbestall_->Fill(drpairbestall,PTOT);
	    M3a_->Fill(m3a,PTOT);
	    Mwa_->Fill(mwa,PTOT);
	    Scprod_->Fill(scprodthbest,PTOT);
	    Thdeta_->Fill(thdetabest,PTOT);
	    M5_->Fill(m5,PTOT);
	    TTMS1_->Fill(ttms1,PTOT);
	    TTMS2_->Fill(ttms2,PTOT);
	    TTMS3_->Fill(ttms3,PTOT);
    
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
	    M3bestN_->Fill(tbestcomb,weight_N);
	    MwbestN_->Fill(wbestcomb,weight_N);
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
	    Et6N_->Fill(Et6,weight_N);
	    MwminN_->Fill(mwmin,weight_N);
	    HbestcombN_->Fill(hbestcomb,weight_N);
	    DrpairbestallN_->Fill(drpairbestall,weight_N);
	    M3aN_->Fill(m3a,weight_N);
	    MwaN_->Fill(mwa,weight_N);
	    ScprodN_->Fill(scprodthbest,weight_N);
	    ThdetaN_->Fill(thdetabest,weight_N);
	    M5N_->Fill(m5,weight_N);
	    TTMS1N_->Fill(ttms1,weight_N);
	    TTMS2N_->Fill(ttms2,weight_N);
	    TTMS3N_->Fill(ttms3,weight_N);
	    
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
	      M3bestS_->Fill(tbestcomb,PTOT);
	      MwbestS_->Fill(wbestcomb,PTOT);
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
	      Et6S_->Fill(Et6,PTOT);
	      MwminS_->Fill(mwmin,PTOT);
	      HbestcombS_->Fill(hbestcomb,PTOT);
	      DrpairbestallS_->Fill(drpairbestall,PTOT);
	      M3aS_->Fill(m3a,PTOT);
	      MwaS_->Fill(mwa,PTOT);
	      ScprodS_->Fill(scprodthbest,PTOT);
	      ThdetaS_->Fill(thdetabest,PTOT);
	      M5S_->Fill(m5,PTOT);
	      TTMS1S_->Fill(ttms1,PTOT);
	      TTMS2S_->Fill(ttms2,PTOT);
	      TTMS3S_->Fill(ttms3,PTOT);
	      
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
	      M3bestSN_->Fill(tbestcomb,weight_N);
	      MwbestSN_->Fill(wbestcomb,weight_N);
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
	      Et6SN_->Fill(Et6,weight_N);
	      MwminSN_->Fill(mwmin,weight_N);
	      HbestcombSN_->Fill(hbestcomb,weight_N);
	      DrpairbestallSN_->Fill(drpairbestall,weight_N);
	      M3aSN_->Fill(m3a,weight_N);
	      MwaSN_->Fill(mwa,weight_N);
	      ScprodSN_->Fill(scprodthbest,weight_N);
	      ThdetaSN_->Fill(thdetabest,weight_N);
	      M5SN_->Fill(m5,weight_N);
	      TTMS1SN_->Fill(ttms1,weight_N);
	      TTMS2SN_->Fill(ttms2,weight_N);
	      TTMS3SN_->Fill(ttms3,weight_N);
	      
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
		M3bestSS_->Fill(tbestcomb,PTOT);
		MwbestSS_->Fill(wbestcomb,PTOT);
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
		Et6SS_->Fill(Et6,PTOT);
		MwminSS_->Fill(mwmin,PTOT);
		HbestcombSS_->Fill(hbestcomb,PTOT);
		DrpairbestallSS_->Fill(drpairbestall,PTOT);
		M3aSS_->Fill(m3a,PTOT);
		MwaSS_->Fill(mwa,PTOT);
		ScprodSS_->Fill(scprodthbest,PTOT);
		ThdetaSS_->Fill(thdetabest,PTOT);
		M5SS_->Fill(m5,PTOT);
		TTMS1SS_->Fill(ttms1,PTOT);
		TTMS2SS_->Fill(ttms2,PTOT);
		TTMS3SS_->Fill(ttms3,PTOT);
		
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
		M3bestSSN_->Fill(tbestcomb,weight_N);
		MwbestSSN_->Fill(wbestcomb,weight_N);
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
	      Et6SSN_->Fill(Et6,weight_N);
	      MwminSSN_->Fill(mwmin,weight_N);
	      HbestcombSSN_->Fill(hbestcomb,weight_N);
	      DrpairbestallSSN_->Fill(drpairbestall,weight_N);
	      M3aSSN_->Fill(m3a,weight_N);
	      MwaSSN_->Fill(mwa,weight_N);
	      ScprodSSN_->Fill(scprodthbest,weight_N);
	      ThdetaSSN_->Fill(thdetabest,weight_N);
	      M5SSN_->Fill(m5,weight_N);
	      TTMS1SSN_->Fill(ttms1,weight_N);
	      TTMS2SSN_->Fill(ttms2,weight_N);
	      TTMS3SSN_->Fill(ttms3,weight_N);
	      
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
		M3bestSSS_->Fill(tbestcomb,PTOT);
		MwbestSSS_->Fill(wbestcomb,PTOT);
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
		Et6SSS_->Fill(Et6,PTOT);
		MwminSSS_->Fill(mwmin,PTOT);
		HbestcombSSS_->Fill(hbestcomb,PTOT);
		DrpairbestallSSS_->Fill(drpairbestall,PTOT);
		M3aSSS_->Fill(m3a,PTOT);
		MwaSSS_->Fill(mwa,PTOT);
		ScprodSSS_->Fill(scprodthbest,PTOT);
		ThdetaSSS_->Fill(thdetabest,PTOT);
		M5SSS_->Fill(m5,PTOT);
		TTMS1SSS_->Fill(ttms1,PTOT);
		TTMS2SSS_->Fill(ttms2,PTOT);
		TTMS3SSS_->Fill(ttms3,PTOT);
		
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
		M3bestSSSN_->Fill(tbestcomb,weight_N);
		MwbestSSSN_->Fill(wbestcomb,weight_N);
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
		Et6SSSN_->Fill(Et6,weight_N);
		MwminSSSN_->Fill(mwmin,weight_N);
		HbestcombSSSN_->Fill(hbestcomb,weight_N);
		DrpairbestallSSSN_->Fill(drpairbestall,weight_N);
		M3aSSSN_->Fill(m3a,weight_N);
		MwaSSSN_->Fill(mwa,weight_N);
		ScprodSSSN_->Fill(scprodthbest,weight_N);
		ThdetaSSSN_->Fill(thdetabest,weight_N);
		M5SSSN_->Fill(m5,weight_N);
		TTMS1SSSN_->Fill(ttms1,weight_N);
		TTMS2SSSN_->Fill(ttms2,weight_N);
		TTMS3SSSN_->Fill(ttms3,weight_N);
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
	    M3bestW_->Fill(tbestcomb,PTOTE2);
	    MwbestW_->Fill(wbestcomb,PTOTE2);
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
	    Et6W_->Fill(Et6,PTOTE2);
	    MwminW_->Fill(mwmin,PTOTE2);
	    HbestcombW_->Fill(hbestcomb,PTOTE2);
	    DrpairbestallW_->Fill(drpairbestall,PTOTE2);
	    M3aW_->Fill(m3a,PTOTE2);
	    MwaW_->Fill(mwa,PTOTE2);
	    ScprodW_->Fill(scprodthbest,PTOTE2);
	    ThdetaW_->Fill(thdetabest,PTOTE2);
	    M5W_->Fill(m5,PTOTE2);
	    TTMS1W_->Fill(ttms1,PTOTE2);
	    TTMS2W_->Fill(ttms2,PTOTE2);
	    TTMS3W_->Fill(ttms3,PTOTE2);
	    
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
	      M3bestSW_->Fill(tbestcomb,PTOTE2);
	      MwbestSW_->Fill(wbestcomb,PTOTE2);
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
	      Et6SW_->Fill(Et6,PTOTE2);
	      MwminSW_->Fill(mwmin,PTOTE2);
	      HbestcombSW_->Fill(hbestcomb,PTOTE2);
	      DrpairbestallSW_->Fill(drpairbestall,PTOTE2);
	      M3aSW_->Fill(m3a,PTOTE2);
	      MwaSW_->Fill(mwa,PTOTE2);
	      ScprodSW_->Fill(scprodthbest,PTOTE2);
	      ThdetaSW_->Fill(thdetabest,PTOTE2);
	      M5SW_->Fill(m5,PTOTE2);
	      TTMS1SW_->Fill(ttms1,PTOTE2);
	      TTMS2SW_->Fill(ttms2,PTOTE2);
	      TTMS3SW_->Fill(ttms3,PTOTE2);
	      
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
		M3bestSSW_->Fill(tbestcomb,PTOTE2);
		MwbestSSW_->Fill(wbestcomb,PTOTE2);
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
		Et6SSW_->Fill(Et6,PTOTE2);
		MwminSSW_->Fill(mwmin,PTOTE2);
		HbestcombSSW_->Fill(hbestcomb,PTOTE2);
		DrpairbestallSSW_->Fill(drpairbestall,PTOTE2);
		M3aSSW_->Fill(m3a,PTOTE2);
		MwaSSW_->Fill(mwa,PTOTE2);
		ScprodSSW_->Fill(scprodthbest,PTOTE2);
		ThdetaSSW_->Fill(thdetabest,PTOTE2);
		M5SSW_->Fill(m5,PTOTE2);
		TTMS1SSW_->Fill(ttms1,PTOTE2);
		TTMS2SSW_->Fill(ttms2,PTOTE2);
		TTMS3SSW_->Fill(ttms3,PTOTE2);
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
		M3bestSSSW_->Fill(tbestcomb,PTOTE2);
		MwbestSSSW_->Fill(wbestcomb,PTOTE2);
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
		Et6SSSW_->Fill(Et6,PTOTE2);
		MwminSSSW_->Fill(mwmin,PTOTE2);
		HbestcombSSSW_->Fill(hbestcomb,PTOTE2);
		DrpairbestallSSSW_->Fill(drpairbestall,PTOTE2);
		M3aSSSW_->Fill(m3a,PTOTE2);
		MwaSSSW_->Fill(mwa,PTOTE2);
		ScprodSSSW_->Fill(scprodthbest,PTOTE2);
		ThdetaSSSW_->Fill(thdetabest,PTOTE2);
		M5SSSW_->Fill(m5,PTOTE2);
		TTMS1SSSW_->Fill(ttms1,PTOTE2);
		TTMS2SSSW_->Fill(ttms2,PTOTE2);
		TTMS3SSSW_->Fill(ttms3,PTOTE2);
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
	  bx[1]= (int)((m8/2500.)*50)+1;
	  bx[2]= (int)((c6/1.2)*50)+1;
	  bx[3]= (int)((GoodHt/4000.)*50)+1;
	  bx[4]= (int)((metsig/20.)*50.)+1;
	  bx[5]= (int)(((thdetabest-5.)/5.)*50)+1;
	  bx[6]= (int)((hbestcomb/400.)*50.)+1;
	  bx[7]= (int)((dp2nd/3.2)*50)+1;
	  bx[8]= (int)((m5/2000.)*50)+1;
	  bx[9]= (int)((m_others/1500.)*50.)+1;
	  bx[10]= (int)((Et6/150.)*50)+1;
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

      } // end for tag combinatorial
      
    } else { // if not QCD /////////////////////////////////////////////

      int iJmax=NHSJ;  // Search for tags in the first NHSJ jets only!
      if ( iJ<iJmax ) iJmax=iJ;

      sumhed4=0.;
      sumhpd4=0.;
      sumhed6=0.;
      sumhpd6=0.;
      ttms1=0.;

      // Compute sum of discriminants
      // ----------------------------
      for ( int i=0; i<iJmax; i++ ) {
	hemax[i]=JHET[i];
	hpmax[i]=JHPT[i];
      }
      for ( int times=0; times<iJmax; times++ ) {
	for ( int j=iJmax-1; j>=1; j-- ) {
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
      for ( int i=0; i<4 && i<iJmax; i++ ) {
	sumhed4+=hemax[i];
	sumhpd4+=hpmax[i];
      }
      for ( int i=0; i<6 && i<iJmax; i++ ) {
	sumhed6+=hemax[i];
	sumhpd6+=hpmax[i];
      }

      // Compute sum of tag masses
      // -------------------------
      for ( int i=0; i<iJmax; i++ ) {
	ttms1+=JMT[i];
      }

      // Fill HED and HPD histograms
      // ---------------------------
      for ( int i=0; i<iJmax; i++ ) {
	int ns1bin=-1;
	for ( int int1=0; int1<8; int1++ ) {
	  if ( JN1[i]>=ns1bounds[int1] && JN1[i]<ns1bounds[int1+1] ) {
	    ns1bin=int1;
	  }
	}
	if ( ns1bin==0 ) HED1_->Fill(JHET[i]);
	if ( ns1bin==0 ) HPD1_->Fill(JHPT[i]);
	if ( ns1bin==1 ) HED2_->Fill(JHET[i]);
	if ( ns1bin==1 ) HPD2_->Fill(JHPT[i]);
	if ( ns1bin==2 ) HED3_->Fill(JHET[i]);
	if ( ns1bin==2 ) HPD3_->Fill(JHPT[i]);
	if ( ns1bin==3 ) HED4_->Fill(JHET[i]);
	if ( ns1bin==3 ) HPD4_->Fill(JHPT[i]);
	if ( ns1bin==4 ) HED5_->Fill(JHET[i]);
	if ( ns1bin==4 ) HPD5_->Fill(JHPT[i]);
	if ( ns1bin==5 ) HED6_->Fill(JHET[i]);
	if ( ns1bin==5 ) HPD6_->Fill(JHPT[i]);
	if ( ns1bin==6 ) HED7_->Fill(JHET[i]);
	if ( ns1bin==6 ) HPD7_->Fill(JHPT[i]);
	if ( ns1bin==7 ) HED8_->Fill(JHET[i]);
	if ( ns1bin==7 ) HPD8_->Fill(JHPT[i]);
      }

      // ////////////////////////////////////////////////////////////////////////
      // Top and H kinematics
      // --------------------

      // Find best b-jet doublet
      // -----------------------
      double px4,px5,py4,py5,pz4,pz5,e4,e5;
      double m45=0.;
      double maxr=-1.;
      int ih1a=-1;
      int ih2a=-1;
      double ph2=0.;
      double phih=0.;
      double etah=0.;
      double pxh=0.;
      double pyh=0.;
      double pzh=0.;
      for ( int ii=0; ii<iJ-1 && ii<NHSJ-1; ii++ ) {
	for ( int jj=ii+1; jj<iJ && jj<NHSJ; jj++ ) { // limit to NHSJ the number of jets from h.s.
	    
	  // Demand that both b-jets from H be tagged
	  // ----------------------------------------
	  if ( JHEM[ii] && JHEM[jj] ) {
	    px4=JEtc[ii]*cos(Jphi[ii]);
	    px5=JEtc[jj]*cos(Jphi[jj]);
	    py4=JEtc[ii]*sin(Jphi[ii]);
	    py5=JEtc[jj]*sin(Jphi[jj]);
	    pz4=JEtc[ii]/tan(2*atan(exp(-Jeta[ii])));
	    pz5=JEtc[jj]/tan(2*atan(exp(-Jeta[jj])));
	    e4 =sqrt(px4*px4+py4*py4+pz4*pz4);
	    e5 =sqrt(px5*px5+py5*py5+pz5*pz5);
	    m45=(e4+e5)*(e4+e5)-(px4+px5)*(px4+px5)-(py4+py5)*(py4+py5)-(pz4+pz5)*(pz4+pz5);
	    if ( m45>0 ) m45=sqrt(m45);
	    // Now define best bb combination based on H matrices
	    // --------------------------------------------------
	    double ptpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5));
	    double deta = Jeta[ii]-Jeta[jj];
	    double dphi = 3.1415926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
	    double drpair = sqrt(deta*deta+dphi*dphi);
	    double ptotpair = sqrt((px4+px5)*(px4+px5)+(py4+py5)*(py4+py5)+(pz4+pz5)*(pz4+pz5));
	    double etapair=0.;
	    if ( ptotpair-pz4-pz5!=0 ) etapair = 0.5*log((ptotpair+pz4+pz5)/(ptotpair-pz4-pz5));
	    double phipair = atan2(py4+py5,px4+px5);
	    int ipt=(int)((ptpair/500.)*10);
	    int ieta=(int)((fabs(etapair)/4.)*10);
	    int idr=(int)((drpair/4.)*10);
	    if ( ipt<0 ) ipt=0;
	    if ( ipt>9 ) ipt=9;
	    if ( ieta<0 ) ieta=0;
	    if ( ieta>9) ieta=9;
	    if ( idr<0 ) idr=0;
	    if ( idr>9 ) idr=9;
	    double r = Hread[100*idr+10*ieta+ipt];
	    if ( Hnotread[100*idr+10*ieta+ipt]>0 ) r=r/Hnotread[100*idr+10*ieta+ipt];
	    if ( r>maxr ) {
	      maxr = r;
	      hbestcomb=m45;
	      ih1a=ii;
	      ih2a=jj;
	      ph2=ptotpair*ptotpair;
	      phih=phipair;
	      etah=etapair;
	      pxh=px4+px5;
	      pyh=py4+py5;
	      pzh=pz4+pz5;
	      dpbball=3.145926-fabs(fabs(Jphi[ii]-Jphi[jj])-3.1415926);
	      drpairbestall = drpair;
	    }
	  }
	}
      }
      double px1,px2,px3,py1,py2,py3,pz1,pz2,pz3,e1,e2,e3;
      double px,py,pz;
      double spx,spy,spz,se;
      double m3a=0.;
      double mwa=0.;
      int it1a=-1;
      int it2a=-1;
      int it3a=-1;
      double tbestcomb=0.;
      double wbestcomb=0.;
      maxr=-1.;
      // First triplet with one b-jet
      // ----------------------------
      for ( int i=0; i<iJ-2 && i<NHSJ-2; i++ ) {
	for ( int j=i+1; j<iJ-1 && j<NHSJ-1; j++ ) {
	  for ( int k=j+1; k<iJ && k<NHSJ; k++ ) { // limit to 8 the number of jets 
	    // in hard subprocess
	    // Require that the triplet owns at least one b-tag
	    // ------------------------------------------------
	    int ib=-1;	      
	    if ( JHEM[i] ) ib=i;
	    if ( JHEM[j] ) ib=j;
	    if ( JHEM[k] ) ib=k;
	    if ( (JHEM[i]||JHEM[j]||JHEM[k]) && i!=ih1a && i!=ih2a && j!=ih1a && j!=ih2a 
		 && k!=ih1a && k!=ih2a ) {
	      px1=JEtc[i]*cos(Jphi[i]);
	      px2=JEtc[j]*cos(Jphi[j]);
	      px3=JEtc[k]*cos(Jphi[k]);
	      py1=JEtc[i]*sin(Jphi[i]);
	      py2=JEtc[j]*sin(Jphi[j]);
	      py3=JEtc[k]*sin(Jphi[k]);
	      pz1=JEtc[i]/tan(2*atan(exp(-Jeta[i])));
	      pz2=JEtc[j]/tan(2*atan(exp(-Jeta[j])));
	      pz3=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
	      e1 =sqrt(px1*px1+py1*py1+pz1*pz1);
	      e2 =sqrt(px2*px2+py2*py2+pz2*pz2);
	      e3 =sqrt(px3*px3+py3*py3+pz3*pz3);
	      if ( m3a==0 ) {
		m3a=(e1+e2+e3)*(e1+e2+e3)-(px1+px2+px3)*(px1+px2+px3)-
		  (py1+py2+py3)*(py1+py2+py3)-(pz1+pz2+pz3)*(pz1+pz2+pz3);
		if ( m3a>0 ) m3a=sqrt(m3a);
	      }
	      double mw=0.;
	      if ( ib==k ) {
		mw=(e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2);
	      } else if ( ib==j ) {
		mw=(e1+e3)*(e1+e3)-(px1+px3)*(px1+px3)-(py1+py3)*(py1+py3)-(pz1+pz3)*(pz1+pz3);
	      } else if ( ib==i ) {
		mw=(e3+e2)*(e3+e2)-(px3+px2)*(px3+px2)-(py3+py2)*(py3+py2)-(pz3+pz2)*(pz3+pz2);
	      }
	      if ( mw>0 ) mw=sqrt(mw);
	      if ( mwa==0 ) mwa==mw;
	      double spx=px1+px2+px3;
	      double spy=py1+py2+py3;
	      double spz=pz1+pz2+pz3;
	      double m123=pow(e1+e2+e3,2)-pow(spx,2)-pow(spy,2)-pow(spz,2);
	      if ( m123>0 ) m123=sqrt(m123);		
	      double pt3 = sqrt(spx*spx+spy*spy);
	      double scprodth=0.;          // projection of top momentum in h direction
	      if ( ph2>0 ) scprodth = (spx*pxh+spy*pyh+spz*pzh)/sqrt(ph2);
	      double ptot3 = sqrt(spx*spx+spy*spy+spz*spz);
	      double eta3=0.;
	      if ( ptot3-spz!=0 ) eta3 = 0.5*log((ptot3+spz)/(ptot3-spz));
	      double phi3 = atan2(spy,spx);
	      double thdphi = 3.141592-fabs(fabs(phi3-phih)-3.141592);
	      double thdeta = eta3-etah; // angle between momenta of t and h		
	      int ipt=(int)((pt3/600.)*10);
	      int idp=(int)((thdphi/3.15)*10);
	      int iet=(int)((fabs(eta3)/4.)*10);
	      if ( ipt<0 ) ipt=0;
	      if ( ipt>9 ) ipt=9;
	      if ( idp<0 ) idp=0;
	      if ( idp>9 ) idp=9;
	      if ( iet<0 ) iet=0;
	      if ( iet>9 ) iet=9;
	      double r = Tread[100*iet+10*idp+ipt];
	      if ( Tnotread[100*iet+10*idp+ipt]>0 ) r=r/Tnotread[100*iet+10*idp+ipt];
	      if ( r>maxr ) {
		maxr = r;
		tbestcomb=m123;
		wbestcomb=mw;
		scprodthbest=scprodth;
		thdetabest=thdeta;
		it1a=i;
		it2a=j;
		it3a=k;
	      }

	    }
	  }
	}
      }
	  
      // Compute mass of jets not selected in triplet of top and doublet of H decay
      // and mass of 5 jets chosen to be H and T decay products
      // --------------------------------------------------------------------------
      if ( ih1a+ih2a>0 && it1a+it2a+it3a>0 ) {
	double spxo=0.;
	double spyo=0.;
	double spzo=0.;
	double seo=0.;
	double spx5=0.;
	double spy5=0.;
	double spz5=0.;
	double se5=0.;
	for ( int kk=0; kk<iJ && kk<NHSJ ; kk++ ) { // limit to NHSJ the number of 
	  // jets from hard subprocess
	  if ( kk!= ih1a && kk!=ih2a && kk!=it1a && kk!=it2a && kk!=it3a ) { 
	    px=JEtc[kk]*cos(Jphi[kk]);
	    py=JEtc[kk]*sin(Jphi[kk]);
	    pz=JEtc[kk]/tan(2*atan(exp(-Jeta[kk])));
	    spxo+=px;
	    spyo+=py;
	    spzo+=pz;
	    seo+=sqrt(px*px+py*py+pz*pz);
	  } else {
	    px=JEtc[kk]*cos(Jphi[kk]);
	    py=JEtc[kk]*sin(Jphi[kk]);
	    pz=JEtc[kk]/tan(2*atan(exp(-Jeta[kk])));
	    spx5+=px;
	    spy5+=py;
	    spz5+=pz;
	    se5+=sqrt(px*px+py*py+pz*pz);
	  }
	}
	m_others=seo*seo-spxo*spxo-spyo*spyo-spzo*spzo;
	if ( m_others>0 ) m_others=sqrt(m_others);
	m5=se5*se5-spx5*spx5-spy5*spy5-spz5*spz5;
	if ( m5>0 ) m5=sqrt(m5);
	  
	// Now compute mass and angle of two b-jets not chosen as H decay
	// --------------------------------------------------------------
	spx=0;
	spy=0;
	spz=0;
	se=0;
	double mbbnoh;
	for ( int k1=0; k1<iJ-1 && k1<NHSJ-1; k1++ ) {
	  if ( k1!=ih1a && k1!=ih2a && JHEM[k1] ) {
	    for ( int k2=k1+1; k2<iJ && k2<NHSJ; k2++ ) { // limit to 8 the number of jets from h.s.
	      if ( k2!=ih1a && k2!=ih2a && JHEM[k2] ) {
		px1=JEtc[k1]*cos(Jphi[k1]);
		py1=JEtc[k1]*sin(Jphi[k1]);
		pz1=JEtc[k1]/tan(2*atan(exp(-Jeta[k1])));
		px2=JEtc[k2]*cos(Jphi[k2]);
		py2=JEtc[k2]*sin(Jphi[k2]);
		pz2=JEtc[k2]/tan(2*atan(exp(-Jeta[k2])));
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

	// Now compute global mass chisquared
	// ----------------------------------
	chi2 = pow((tbestcomb-mtref)/40.,2)+pow((wbestcomb-mwref)/20.,2)+pow((hbestcomb-mhref)/25.,2);

      }      // all 5 jets from h and t defined
	
      // Now compute mass of first 8 jets and their centrality
      // -----------------------------------------------------
      spx=0;
      spy=0;
      spz=0;
      se=0;
      set=0;
      for ( int k=0; k<iJ && k<NHSJ; k++ ) {
	px=JEtc[k]*cos(Jphi[k]);
	py=JEtc[k]*sin(Jphi[k]);
	pz=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
	spx+=px;
	spy+=py;
	spz+=pz;
	se+=sqrt(px*px+py*py+pz*pz);
	set+=JEtc[k];
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
	px=JEtc[k]*cos(Jphi[k]);
	py=JEtc[k]*sin(Jphi[k]);
	pz=JEtc[k]/tan(2*atan(exp(-Jeta[k])));
	spx+=px;
	spy+=py;
	spz+=pz;
	se+=sqrt(px*px+py*py+pz*pz);
	set+=JEtc[k];
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
      
      if ( iJ>5 ) Et6=JEtc[5];

      double UncorrMEtSig=0;
      if ( UncorrSumEt>0 ) UncorrMEtSig = UncorrMEt/sqrt(UncorrSumEt);
      double CorrMEtSig = 0;
      if ( CorrSumEt>0 ) CorrMEtSig = CorrMEt/sqrt(CorrSumEt);  
      
      // Missing Et significance computed with tuned resolution function
      // ---------------------------------------------------------------
      double metsignew=0;
      if ( sumet>0 ) 
	metsignew = met/(sqrt(2)*0.8033*pow(sumet,0.5004));
      
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
	M3best_->Fill(tbestcomb);
	Mwbest_->Fill(wbestcomb);
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
	Et6_->Fill(Et6);
	Mwmin_->Fill(mwmin);
	Hbestcomb_->Fill(hbestcomb);
	Drpairbestall_->Fill(drpairbestall);
	M3a_->Fill(m3a);
	Mwa_->Fill(mwa);
	Scprod_->Fill(scprodthbest);
	Thdeta_->Fill(thdetabest);
	M5_->Fill(m5);
	TTMS1_->Fill(ttms1);
	TTMS2_->Fill(ttms2);
	TTMS3_->Fill(ttms3);
	
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
	  M3bestS_->Fill(tbestcomb);
	  MwbestS_->Fill(wbestcomb);
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
	  Et6S_->Fill(Et6);
	  MwminS_->Fill(mwmin);
	  HbestcombS_->Fill(hbestcomb);
	  DrpairbestallS_->Fill(drpairbestall);
	  M3aS_->Fill(m3a);
	  MwaS_->Fill(mwa);
	  ScprodS_->Fill(scprodthbest);
	  ThdetaS_->Fill(thdetabest);
	  M5S_->Fill(m5);
	  TTMS1S_->Fill(ttms1);
	  TTMS2S_->Fill(ttms2);
	  TTMS3S_->Fill(ttms3);
	  
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
	    M3bestSS_->Fill(tbestcomb);
	    MwbestSS_->Fill(wbestcomb);
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
	    Et6SS_->Fill(Et6);
	    MwminSS_->Fill(mwmin);
	    HbestcombSS_->Fill(hbestcomb);
	    DrpairbestallSS_->Fill(drpairbestall);
	    M3aSS_->Fill(m3a);
	    MwaSS_->Fill(mwa);
	    ScprodSS_->Fill(scprodthbest);
	    ThdetaSS_->Fill(thdetabest);
	    M5SS_->Fill(m5);
	    TTMS1SS_->Fill(ttms1);
	    TTMS2SS_->Fill(ttms2);
	    TTMS3SS_->Fill(ttms3);
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
	    M3bestSSS_->Fill(tbestcomb);
	    MwbestSSS_->Fill(wbestcomb);
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
	    Et6SSS_->Fill(Et6);
	    MwminSSS_->Fill(mwmin);
	    HbestcombSSS_->Fill(hbestcomb);
	    DrpairbestallSSS_->Fill(drpairbestall);
	    M3aSSS_->Fill(m3a);
	    MwaSSS_->Fill(mwa);
	    ScprodSSS_->Fill(scprodthbest);
	    ThdetaSSS_->Fill(thdetabest);
	    M5SSS_->Fill(m5);
	    TTMS1SSS_->Fill(ttms1);
	    TTMS2SSS_->Fill(ttms2);
	    TTMS3SSS_->Fill(ttms3);

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
//       bx[0]= (int)((c8/1.2)*50)+1;
//       bx[1]= (int)((m45bestall/300.)*50)+1;
//       bx[2]= (int)(chi2extall)+1;
//       bx[3]= (int)((metsig/20.)*50)+1;
//       bx[4]= (int)((m8/2000)*50.)+1;
//       bx[5]= (int)((met/500.)*50)+1;
//       bx[6]= (int)((dp12/3.2)*50.)+1;
//       bx[7]= (int)((dp2nd/3.2)*50)+1;
//       bx[8]= (int)((c6/1.2)*50)+1;
//       bx[9]= (int)((m6/2500.)*50.)+1;
//       bx[10]= (int)((sumet/4000.)*50)+1;
      bx[0]= (int)((c8/1.2)*50)+1;
      bx[1]= (int)((m8/2500.)*50)+1;
      bx[2]= (int)((c6/1.2)*50)+1;
      bx[3]= (int)((GoodHt/4000.)*50)+1;
      bx[4]= (int)((metsig/20.)*50.)+1;
      bx[5]= (int)(((thdetabest-5.)/5.)*50)+1;
      bx[6]= (int)((hbestcomb/400.)*50.)+1;
      bx[7]= (int)((dp2nd/3.2)*50)+1;
      bx[8]= (int)((m5/2000.)*50)+1;
      bx[9]= (int)((m_others/1500.)*50.)+1;
      bx[10]= (int)((Et6/150.)*50)+1;
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

    } // if QCD or not QCD

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

  // Write H matrix containing kinematics of b-pair
  // ----------------------------------------------
  ofstream Hfile("H.asc");
  for ( int i=0; i<1000; i++ ) {
    Hfile << H[i] << " " << Hnot[i] << endl;
  }
  Hfile.close();

  // Write T matrix containing kinematics of triplet
  // -----------------------------------------------
  ofstream Tfile("T.asc");
  for ( int i=0; i<1000; i++ ) {
    Tfile << T[i] << " " << Tnot[i] << endl;
  }
  Tfile.close();

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
    cout << Hname[hdecay] << setprecision(7) << grandtotalttpass[hdecay] << "/" 
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
  cout << setprecision(7) << grandgrandtotalpass << "/" << grandgrandtotal << " = " 
       << "(" <<setprecision(5) << f0t*100. 
       << "+-" << setprecision(5) << sf0t*100. << ") %" << endl;

  decayfile << "Totals:     " 
	    << setprecision(6) << setw(7) << grandtotalttpass[1] << "/" 
	    <<  setprecision(7) << grandtotaltt[1] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[2] << "/" 
	    <<  setprecision(7) << grandtotaltt[2] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[3] << "/" 
	    <<  setprecision(7) << grandtotaltt[3] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[4] << "/" 
	    <<  setprecision(7) << grandtotaltt[4] << " " 
	    << setprecision(6) << setw(8) << grandtotalttpass[0] << "/" 
	    <<  setprecision(7) << grandtotaltt[0] << " " 
	    << setprecision(6) << setw(8) << grandgrandtotalpass << "/" 
	    <<  setprecision(7) << grandgrandtotal << endl;
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
  cout << "Writing histograms..." << endl;
  if ( Nsltt_hjj>0 ) {
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
    Trecfrac_->Scale(1./Nsltt_hjj);

  } else cout << "N sl tt decays=0" << endl;

  DEtb_prof_->Write();
  DEtq_prof_->Write();
  DEtcb_prof_->Write();
  DEtcq_prof_->Write();

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
  Trecfrac_->Write();

  MHbest_->Write();
  MTbest_->Write();
  MWbest_->Write();
  HBJ_etrank_->Write();
  Hpt_->Write();
  Heta_->Write();
  Hdr_->Write();
  MHnot_->Write();
  Hnotpt_->Write();
  Hnoteta_->Write();
  Hnotdr_->Write();
  MTnotbest_->Write();
  Tpt_->Write();
  Teta_->Write();
  THdeta_->Write();
  THdphi_->Write();
  THproj_->Write();
  Tnotpt_->Write();
  Tnoteta_->Write();
  THnotdphi_->Write();
  THnotdeta_->Write();
  THnotproj_->Write();

  HED1_->Write();
  HPD1_->Write();
  HED2_->Write();
  HPD2_->Write();
  HED3_->Write();
  HPD3_->Write();
  HED4_->Write();
  HPD4_->Write();
  HED5_->Write();
  HPD5_->Write();
  HED6_->Write();
  HPD6_->Write();
  HED7_->Write();
  HPD7_->Write();
  HED8_->Write();
  HPD8_->Write();

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
  Et6_->Write(); 
  Mwmin_->Write();
  Hbestcomb_->Write();
  Drpairbestall_->Write();
  M3a_->Write();
  Mwa_->Write();
  Scprod_->Write();
  Thdeta_->Write();
  M5_->Write();
  TTMS1_->Write();
  TTMS2_->Write();
  TTMS3_->Write();

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
  Et6N_->Write(); 
  MwminN_->Write();
  HbestcombN_->Write();
  DrpairbestallN_->Write();
  M3aN_->Write();
  MwaN_->Write();
  ScprodN_->Write();
  ThdetaN_->Write();
  M5N_->Write();
  TTMS1N_->Write();
  TTMS2N_->Write();
  TTMS3N_->Write();

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
  Et6S_->Write(); 
  MwminS_->Write();
  HbestcombS_->Write();
  DrpairbestallS_->Write();
  M3aS_->Write();
  MwaS_->Write();
  ScprodS_->Write();
  ThdetaS_->Write();
  M5S_->Write();
  TTMS1S_->Write();
  TTMS2S_->Write();
  TTMS3S_->Write();

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
  Et6SN_->Write(); 
  MwminSN_->Write();
  HbestcombSN_->Write();
  DrpairbestallSN_->Write();
  M3aSN_->Write();
  MwaSN_->Write();
  ScprodSN_->Write();
  ThdetaSN_->Write();
  M5SN_->Write();
  TTMS1SN_->Write();
  TTMS2SN_->Write();
  TTMS3SN_->Write();

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
  Et6SS_->Write(); 
  MwminSS_->Write();
  HbestcombSS_->Write();
  DrpairbestallSS_->Write();
  M3aSS_->Write();
  MwaSS_->Write();
  ScprodSS_->Write();
  ThdetaSS_->Write();
  M5SS_->Write();
  TTMS1SS_->Write();
  TTMS2SS_->Write();
  TTMS3SS_->Write();

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
  Et6SSN_->Write(); 
  MwminSSN_->Write();
  HbestcombSSN_->Write();
  DrpairbestallSSN_->Write();
  M3aSSN_->Write();
  MwaSSN_->Write();
  ScprodSSN_->Write();
  ThdetaSSN_->Write();
  M5SSN_->Write();
  TTMS1SSN_->Write();
  TTMS2SSN_->Write();
  TTMS3SSN_->Write();

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
  Et6SSS_->Write(); 
  MwminSSS_->Write();
  HbestcombSSS_->Write();
  DrpairbestallSSS_->Write();
  M3aSSS_->Write();
  MwaSSS_->Write();
  ScprodSSS_->Write();
  ThdetaSSS_->Write();
  M5SSS_->Write();
  TTMS1SSS_->Write();
  TTMS2SSS_->Write();
  TTMS3SSS_->Write();

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
  Et6SSSN_->Write(); 
  MwminSSSN_->Write();
  HbestcombSSSN_->Write();
  DrpairbestallSSSN_->Write();
  M3aSSSN_->Write();
  MwaSSSN_->Write();
  ScprodSSSN_->Write();
  ThdetaSSSN_->Write();
  M5SSSN_->Write();
  TTMS1SSSN_->Write();
  TTMS2SSSN_->Write();
  TTMS3SSSN_->Write();

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
  Et6W_->Write(); 
  MwminW_->Write();
  HbestcombW_->Write();
  DrpairbestallW_->Write();
  M3aW_->Write();
  MwaW_->Write();
  ScprodW_->Write();
  ThdetaW_->Write();
  M5W_->Write();
  TTMS1W_->Write();
  TTMS2W_->Write();
  TTMS3W_->Write();

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
  Et6SW_->Write(); 
  MwminSW_->Write();
  HbestcombSW_->Write();
  DrpairbestallSW_->Write();
  M3aSW_->Write();
  MwaSW_->Write();
  ScprodSW_->Write();
  ThdetaSW_->Write();
  M5SW_->Write();
  TTMS1SW_->Write();
  TTMS2SW_->Write();
  TTMS3SW_->Write();

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
  Et6SSW_->Write(); 
  MwminSSW_->Write();
  HbestcombSSW_->Write();
  DrpairbestallSSW_->Write();
  M3aSSW_->Write();
  MwaSSW_->Write();
  ScprodSSW_->Write();
  ThdetaSSW_->Write();
  M5SSW_->Write();
  TTMS1SSW_->Write();
  TTMS2SSW_->Write();
  TTMS3SSW_->Write();

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
  Et6SSSW_->Write(); 
  MwminSSSW_->Write();
  HbestcombSSSW_->Write();
  DrpairbestallSSSW_->Write();
  M3aSSSW_->Write();
  MwaSSSW_->Write();
  ScprodSSSW_->Write();
  ThdetaSSSW_->Write();
  M5SSSW_->Write();
  TTMS1SSSW_->Write();
  TTMS2SSSW_->Write();
  TTMS3SSSW_->Write();

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
