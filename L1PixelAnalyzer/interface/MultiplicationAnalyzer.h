// -*- C++ -*-
//
// Package:    MultiplicationAnalyzer
// Class:      MultiplicationAnalyzer
// 
/**\class MultiplicationAnalyzer MultiplicationAnalyzer.cc AnalysisExamples/MultiplicationAnalyzer/src/MultiplicationAnalyzer.cc
 *
 * Description:
 * Analyzer for the multiplication.
 *
 * Implementation:
 *
 * Evaluates:
 *
 */
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: MultiplicationAnalyzer.h,v 1.1 2008/01/06 18:43:47 demattia Exp $
//
//

// System include files
// --------------------
#include <memory>
#include <vector>

// User include files
// ------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

// Root includes
// -------------
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLegend.h"

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"
// EtSorter
#include "AnalysisExamples/AnalysisClasses/interface/EtSort.h"

// L1Trigger evaluator
// -------------------
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiStack.h"
#include "AnalysisExamples/AnalysisClasses/interface/L1PixelTrig.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Multiplier
#include "AnalysisExamples/L1PixelAnalyzer/interface/Multiplier.h"

// Class declaration
// -----------------

class MultiplicationAnalyzer : public edm::EDAnalyzer {
public:
  explicit MultiplicationAnalyzer(const edm::ParameterSet&);
  ~MultiplicationAnalyzer();

  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one MultiplicationAnalyzer object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Histos
  // ------
  // The histograms must be created after the TFile is opened.
  // ---------------------------------------------------------
  // Histograms
  // ----------
  TH1F * NJET;
  TH1F * JETET;
  TH1F * JETET_5j;
  TH1F * JETET_5j_s3;
  TH1F * MET;
  TH1F * MET_5j;
  TH1F * MET_5j_s3;
  TH1F * METX;
  TH1F * METX_5j;
  TH1F * METX_5j_s3;
  TH1F * SET;
  TH1F * SET_5j;
  TH1F * SET_5j_s3;
  TH1F * SMET;
  TH1F * SMET_5j;
  TH1F * SMET_5j_s3;

  TH1F * NJETC;
  TH1F * JETETC;
  TH1F * JETETC_5j;
  TH1F * JETETC_5j_s3;
  TH1F * METC;
  TH1F * METC_5j;
  TH1F * METC_5j_s3;
  TH1F * METXC;
  TH1F * METXC_5j;
  TH1F * METXC_5j_s3;
  TH1F * SETC;
  TH1F * SETC_5j;
  TH1F * SETC_5j_s3;
  TH1F * SMETC;
  TH1F * SMETC_5j;
  TH1F * SMETC_5j_s3;

  TH1F * DPtmod;
  TH2F * DSET_DMET;
  TH1F * DSET;
  TH1F * DMET;

  TH1F * ChangedNum;
  TH1F * ChangedRatio;

  // For Kolmogorov test
  TH1F * KolmogorovTestMEt;
  TH1F * KolMEt[50];
  TH1F * KolmogorovTestSumEt;
  TH1F * KolSumEt[50];

  // Collections labels
  // ------------------
  edm::InputTag cenJetLabel_;
  edm::InputTag forJetLabel_;
  edm::InputTag tauJetLabel_;
  edm::InputTag l1MEtLabel_;
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag MCParticleLabel_;
  edm::InputTag simplePixelJetLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag summaryLabel_;
  edm::InputTag genJetLabel_;

  // For multiplication
  Multiplier * multiplier;
  double minMultiplicationEt_;
  double mEtAlpha_;
  // Total number of events to generate (includes the original event as the first one)
  int Ntot_;
};
