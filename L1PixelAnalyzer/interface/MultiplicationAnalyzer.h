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
// $Id: MultiplicationAnalyzer.h,v 1.2 2008/01/06 11:43:48 demattia Exp $
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
};
