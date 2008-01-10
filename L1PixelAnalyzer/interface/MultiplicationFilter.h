// -*- C++ -*-
//
// Package:    MultiplicationFilter
// Class:      MultiplicationFilter
// 
/**\class MultiplicationFilter MultiplicationFilter.cc AnalysisExamples/MultiplicationFilter/src/MultiplicationFilter.cc

 Description: <one line class summary>

 Implementation:

 Evaluates:

*/
//
// Original Author:  Marco De Mattia
//         Created:  Wed Gen 9 16:00:00 CEST 2007
// $Id: MultiplicationFilter.h,v 1.2 2007/11/12 18:01:49 demattia Exp $
//
//

// System include files
// --------------------
#include <memory>
#include <vector>

// User include files
// ------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

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

//
// class declaration
//

class MultiplicationFilter : public edm::EDFilter {
 public:
  explicit MultiplicationFilter(const edm::ParameterSet&);
  ~MultiplicationFilter();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  int eventCounter_;

  edm::ParameterSet conf_;

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
  // Seed for the random number generator for the histograms
  int seed_;
  string outputFileName_;
  // Total number of events to generate (includes the original event as the first one)
  int nTotChanged_;
  int nTotWritten_;
  int nTotGoodDiscarded_;

  int eventType_;
};
