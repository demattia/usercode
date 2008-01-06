// -*- C++ -*-
//
// Package:    JetResMaker
// Class:      JetResMaker
// 
/**\class JetResMaker JetResMaker.cc AnalysisExamples/JetResMaker/src/JetResMaker.cc
 *
 * Description:
 * Creates the files with the OfflineJet vs GenJet resolutions.
 *
 * Implementation:
 *
 * Evaluates:
 *
 */
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: JetResMaker.h,v 1.1 2008/01/06 10:44:36 demattia Exp $
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

class JetResMaker : public edm::EDAnalyzer {
public:
  explicit JetResMaker(const edm::ParameterSet&);
  ~JetResMaker();

  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one JetResMaker object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Histos
  // ------
  // The histograms must be created after the TFile is opened.
  // ---------------------------------------------------------

  TH1F * Ptres[40];
  TH1F * Dr[40];
  TH1F * Ave_vs_pt[4];
  TH1F * Res_vs_pt[4];
  TH1F * l1Ptres[40];
  TH1F * l1Dr[40];

  // Debug histos
  // ------------
  TH1F * Ipt;
  TH1F * Ieta;
  TH1F * Totptleft;
  TH1F * Totptass;
  TH1F * Nleft;
  TH1F * Nass;
  TH1F * Njtot;
  TH1F * Njgood;
  TH1F * l1Ipt;
  TH1F * l1Ieta;
  TH1F * l1Totptleft;
  TH1F * l1Totptass;
  TH1F * l1Nleft;
  TH1F * l1Nass;
  TH1F * l1Njtot;
  TH1F * l1Njgood;

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
