// -*- C++ -*-
//
// Package:    OfflineAnalyzer
// Class:      OfflineAnalyzer
// 
/**\class OfflineAnalyzer OfflineAnalyzer.cc AnalysisExamples/OfflineAnalyzer/src/OfflineAnalyzer.cc

 Description: <one line class summary>

 Implementation:
 This class shows how to access:
 - level 1 calorimetric quantities
 - offline corrected jets (calibration performed here)
 - offline corrected MET, depending on jets corrections
 - MC informations         <---------------------------------- to do
 - B tagging               <---------------------------------- to do

 Evaluates:
 - DPhimin between MET and closest (in phi) offline jet
 - association of MC partons to offline jets            <----- to do
 - association of btags to offline jets                 <----- to do

*/
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: OfflineAnalyzer.h,v 1.5 2008/01/14 13:27:00 demattia Exp $
//
//

// system include files
#include <memory>
#include <vector>

// user include files
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
#include "TProfile.h"
#include "TFile.h"

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"

// L1Trigger evaluator
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"

//
// class declaration
//

class OfflineAnalyzer : public edm::EDAnalyzer {
 public:
  explicit OfflineAnalyzer(const edm::ParameterSet&);
  ~OfflineAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one OfflineAnalyzer object is created
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the OfflineAnalyzer constructor ).
  //  HiVariables HiVar;

  // Use a dynamic construction, or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.
  HiVariables * HiVar;

  edm::InputTag cenJetLabel_;
  edm::InputTag forJetLabel_;
  edm::InputTag tauJetLabel_;
  edm::InputTag l1MEtLabel_;
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag MCParticleLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag summaryLabel_;

  unsigned int numTkCut_;
  double minDz_;
  double maxDz_;
  bool doTrigger_;
  bool extendedInfo_;
  std::string OutputEffFileName;

  TH1F * uncorr_JetPt_IC5_;
  TH1F * corr_JetPt_IC5_;
  TH1F * JetNumber_IC5_;

  int corrTypeNum_;
  vector<TH1F *> MEtPt_;
  vector<TH1F *> MEtPhi_;
  vector<TH1F *> MEtSumEt_;
  vector<TH1F *> MEtmEtSig_;

  // Means
  MultiTH1F * Multi_Vertex_Dz_;
  MultiTH1F * Multi_Prim_Second_Vertex_Dz_;
  MultiTH1F * Multi_Vertex_Num_;

  MultiTH1F * Multi_PrimVNum_;
  MultiTH1F * Multi_PrimVPt_;
  MultiTH1F * Multi_PrimVEta_;
  MultiTH1F * Multi_PrimVPhi_;

  MultiTH1F * Multi_SecVNum_;
  MultiTH1F * Multi_SecVPt_;
  MultiTH1F * Multi_SecVEta_;
  MultiTH1F * Multi_SecVPhi_;

  MultiTH1F * Multi_AllSecVNum_;
  MultiTH1F * Multi_AllSecVPt_;
  MultiTH1F * Multi_AllSecVEta_;
  MultiTH1F * Multi_AllSecVPhi_;

  TH1F * DPhimin_;

  // Trigger efficiency counters
  // Multijet
  int Eff_;
  int Eff_et1_;
  int Eff_et2_;
  int Eff_et3_;
  int Eff_et4_;

  int Eff_cen_;
  int Eff_cen_et1_;
  int Eff_cen_et2_;
  int Eff_cen_et3_;
  int Eff_cen_et4_;

  int Eff_tau_;
  int Eff_tau_et1_;
  int Eff_tau_et2_;
  int Eff_tau_et3_;
  int Eff_tau_et4_;

  int Eff_for_;
  int Eff_for_et1_;
  int Eff_for_et2_;
  int Eff_for_et3_;
  int Eff_for_et4_;

  int Eff_nofor_;
  int Eff_nofor_et1_;
  int Eff_nofor_et2_;
  int Eff_nofor_et3_;
  int Eff_nofor_et4_;

  // MEt+Jet
  int Eff_MEtJet_;
  int Eff_MEtJet_cen_;
  int Eff_MEtJet_tau_;
  int Eff_MEtJet_for_;
  int Eff_MEtJet_nofor_;
  // Tau
  int Eff_tautrig_;
  int Eff_tautrig_single_;
  int Eff_tautrig_ditau_;

  // MultiJet Or MEtJet
  int Eff_Multi_Or_MEtJet_;
  int Eff_Multi_Or_MEtJet_nofor_;

  // Offline
  int offlineEffMultijet_;
  int offlineEffMEtJet_;
  int offlineEffTauTrig_;
  int offlineEff_Multi_Or_MEtJet_;
  int offlineEff_Multi_Or_MEtJet_nofor_;

  double dz_;
  double dzmax_;
  int bins_;

  // Directory in the root file to hold the multiple histograms
  TDirectory *DirVertexDz_;

  int *numgoodpjeff_;
  int *numgoodpjeff_3_;
  int *numgoodpjeff_4_;
  int *numgoodpjeff_5_;
  int *numgoodpjeff_6_;

  TH1F * EffNumGoodPj_;
  TH1F * EffNumGoodPj_3_;
  TH1F * EffNumGoodPj_4_;
  TH1F * EffNumGoodPj_5_;
  TH1F * EffNumGoodPj_6_;

  // ----------member data ---------------------------
};
