// -*- C++ -*-
//
// Package:    L1TrigPixelAnalyzer
// Class:      L1TrigPixelAnalyzer
// 
/**\class L1TrigPixelAnalyzer L1TrigPixelAnalyzer.cc AnalysisExamples/L1TrigPixelAnalyzer/src/L1TrigPixelAnalyzer.cc

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
// $Id: L1TrigPixelAnalyzer.h,v 1.2 2007/07/23 14:57:30 demattia Exp $
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

// Data includes
// -------------

// // GCT and RCT data formats
// #include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
// #include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"
// #include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

// L1Extra
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"

// L1 Pixel
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
// #include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
// #include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
// #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// #include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
// #include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "AnalysisExamples/PixelJet/interface/PixelJet.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetfwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Associator for the jets
#include "../../PJVERTEX_CMSSW/Classes/Associator/Associator.h"

// L1Trigger evaluator
#include "../../PJVERTEX_CMSSW/Classes/L1Trig/L1Trig.C"
#include "../../PJVERTEX_CMSSW/Classes/HiVariables/HiVariables.cc"
#include "../../PJVERTEX_CMSSW/Classes/MultiTH1F/MultiTH1F.h"

//
// class declaration
//

class L1TrigPixelAnalyzer : public edm::EDAnalyzer {
 public:
  explicit L1TrigPixelAnalyzer(const edm::ParameterSet&);
  ~L1TrigPixelAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one L1TrigPixelAnalyzer object is created
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the L1TrigPixelAnalyzer constructor ).
  //  HiVariables HiVar;

  // Use a dynamic construction, or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.
  HiVariables * HiVar;

  std::string CaloJetAlgorithm, JetCorrectionService, METCollection;
  std::string genParticleCandidates;
  unsigned int numTkCut;
  std::string OutputEffFileName;

  TH1F * uncorr_JetPt_IC5_;
  TH1F * corr_JetPt_IC5_;
  TH1F * JetNumber_IC5_;

  TH1F * MEt_CorrIC5_Pt_;
  TH1F * MEt_CorrIC5_Phi_;
  TH1F * MEt_CorrIC5_SumEt_;
  TH1F * MEt_CorrIC5_mEtSig_;

  TH1F * PixelJet_dz_;
  TH1F * PixelJet_Num_;
  TH1F * PixelJet_Track_Num_;

  std::vector<TH1F *> Vertex_Dz_;
  std::vector<TH1F *> Vertex_Num_;
  std::vector<TH1F *> Prim_Second_Vertex_Dz_;

  std::vector<TH1F *> PrimVNum_;
  std::vector<TH1F *> PrimVPt_;
  std::vector<TH1F *> PrimVEta_;
  std::vector<TH1F *> PrimVPhi_;

  std::vector<TH1F *> SecVNum_;
  std::vector<TH1F *> SecVPt_;
  std::vector<TH1F *> SecVEta_;
  std::vector<TH1F *> SecVPhi_;

  std::vector<TH1F *> AllSecVNum_;
  std::vector<TH1F *> AllSecVPt_;
  std::vector<TH1F *> AllSecVEta_;
  std::vector<TH1F *> AllSecVPhi_;

  // Means
  TH1F * Vertex_Dz_Mean_;
  TH1F * Vertex_Num_Mean_;
  TH1F * Prim_Second_Vertex_Dz_Mean_;

  MultiTH1F * Multi_Vertex_Dz_;

  TH1F * PrimVNum_Mean_;
  TH1F * PrimVPt_Mean_;
  TH1F * PrimVEta_Mean_;
  TH1F * PrimVPhi_Mean_;

  TH1F * SecVNum_Mean_;
  TH1F * SecVPt_Mean_;
  TH1F * SecVEta_Mean_;
  TH1F * SecVPhi_Mean_;

  TH1F * AllSecVNum_Mean_;
  TH1F * AllSecVPt_Mean_;
  TH1F * AllSecVEta_Mean_;
  TH1F * AllSecVPhi_Mean_;

  TH1F * DPhimin_;

  int Eff_;
  double dz_;
  double dzmax_;
  int bins_;

  // Directory in the root file to hold the multiple histograms
  TDirectory *DirVertexDz_;

  // ----------member data ---------------------------
};
