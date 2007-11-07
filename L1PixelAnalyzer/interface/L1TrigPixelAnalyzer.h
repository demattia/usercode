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
// #include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "FastSimulation/L1CaloTriggerProducer/interface/FastL1Region.h"
// No BitInfos for release versions
#include "FastSimDataFormats/External/interface/FastL1BitInfo.h"

#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "DataFormats/Math/interface/Vector3D.h"

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
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Calo MEt
#include "DataFormats/METReco/interface/CaloMET.h"


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
  int Eff_cen_;
  int Eff_tau_;
  int Eff_for_;
  int Eff_nofor_;
  // MEt+Jet
  int Eff_MEtJet_;
  int Eff_MEtJet_cen_;
  int Eff_MEtJet_tau_;
  int Eff_MEtJet_for_;
  int Eff_MEtJet_nofor_;
  // Tau
  int Eff_tautrig_;

  double dz_;
  double dzmax_;
  int bins_;

  // Pixel trigger efficiency
  TH1F ** EffMultijetPixel_;
  TH1F ** EffMEtJetPixel_;
  int EffMultijetPixelSize_;
  int EffMEtJetPixelSize_;
  int ** EffMultijetPixelArray_;
  int ** EffMEtJetPixelArray_;

  // Directory in the root file to hold the multiple histograms
  TDirectory *DirVertexDz_;

  // ----------member data ---------------------------
};
