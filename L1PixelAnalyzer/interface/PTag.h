// -*- C++ -*-
//
// Package:    PTag
// Class:      PTag
// 
/**\class PTag PTag.cc AnalysisExamples/PTag/src/PTag.cc

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
// $Id: PTag.h,v 1.3 2007/12/10 12:27:41 tosi Exp $
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

// // // GCT and RCT data formats
// // #include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
// // #include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"
// // #include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

// // L1Extra
// // #include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
// #include "DataFormats/Candidate/interface/Candidate.h"
// #include "DataFormats/Candidate/interface/CandidateFwd.h"

// #include "DataFormats/L1Trigger/interface/L1EmParticle.h"
// #include "DataFormats/L1Trigger/interface/L1JetParticle.h"
// #include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
// #include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
// #include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

// #include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
// #include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// #include "FastSimulation/L1CaloTriggerProducer/interface/FastL1Region.h"
// // No BitInfos for release versions
// #include "FastSimDataFormats/External/interface/FastL1BitInfo.h"

// #include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
// #include "DataFormats/Math/interface/Vector3D.h"

// // L1 Pixel
// // #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
// // #include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
// // #include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
// // #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// // #include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
// // #include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
// #include "FWCore/Framework/interface/ESHandle.h"
// #include "DataFormats/GeometryVector/interface/GlobalPoint.h"
// #include "DataFormats/GeometryVector/interface/GlobalVector.h"
// #include "DataFormats/GeometryVector/interface/LocalVector.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "TrackingTools/Records/interface/TransientRecHitRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
// #include "DataFormats/DetId/interface/DetId.h" 
// #include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
// #include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
// #include "DataFormats/GeometryVector/interface/LocalPoint.h"
// #include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// #include "DataFormats/TrackReco/interface/TrackFwd.h"
// #include "DataFormats/TrackReco/interface/Track.h"

// #include "AnalysisExamples/PixelJet/interface/PixelJet.h"

// // GenJets
// #include "DataFormats/JetReco/interface/GenJetCollection.h"
// #include "DataFormats/JetReco/interface/GenJet.h"

// // Calo MEt
// #include "DataFormats/METReco/interface/CaloMET.h"


// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"

// L1Trigger evaluator
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiStack.h"
#include "AnalysisExamples/AnalysisClasses/interface/L1PixelTrig.h"

//
// class declaration
//

class PTag : public edm::EDAnalyzer {
 public:
  explicit PTag(const edm::ParameterSet&);
  ~PTag();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one PTag object is created
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the PTag constructor ).
  //  HiVariables HiVar;

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

  unsigned int numTkCut_;
  std::string OutputEffFileName;

  // Counters for Ptag
  // -----------------
  double N0HETL[1000];
  double N1HETL[1000];
  double N0HETM[1000];
  double N1HETM[1000];
  double N0HETT[1000];
  double N1HETT[1000];
  double N0HPTL[1000];
  double N1HPTL[1000];
  double N0HPTM[1000];
  double N1HPTM[1000];
  double N0HPTT[1000];
  double N1HPTT[1000];

  double loose_;
  double medium_;
  double tight_;

  // ----------member data ---------------------------
};
