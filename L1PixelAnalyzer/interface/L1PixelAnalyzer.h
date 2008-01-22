// -*- C++ -*-
//
// Package:    L1PixelAnalyzer
// Class:      L1PixelAnalyzer
// 
/**\class L1PixelAnalyzer L1PixelAnalyzer.cc AnalysisExamples/L1PixelAnalyzer/src/L1PixelAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: L1PixelAnalyzer.h,v 1.6 2007/11/09 17:25:29 demattia Exp $
//
//

// system include files
#include <memory>

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

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"

// L1Trigger evaluator
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"

//
// class declaration
//

class L1PixelAnalyzer : public edm::EDAnalyzer {
 public:
  explicit L1PixelAnalyzer(const edm::ParameterSet&);
  ~L1PixelAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int eventcounter;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  TH1F* L1ExtraCenJetsEt_;
  TH1F* L1ExtraCenJetsEta_;
  TH1F* L1ExtraCenJetsPhi_;

  TH1F* L1ExtraForJetsEt_;
  TH1F* L1ExtraForJetsEta_;
  TH1F* L1ExtraForJetsPhi_;

  TH1F* L1ExtraTauJetsEt_;
  TH1F* L1ExtraTauJetsEta_;
  TH1F* L1ExtraTauJetsPhi_;

  TH1F* L1ExtraIsoEmEt_;
  TH1F* L1ExtraIsoEmEta_;
  TH1F* L1ExtraIsoEmPhi_;

  TH1F* L1ExtraNonIsoEmEt_;
  TH1F* L1ExtraNonIsoEmEta_;
  TH1F* L1ExtraNonIsoEmPhi_;

  TH1F* L1ExtraEtMiss_;
  TH1F* L1ExtraEtMissPhi_;
  TH1F* L1ExtraEtTotal_;
  TH1F* L1ExtraEtHad_;

  TH1F* PixelTrack_P_X_;
  TH1F* PixelTrack_P_Y_;
  TH1F* PixelTrack_P_Z_;
  TH1F* PixelTrack_Pt_;
  TH1F* PixelTrack_Eta_;
  TH1F* PixelTrack_Phi_;
  TH1F* PixelTrack_Charge_;
  TH1F* PixelTrack_Chi2_;
  TH1F* PixelTrack_Zip_;
  TH1F* PixelTrack_Tip_;
  TH1F* PixelTrack_Vertex_X_; 
  TH1F* PixelTrack_Vertex_Y_; 
  TH1F* PixelTrack_Vertex_Z_; 

  TH1F* PixelHit_X_;
  TH1F* PixelHit_Y_;
  TH1F* PixelHit_Z_;
  TH2F* PixelHit_XY_;

  TH1F* PixelJet_Pt_;
  TH1F* PixelJet_Eta_;
  TH1F* PixelJet_Phi_;
  TH1F* PixelJet_NumTk_;
  TH1F* PixelJet_Vertex_Z_; 

  TH1F* GenJet_Pt_;
  TH1F* GenJet_Eta_;
  TH1F* GenJet_Phi_;

  TProfile* PJ_PtRes_;
  TH1F* PJ_EtaRes_;
  TH1F* PJ_PhiRes_;
  TH1F* PJ_DeltaRres_;
  TProfile* PJ_L1J_PtRes_;
  TProfile* L1J_PtRes_;
  TH1F* L1J_EtaRes_;
  TH1F* L1J_PhiRes_;
  TH1F* L1J_DeltaRres_;

  // ----------member data ---------------------------
};
