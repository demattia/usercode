// -*- C++ -*-
//
// Package:    BaseTDAna
// Class:      BaseTDAna
// 
/**\class BaseTDAna BaseTDAna.cc AnalysisExamples/BaseTDAna/src/BaseTDAna.cc

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
// $Id: BaseTDAna.h,v 1.2 2007/12/06 12:08:22 tosi Exp $
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
#include "AnalysisExamples/AnalysisClasses/interface/L1PixelTrig.h"

//
// class declaration
//

class BaseTDAna : public edm::EDAnalyzer {
 public:
  explicit BaseTDAna(const edm::ParameterSet&);
  ~BaseTDAna();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one BaseTDAna object is created
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the BaseTDAna constructor ).
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
  edm::InputTag simplePixelJetLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag summaryLabel_;

  unsigned int numTkCut_;
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

  // Offline
  int offlineEffMultijet_;
  int offlineEffMEtJet_;
  int offlineEffTauTrig_;

  double dz_;
  double dzmax_;
  int bins_;
  double dR_;
  double dRmax_;

  // Pixel trigger efficiency
  TH1F ** EffMultijetPixel_;
  TH1F ** EffMultijetPixelEt1_;
  TH1F ** EffMultijetPixelEt2_;
  TH1F ** EffMultijetPixelEt3_;
  TH1F ** EffMultijetPixelEt4_;
  TH1F ** EffMEtJetPixel_;
  int EffMultijetPixelSize_;
  int EffMultijetPixelSizeEt1_;
  int EffMultijetPixelSizeEt2_;
  int EffMultijetPixelSizeEt3_;
  int EffMultijetPixelSizeEt4_;
  int EffMEtJetPixelSize_;
  int ** EffMultijetPixelArray_;
  int ** EffMultijetPixelArrayEt1_;
  int ** EffMultijetPixelArrayEt2_;
  int ** EffMultijetPixelArrayEt3_;
  int ** EffMultijetPixelArrayEt4_;
  int ** EffMEtJetPixelArray_;

  // Offline efficiency
  TH1F ** offlineEffMultijetPixel_;
  TH1F ** offlineEffMEtJetPixel_;
  int offlineEffMultijetPixelSize_;
  int offlineEffMEtJetPixelSize_;
  int ** offlineEffMultijetPixelArray_;
  int ** offlineEffMEtJetPixelArray_;

  // Directory in the root file to hold the multiple histograms
  TDirectory *DirVertexDz_;

  // PixelTrigger alone efficiency
  int pixelTrig_3_;
  int pixelTrig_4_;
  int pixelTrig_5_;
  int pixelTrig_6_;

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


  // b tag discriminator histograms
  MultiTH1F * Multi_bDiscrHighEff;
  MultiTH1F * Multi_bDiscrHighPur;
  MultiTH1F * Multi_nobDiscrHighEff;
  MultiTH1F * Multi_nobDiscrHighPur;
  MultiTH1F * Multi_tagMassS1;
  MultiTH1F * Multi_tagMassS2;
  MultiTH1F * Multi_tagMassS3;
  MultiTProfile * Multi_jetEtVSbPt;

  // b tag discriminator histograms
  TH1F * bHighEff;   
  TH1F * bHighPur;   
  TH1F * nobHighEff; 
  TH1F * nobHighPur; 
  TH1F * deltaR;
  TH2F * jetEtVSbParton;
  TH2F * jetUncorrEtVSbParton;
TH1F * tagTkMassS1;
TH1F * tagTkMassS2;
TH1F * tagTkMassS3;

TH2F * mPdiVsDiscHighEff;
TH2F * uncorrEtVsDiscHighEff;
TH2F * emEnergyFractionVsDiscHighEff;
TH2F * jetMassVsDiscHighEff;
TH2F * tkNumS1VsDiscHighEff;
TH2F * tkSumPtS1VsDiscHighEff;
TH2F * tagTkMassS1VsDiscHighEff;
TH2F * tkNumS2VsDiscHighEff;
TH2F * tkSumPtS2VsDiscHighEff;
TH2F * tagTkMassS2VsDiscHighEff;
TH2F * tkNumS3VsDiscHighEff;
TH2F * tkSumPtS3VsDiscHighEff;
TH2F * tagTkMassS3VsDiscHighEff;

  // ----------member data ---------------------------
};
