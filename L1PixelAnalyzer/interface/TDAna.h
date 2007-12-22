// -*- C++ -*-
//
// Package:    TDAna
// Class:      TDAna
// 
/**\class TDAna TDAna.cc AnalysisExamples/TDAna/src/TDAna.cc

 Description: <one line class summary>

 Implementation:
 This class shows how to access:
 - level 1 calorimetric quantities
 - offline corrected jets (calibration performed here)
 - offline corrected MET, depending on jets corrections
 - MC informations         <---------------------------------- to do
 - B tagging               <---------------------------------- to do

*/
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: TDAna.h,v 1.10 2007/12/22 14:01:36 dorigo Exp $
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

// Root include files
// ------------------
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"

// Data include files
// ------------------

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
// -----------------------
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"

// L1Trigger evaluator
// -------------------
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/HiVariables.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiStack.h"
#include "AnalysisExamples/AnalysisClasses/interface/L1PixelTrig.h"

// Class declaration
// -----------------
class TDAna : public edm::EDAnalyzer {
 public:
  explicit TDAna(const edm::ParameterSet&);
  ~TDAna();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one TDAna object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;
  
  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the TDAna constructor ).
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
  bool QCD_;
  std::string OutputEffFileName;

  TH1D * uncorr_JetPt_IC5_;
  TH1D * corr_JetPt_IC5_;
  TH1D * JetNumber_IC5_;

  TH1D * MEt_CorrIC5_Pt_;
  TH1D * MEt_CorrIC5_Phi_;
  TH1D * MEt_CorrIC5_SumEt_;
  TH1D * MEt_CorrIC5_mEtSig_;


  // Function histograms
  // -------------------
  TH1D * HSS_sig[8];
  TH1D * HSS_bgr[8];
  
  // My histograms
  // -------------
  TH1D * NJets_;
  TH1D * UncorrSumEt_;
  TH1D * UncorrHt_;
  TH1D * CorrSumEt_;
  TH1D * CorrHt_;
  TH1D * GoodSumEt_;
  TH1D * GoodHt_;
  TH1D * GoodHt2_;
  TH1D * MEt_;
  TH1D * MEtSig_;
  TH1D * MEtSigNew_;
  TH1D * MEtDPM_;
  TH1D * MEtDP1_;
  TH1D * MEtDP2_;
  TH1D * MEtDP3_;
  TH1D * UncorrMEtSig_;
  TH1D * CorrMEtSig_;
  TH1D * M3best_;
  TH1D * Mwbest_;
  TH1D * Chi2mass_;
  TH1D * M45best_;
  TH1D * Chi2ext_;
  TH2D * MEx_SumEt_;
  TH1D * DP12_;
  TH1D * DPbb_;
  TH1D * M_others_;
  TH1D * Mbbnoh_;
  TH1D * DPbbnoh_;
  TH1D * M8_;
  TH1D * C8_;
  TH1D * M45bestall_;
  TH1D * Chi2extall_;
  TH1D * DPbball_;
  TH1D * SumHED4_;
  TH1D * SumHPD4_;
  TH1D * SumHED6_;
  TH1D * SumHPD6_;

  TH1D * NJetsS_;
  TH1D * UncorrSumEtS_;
  TH1D * UncorrHtS_;
  TH1D * CorrSumEtS_;
  TH1D * CorrHtS_;
  TH1D * GoodSumEtS_;
  TH1D * GoodHtS_;
  TH1D * GoodHt2S_;
  TH1D * MEtS_;
  TH1D * MEtSigS_;
  TH1D * MEtSigNewS_;
  TH1D * MEtDPMS_;
  TH1D * MEtDP1S_; 
  TH1D * MEtDP2S_;
  TH1D * MEtDP3S_;
  TH1D * UncorrMEtSigS_;
  TH1D * CorrMEtSigS_;
  TH1D * M3bestS_;
  TH1D * MwbestS_;
  TH1D * Chi2massS_;
  TH1D * M45bestS_;
  TH1D * Chi2extS_;
  TH2D * MEx_SumEtS_;
  TH1D * DP12S_;
  TH1D * DPbbS_;
  TH1D * M_othersS_;
  TH1D * MbbnohS_;
  TH1D * DPbbnohS_;
  TH1D * M8S_;
  TH1D * C8S_;
  TH1D * M45bestallS_;
  TH1D * Chi2extallS_;
  TH1D * DPbballS_;
  TH1D * SumHED4S_;
  TH1D * SumHPD4S_;
  TH1D * SumHED6S_;
  TH1D * SumHPD6S_;

  TH1D * NJetsSS_;
  TH1D * UncorrSumEtSS_;
  TH1D * UncorrHtSS_;
  TH1D * CorrSumEtSS_;
  TH1D * CorrHtSS_;
  TH1D * GoodSumEtSS_;
  TH1D * GoodHtSS_;
  TH1D * GoodHt2SS_;
  TH1D * MEtSS_;
  TH1D * MEtSigSS_;
  TH1D * MEtSigNewSS_;
  TH1D * MEtDPMSS_;
  TH1D * MEtDP1SS_;
  TH1D * MEtDP2SS_;
  TH1D * MEtDP3SS_;
  TH1D * UncorrMEtSigSS_;
  TH1D * CorrMEtSigSS_;
  TH1D * M3bestSS_;
  TH1D * MwbestSS_;
  TH1D * Chi2massSS_;
  TH1D * M45bestSS_;
  TH1D * Chi2extSS_;
  TH2D * MEx_SumEtSS_;
  TH1D * DP12SS_;
  TH1D * DPbbSS_;
  TH1D * M_othersSS_;
  TH1D * MbbnohSS_;
  TH1D * DPbbnohSS_;
  TH1D * M8SS_;
  TH1D * C8SS_;
  TH1D * M45bestallSS_;
  TH1D * Chi2extallSS_;
  TH1D * DPbballSS_;
  TH1D * SumHED4SS_;
  TH1D * SumHPD4SS_;
  TH1D * SumHED6SS_;
  TH1D * SumHPD6SS_;

  TH1D * NJetsSSS_;
  TH1D * UncorrSumEtSSS_;
  TH1D * UncorrHtSSS_;
  TH1D * CorrSumEtSSS_;
  TH1D * CorrHtSSS_;
  TH1D * GoodSumEtSSS_;
  TH1D * GoodHtSSS_;
  TH1D * GoodHt2SSS_;
  TH1D * MEtSSS_;
  TH1D * MEtSigSSS_;
  TH1D * MEtSigNewSSS_;
  TH1D * MEtDPMSSS_;
  TH1D * MEtDP1SSS_;
  TH1D * MEtDP2SSS_;
  TH1D * MEtDP3SSS_;
  TH1D * UncorrMEtSigSSS_;
  TH1D * CorrMEtSigSSS_;
  TH1D * M3bestSSS_;
  TH1D * MwbestSSS_;
  TH1D * Chi2massSSS_;
  TH1D * M45bestSSS_;
  TH1D * Chi2extSSS_;
  TH2D * MEx_SumEtSSS_;
  TH1D * DP12SSS_;
  TH1D * DPbbSSS_;
  TH1D * M_othersSSS_;
  TH1D * MbbnohSSS_;
  TH1D * DPbbnohSSS_;
  TH1D * M8SSS_;
  TH1D * C8SSS_;
  TH1D * M45bestallSSS_;
  TH1D * Chi2extallSSS_;
  TH1D * DPbballSSS_;
  TH1D * SumHED4SSS_;
  TH1D * SumHPD4SSS_;
  TH1D * SumHED6SSS_;
  TH1D * SumHPD6SSS_;

  TH1D * N4NJSSS_;
  TH1D * E4NJSSS_;

  // Histograms with number of entries
  // ---------------------------------
  TH1D * NJetsN_;
  TH1D * UncorrSumEtN_;
  TH1D * UncorrHtN_;
  TH1D * CorrSumEtN_;
  TH1D * CorrHtN_;
  TH1D * GoodSumEtN_;
  TH1D * GoodHtN_;
  TH1D * GoodHt2N_;
  TH1D * MEtN_;
  TH1D * MEtSigN_;
  TH1D * MEtSigNewN_;
  TH1D * MEtDPMN_;
  TH1D * MEtDP1N_;
  TH1D * MEtDP2N_;
  TH1D * MEtDP3N_;
  TH1D * UncorrMEtSigN_;
  TH1D * CorrMEtSigN_;
  TH1D * M3bestN_;
  TH1D * MwbestN_;
  TH1D * Chi2massN_;
  TH1D * M45bestN_;
  TH1D * Chi2extN_;
  TH2D * MEx_SumEtN_;
  TH1D * DP12N_;
  TH1D * DPbbN_;
  TH1D * M_othersN_;
  TH1D * MbbnohN_;
  TH1D * DPbbnohN_;
  TH1D * M8N_;
  TH1D * C8N_;
  TH1D * M45bestallN_;
  TH1D * Chi2extallN_;
  TH1D * DPbballN_;
  TH1D * SumHED4N_;
  TH1D * SumHPD4N_;
  TH1D * SumHED6N_;
  TH1D * SumHPD6N_;

  TH1D * NJetsSN_;
  TH1D * UncorrSumEtSN_;
  TH1D * UncorrHtSN_;
  TH1D * CorrSumEtSN_;
  TH1D * CorrHtSN_;
  TH1D * GoodSumEtSN_;
  TH1D * GoodHtSN_;
  TH1D * GoodHt2SN_;
  TH1D * MEtSN_;
  TH1D * MEtSigSN_;
  TH1D * MEtSigNewSN_;
  TH1D * MEtDPMSN_;
  TH1D * MEtDP1SN_; 
  TH1D * MEtDP2SN_;
  TH1D * MEtDP3SN_;
  TH1D * UncorrMEtSigSN_;
  TH1D * CorrMEtSigSN_;
  TH1D * M3bestSN_;
  TH1D * MwbestSN_;
  TH1D * Chi2massSN_;
  TH1D * M45bestSN_;
  TH1D * Chi2extSN_;
  TH2D * MEx_SumEtSN_;
  TH1D * DP12SN_;
  TH1D * DPbbSN_;
  TH1D * M_othersSN_;
  TH1D * MbbnohSN_;
  TH1D * DPbbnohSN_;
  TH1D * M8SN_;
  TH1D * C8SN_;
  TH1D * M45bestallSN_;
  TH1D * Chi2extallSN_;
  TH1D * DPbballSN_;
  TH1D * SumHED4SN_;
  TH1D * SumHPD4SN_;
  TH1D * SumHED6SN_;
  TH1D * SumHPD6SN_;

  TH1D * NJetsSSN_;
  TH1D * UncorrSumEtSSN_;
  TH1D * UncorrHtSSN_;
  TH1D * CorrSumEtSSN_;
  TH1D * CorrHtSSN_;
  TH1D * GoodSumEtSSN_;
  TH1D * GoodHtSSN_;
  TH1D * GoodHt2SSN_;
  TH1D * MEtSSN_;
  TH1D * MEtSigSSN_;
  TH1D * MEtSigNewSSN_;
  TH1D * MEtDPMSSN_;
  TH1D * MEtDP1SSN_;
  TH1D * MEtDP2SSN_;
  TH1D * MEtDP3SSN_;
  TH1D * UncorrMEtSigSSN_;
  TH1D * CorrMEtSigSSN_;
  TH1D * M3bestSSN_;
  TH1D * MwbestSSN_;
  TH1D * Chi2massSSN_;
  TH1D * M45bestSSN_;
  TH1D * Chi2extSSN_;
  TH2D * MEx_SumEtSSN_;
  TH1D * DP12SSN_;
  TH1D * DPbbSSN_;
  TH1D * M_othersSSN_;
  TH1D * MbbnohSSN_;
  TH1D * DPbbnohSSN_;
  TH1D * M8SSN_;
  TH1D * C8SSN_;
  TH1D * M45bestallSSN_;
  TH1D * Chi2extallSSN_;
  TH1D * DPbballSSN_;
  TH1D * SumHED4SSN_;
  TH1D * SumHPD4SSN_;
  TH1D * SumHED6SSN_;
  TH1D * SumHPD6SSN_;

  TH1D * NJetsSSSN_;
  TH1D * UncorrSumEtSSSN_;
  TH1D * UncorrHtSSSN_;
  TH1D * CorrSumEtSSSN_;
  TH1D * CorrHtSSSN_;
  TH1D * GoodSumEtSSSN_;
  TH1D * GoodHtSSSN_;
  TH1D * GoodHt2SSSN_;
  TH1D * MEtSSSN_;
  TH1D * MEtSigSSSN_;
  TH1D * MEtSigNewSSSN_;
  TH1D * MEtDPMSSSN_;
  TH1D * MEtDP1SSSN_;
  TH1D * MEtDP2SSSN_;
  TH1D * MEtDP3SSSN_;
  TH1D * UncorrMEtSigSSSN_;
  TH1D * CorrMEtSigSSSN_;
  TH1D * M3bestSSSN_;
  TH1D * MwbestSSSN_;
  TH1D * Chi2massSSSN_;
  TH1D * M45bestSSSN_;
  TH1D * Chi2extSSSN_;
  TH2D * MEx_SumEtSSSN_;
  TH1D * DP12SSSN_;
  TH1D * DPbbSSSN_;
  TH1D * M_othersSSSN_;
  TH1D * MbbnohSSSN_;
  TH1D * DPbbnohSSSN_;
  TH1D * M8SSSN_;
  TH1D * C8SSSN_;
  TH1D * M45bestallSSSN_;
  TH1D * Chi2extallSSSN_;
  TH1D * DPbballSSSN_;
  TH1D * SumHED4SSSN_;
  TH1D * SumHPD4SSSN_;
  TH1D * SumHED6SSSN_;
  TH1D * SumHPD6SSSN_;

  // Error histograms
  // ----------------
  TH1D * NJetsW_;
  TH1D * UncorrSumEtW_;
  TH1D * UncorrHtW_;
  TH1D * CorrSumEtW_;
  TH1D * CorrHtW_;
  TH1D * GoodSumEtW_;
  TH1D * GoodHtW_;
  TH1D * GoodHt2W_;
  TH1D * MEtW_;
  TH1D * MEtSigW_;
  TH1D * MEtSigNewW_;
  TH1D * MEtDPMW_;
  TH1D * MEtDP1W_;
  TH1D * MEtDP2W_;
  TH1D * MEtDP3W_;
  TH1D * UncorrMEtSigW_;
  TH1D * CorrMEtSigW_;
  TH1D * M3bestW_;
  TH1D * MwbestW_;
  TH1D * Chi2massW_;
  TH1D * M45bestW_;
  TH1D * Chi2extW_;
  TH2D * MEx_SumEtW_;
  TH1D * DP12W_;
  TH1D * DPbbW_;
  TH1D * M_othersW_;
  TH1D * MbbnohW_;
  TH1D * DPbbnohW_;
  TH1D * M8W_;
  TH1D * C8W_;
  TH1D * M45bestallW_;
  TH1D * Chi2extallW_;
  TH1D * DPbballW_;
  TH1D * SumHED4W_;
  TH1D * SumHPD4W_;
  TH1D * SumHED6W_;
  TH1D * SumHPD6W_;

  TH1D * NJetsSW_;
  TH1D * UncorrSumEtSW_;
  TH1D * UncorrHtSW_;
  TH1D * CorrSumEtSW_;
  TH1D * CorrHtSW_;
  TH1D * GoodSumEtSW_;
  TH1D * GoodHtSW_;
  TH1D * GoodHt2SW_;
  TH1D * MEtSW_;
  TH1D * MEtSigSW_;
  TH1D * MEtSigNewSW_;
  TH1D * MEtDPMSW_;
  TH1D * MEtDP1SW_; 
  TH1D * MEtDP2SW_;
  TH1D * MEtDP3SW_;
  TH1D * UncorrMEtSigSW_;
  TH1D * CorrMEtSigSW_;
  TH1D * M3bestSW_;
  TH1D * MwbestSW_;
  TH1D * Chi2massSW_;
  TH1D * M45bestSW_;
  TH1D * Chi2extSW_;
  TH2D * MEx_SumEtSW_;
  TH1D * DP12SW_;
  TH1D * DPbbSW_;
  TH1D * M_othersSW_;
  TH1D * MbbnohSW_;
  TH1D * DPbbnohSW_;
  TH1D * M8SW_;
  TH1D * C8SW_;
  TH1D * M45bestallSW_;
  TH1D * Chi2extallSW_;
  TH1D * DPbballSW_;
  TH1D * SumHED4SW_;
  TH1D * SumHPD4SW_;
  TH1D * SumHED6SW_;
  TH1D * SumHPD6SW_;

  TH1D * NJetsSSW_;
  TH1D * UncorrSumEtSSW_;
  TH1D * UncorrHtSSW_;
  TH1D * CorrSumEtSSW_;
  TH1D * CorrHtSSW_;
  TH1D * GoodSumEtSSW_;
  TH1D * GoodHtSSW_;
  TH1D * GoodHt2SSW_;
  TH1D * MEtSSW_;
  TH1D * MEtSigSSW_;
  TH1D * MEtSigNewSSW_;
  TH1D * MEtDPMSSW_;
  TH1D * MEtDP1SSW_;
  TH1D * MEtDP2SSW_;
  TH1D * MEtDP3SSW_;
  TH1D * UncorrMEtSigSSW_;
  TH1D * CorrMEtSigSSW_;
  TH1D * M3bestSSW_;
  TH1D * MwbestSSW_;
  TH1D * Chi2massSSW_;
  TH1D * M45bestSSW_;
  TH1D * Chi2extSSW_;
  TH2D * MEx_SumEtSSW_;
  TH1D * DP12SSW_;
  TH1D * DPbbSSW_;
  TH1D * M_othersSSW_;
  TH1D * MbbnohSSW_;
  TH1D * DPbbnohSSW_;
  TH1D * M8SSW_;
  TH1D * C8SSW_;
  TH1D * M45bestallSSW_;
  TH1D * Chi2extallSSW_;
  TH1D * DPbballSSW_;
  TH1D * SumHED4SSW_;
  TH1D * SumHPD4SSW_;
  TH1D * SumHED6SSW_;
  TH1D * SumHPD6SSW_;

  TH1D * NJetsSSSW_;
  TH1D * UncorrSumEtSSSW_;
  TH1D * UncorrHtSSSW_;
  TH1D * CorrSumEtSSSW_;
  TH1D * CorrHtSSSW_;
  TH1D * GoodSumEtSSSW_;
  TH1D * GoodHtSSSW_;
  TH1D * GoodHt2SSSW_;
  TH1D * MEtSSSW_;
  TH1D * MEtSigSSSW_;
  TH1D * MEtSigNewSSSW_;
  TH1D * MEtDPMSSSW_;
  TH1D * MEtDP1SSSW_;
  TH1D * MEtDP2SSSW_;
  TH1D * MEtDP3SSSW_;
  TH1D * UncorrMEtSigSSSW_;
  TH1D * CorrMEtSigSSSW_;
  TH1D * M3bestSSSW_;
  TH1D * MwbestSSSW_;
  TH1D * Chi2massSSSW_;
  TH1D * M45bestSSSW_;
  TH1D * Chi2extSSSW_;
  TH2D * MEx_SumEtSSSW_;
  TH1D * DP12SSSW_;
  TH1D * DPbbSSSW_;
  TH1D * M_othersSSSW_;
  TH1D * MbbnohSSSW_;
  TH1D * DPbbnohSSSW_;
  TH1D * M8SSSW_;
  TH1D * C8SSSW_;
  TH1D * M45bestallSSSW_;
  TH1D * Chi2extallSSSW_;
  TH1D * DPbballSSSW_;
  TH1D * SumHED4SSSW_;
  TH1D * SumHPD4SSSW_;
  TH1D * SumHED6SSSW_;
  TH1D * SumHPD6SSSW_;

  TH1D * N4NJSSSW_;
  TH1D * E4NJSSSW_;

  TH2D * UncorrMEt_SumEt_;
  TH2D * CorrMEt_SumEt_; 
  TH2D * MEt_SumEt_; 
  TH2D * UncorrMEt_SumEtC_; 
  TH2D * CorrMEt_SumEtC_; 
  TH2D * MEt_SumEtC_; 
  TH2D * UncorrMEt_SumEtJ_; 
  TH2D * CorrMEt_SumEtJ_; 
  TH2D * MEt_SumEtJ_; 

  float loose_;
  float medium_;
  float tight_;

  double rel_lik;
  double njsss;

  // Likelihood histograms
  // ---------------------
  TH1D * L_;
  TH1D * LS_;
  TH1D * LSS_;
  TH1D * LSSS_;
  TH1D * LW_;
  TH1D * LSW_;
  TH1D * LSSW_;
  TH1D * LSSSW_;
  TH1D * LN_;
  TH1D * LSN_;
  TH1D * LSSN_;
  TH1D * LSSSN_;

  // PTag numbers
  // ------------
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

  double PHETL[1000];
  double PHETM[1000];
  double PHETT[1000];
  double PHPTL[1000];
  double PHPTM[1000];
  double PHPTT[1000];
  double PHETLS[1000];
  double PHETMS[1000];
  double PHETTS[1000];
  double PHPTLS[1000];
  double PHPTMS[1000];
  double PHPTTS[1000];

};
