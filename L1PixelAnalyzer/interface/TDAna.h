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
// $Id: TDAna.h,v 1.26 2008/01/21 08:12:47 dorigo Exp $
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
  TH1D * HSS_sig[11];
  TH1D * HSS_bgr[11];
  TH1D * HEDpdf[8];
  TH1D * HPDpdf[8];
  TH1D * MTS1pdf[8];
  TH1D * MNS1pdf[8];
  TH1D * MTS2pdf[8];
  TH1D * MNS2pdf[8];
  TH1D * MTS3pdf[8];
  TH1D * MNS3pdf[8];


  // My histograms
  // -------------
  TH2D * Drmax_;
  TH2D * Drmedall_;
  TH2D * Drmed07_;
  TH2D * N07_;
  TH2D * N04_;
  TH2D * N02_;
  TH2D * Nlo_;
  TH2D * Detmedall_;
  TH2D * Detmed07_;
  TH2D * Perf07_;
  TH2D * Perf04_;
  TH2D * Perf02_;
  TH2D * Det2medall_;
  TH2D * Det2med07_;
  TH2D * Hrecfrac_;
  TH2D * Trecfrac_;

  TH1D * MHbest_;
  TH1D * MTbest_;
  TH1D * MWbest_;
  TH1D * HBJ_etrank_;
  TH1D * Hpt_;
  TH1D * Heta_;
  TH1D * Hdr_;
  TH1D * MHnot_;
  TH1D * Hnotpt_;
  TH1D * Hnoteta_;
  TH1D * Hnotdr_;
  TH1D * MTnotbest_;
  TH1D * Tpt_;
  TH1D * Teta_;
  TH1D * THdeta_;
  TH1D * THdphi_;
  TH1D * THproj_;
  TH1D * Tnotpt_;
  TH1D * Tnoteta_;
  TH1D * THnotdeta_;
  TH1D * THnotdphi_;
  TH1D * THnotproj_;

  TH1D * HED1_;
  TH1D * HPD1_;
  TH1D * HED2_;
  TH1D * HPD2_;
  TH1D * HED3_;
  TH1D * HPD3_;
  TH1D * HED4_;
  TH1D * HPD4_;
  TH1D * HED5_;
  TH1D * HPD5_;
  TH1D * HED6_;
  TH1D * HPD6_;
  TH1D * HED7_;
  TH1D * HPD7_;
  TH1D * HED8_;
  TH1D * HPD8_;

  TProfile * DEtb_prof_;
  TProfile * DEtq_prof_;
  TProfile * DEtcb_prof_;
  TProfile * DEtcq_prof_;

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
  TH1D * M6_;
  TH1D * C6_;
  TH1D * M8_;
  TH1D * C8_;
  TH1D * M45bestall_;
  TH1D * Chi2extall_;
  TH1D * DPbball_;
  TH1D * SumHED4_;
  TH1D * SumHPD4_;
  TH1D * SumHED6_;
  TH1D * SumHPD6_;
  TH1D * HED_;
  TH1D * HPD_;
  TH1D * Et6_;
  TH1D * Mwmin_;
  TH1D * Hbestcomb_;
  TH1D * Drpairbestall_;
  TH1D * M3a_;
  TH1D * Mwa_;
  TH1D * Scprod_;
  TH1D * Thdeta_;
  TH1D * M5_;
  TH1D * TTMS1_;
  TH1D * TTMS2_;
  TH1D * TTMS3_;

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
  TH1D * M6S_;
  TH1D * C6S_;
  TH1D * M8S_;
  TH1D * C8S_;
  TH1D * M45bestallS_;
  TH1D * Chi2extallS_;
  TH1D * DPbballS_;
  TH1D * SumHED4S_;
  TH1D * SumHPD4S_;
  TH1D * SumHED6S_;
  TH1D * SumHPD6S_;
  TH1D * HEDS_;
  TH1D * HPDS_;
  TH1D * Et6S_;
  TH1D * MwminS_;
  TH1D * HbestcombS_;
  TH1D * DrpairbestallS_;
  TH1D * M3aS_;
  TH1D * MwaS_;
  TH1D * ScprodS_;
  TH1D * ThdetaS_;
  TH1D * M5S_;
  TH1D * TTMS1S_;
  TH1D * TTMS2S_;
  TH1D * TTMS3S_;

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
  TH1D * M6SS_;
  TH1D * C6SS_;
  TH1D * M8SS_;
  TH1D * C8SS_;
  TH1D * M45bestallSS_;
  TH1D * Chi2extallSS_;
  TH1D * DPbballSS_;
  TH1D * SumHED4SS_;
  TH1D * SumHPD4SS_;
  TH1D * SumHED6SS_;
  TH1D * SumHPD6SS_;
  TH1D * HEDSS_;
  TH1D * HPDSS_;
  TH1D * Et6SS_;
  TH1D * MwminSS_;
  TH1D * HbestcombSS_;
  TH1D * DrpairbestallSS_;
  TH1D * M3aSS_;
  TH1D * MwaSS_;
  TH1D * ScprodSS_;
  TH1D * ThdetaSS_;
  TH1D * M5SS_;
  TH1D * TTMS1SS_;
  TH1D * TTMS2SS_;
  TH1D * TTMS3SS_;

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
  TH1D * M6SSS_;
  TH1D * C6SSS_;
  TH1D * M8SSS_;
  TH1D * C8SSS_;
  TH1D * M45bestallSSS_;
  TH1D * Chi2extallSSS_;
  TH1D * DPbballSSS_;
  TH1D * SumHED4SSS_;
  TH1D * SumHPD4SSS_;
  TH1D * SumHED6SSS_;
  TH1D * SumHPD6SSS_;
  TH1D * HEDSSS_;
  TH1D * HPDSSS_;
  TH1D * Et6SSS_;
  TH1D * MwminSSS_;
  TH1D * HbestcombSSS_;
  TH1D * DrpairbestallSSS_;
  TH1D * M3aSSS_;
  TH1D * MwaSSS_;
  TH1D * ScprodSSS_;
  TH1D * ThdetaSSS_;
  TH1D * M5SSS_;
  TH1D * TTMS1SSS_;
  TH1D * TTMS2SSS_;
  TH1D * TTMS3SSS_;

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
  TH1D * M6N_;
  TH1D * C6N_;
  TH1D * M8N_;
  TH1D * C8N_;
  TH1D * M45bestallN_;
  TH1D * Chi2extallN_;
  TH1D * DPbballN_;
  TH1D * SumHED4N_;
  TH1D * SumHPD4N_;
  TH1D * SumHED6N_;
  TH1D * SumHPD6N_;
  TH1D * HEDN_;
  TH1D * HPDN_;
  TH1D * Et6N_;
  TH1D * MwminN_;
  TH1D * HbestcombN_;
  TH1D * DrpairbestallN_;
  TH1D * M3aN_;
  TH1D * MwaN_;
  TH1D * ScprodN_;
  TH1D * ThdetaN_;
  TH1D * M5N_;
  TH1D * TTMS1N_;
  TH1D * TTMS2N_;
  TH1D * TTMS3N_;

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
  TH1D * M6SN_;
  TH1D * C6SN_;
  TH1D * M8SN_;
  TH1D * C8SN_;
  TH1D * M45bestallSN_;
  TH1D * Chi2extallSN_;
  TH1D * DPbballSN_;
  TH1D * SumHED4SN_;
  TH1D * SumHPD4SN_;
  TH1D * SumHED6SN_;
  TH1D * SumHPD6SN_;
  TH1D * HEDSN_;
  TH1D * HPDSN_;
  TH1D * Et6SN_;
  TH1D * MwminSN_;
  TH1D * HbestcombSN_;
  TH1D * DrpairbestallSN_;
  TH1D * M3aSN_;
  TH1D * MwaSN_;
  TH1D * ScprodSN_;
  TH1D * ThdetaSN_;
  TH1D * M5SN_;
  TH1D * TTMS1SN_;
  TH1D * TTMS2SN_;
  TH1D * TTMS3SN_;

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
  TH1D * M6SSN_;
  TH1D * C6SSN_;
  TH1D * M8SSN_;
  TH1D * C8SSN_;
  TH1D * M45bestallSSN_;
  TH1D * Chi2extallSSN_;
  TH1D * DPbballSSN_;
  TH1D * SumHED4SSN_;
  TH1D * SumHPD4SSN_;
  TH1D * SumHED6SSN_;
  TH1D * SumHPD6SSN_;
  TH1D * HEDSSN_;
  TH1D * HPDSSN_;
  TH1D * Et6SSN_;
  TH1D * MwminSSN_;
  TH1D * HbestcombSSN_;
  TH1D * DrpairbestallSSN_;
  TH1D * M3aSSN_;
  TH1D * MwaSSN_;
  TH1D * ScprodSSN_;
  TH1D * ThdetaSSN_;
  TH1D * M5SSN_;
  TH1D * TTMS1SSN_;
  TH1D * TTMS2SSN_;
  TH1D * TTMS3SSN_;

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
  TH1D * M6SSSN_;
  TH1D * C6SSSN_;
  TH1D * M8SSSN_;
  TH1D * C8SSSN_;
  TH1D * M45bestallSSSN_;
  TH1D * Chi2extallSSSN_;
  TH1D * DPbballSSSN_;
  TH1D * SumHED4SSSN_;
  TH1D * SumHPD4SSSN_;
  TH1D * SumHED6SSSN_;
  TH1D * SumHPD6SSSN_;
  TH1D * HEDSSSN_;
  TH1D * HPDSSSN_;
  TH1D * Et6SSSN_;
  TH1D * MwminSSSN_;
  TH1D * HbestcombSSSN_;
  TH1D * DrpairbestallSSSN_;
  TH1D * M3aSSSN_;
  TH1D * MwaSSSN_;
  TH1D * ScprodSSSN_;
  TH1D * ThdetaSSSN_;
  TH1D * M5SSSN_;
  TH1D * TTMS1SSSN_;
  TH1D * TTMS2SSSN_;
  TH1D * TTMS3SSSN_;

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
  TH1D * M6W_;
  TH1D * C6W_;
  TH1D * M8W_;
  TH1D * C8W_;
  TH1D * M45bestallW_;
  TH1D * Chi2extallW_;
  TH1D * DPbballW_;
  TH1D * SumHED4W_;
  TH1D * SumHPD4W_;
  TH1D * SumHED6W_;
  TH1D * SumHPD6W_;
  TH1D * HEDW_;
  TH1D * HPDW_;
  TH1D * Et6W_;
  TH1D * MwminW_;
  TH1D * HbestcombW_;
  TH1D * DrpairbestallW_;
  TH1D * M3aW_;
  TH1D * MwaW_;
  TH1D * ScprodW_;
  TH1D * ThdetaW_;
  TH1D * M5W_;
  TH1D * TTMS1W_;
  TH1D * TTMS2W_;
  TH1D * TTMS3W_;

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
  TH1D * M6SW_;
  TH1D * C6SW_;
  TH1D * M8SW_;
  TH1D * C8SW_;
  TH1D * M45bestallSW_;
  TH1D * Chi2extallSW_;
  TH1D * DPbballSW_;
  TH1D * SumHED4SW_;
  TH1D * SumHPD4SW_;
  TH1D * SumHED6SW_;
  TH1D * SumHPD6SW_;
  TH1D * HEDSW_;
  TH1D * HPDSW_;
  TH1D * Et6SW_;
  TH1D * MwminSW_;
  TH1D * HbestcombSW_;
  TH1D * DrpairbestallSW_;
  TH1D * M3aSW_;
  TH1D * MwaSW_;
  TH1D * ScprodSW_;
  TH1D * ThdetaSW_;
  TH1D * M5SW_;
  TH1D * TTMS1SW_;
  TH1D * TTMS2SW_;
  TH1D * TTMS3SW_;

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
  TH1D * M6SSW_;
  TH1D * C6SSW_;
  TH1D * M8SSW_;
  TH1D * C8SSW_;
  TH1D * M45bestallSSW_;
  TH1D * Chi2extallSSW_;
  TH1D * DPbballSSW_;
  TH1D * SumHED4SSW_;
  TH1D * SumHPD4SSW_;
  TH1D * SumHED6SSW_;
  TH1D * SumHPD6SSW_;
  TH1D * HEDSSW_;
  TH1D * HPDSSW_;
  TH1D * Et6SSW_;
  TH1D * MwminSSW_;
  TH1D * HbestcombSSW_;
  TH1D * DrpairbestallSSW_;
  TH1D * M3aSSW_;
  TH1D * MwaSSW_;
  TH1D * ScprodSSW_;
  TH1D * ThdetaSSW_;
  TH1D * M5SSW_;
  TH1D * TTMS1SSW_;
  TH1D * TTMS2SSW_;
  TH1D * TTMS3SSW_;

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
  TH1D * M6SSSW_;
  TH1D * C6SSSW_;
  TH1D * M8SSSW_;
  TH1D * C8SSSW_;
  TH1D * M45bestallSSSW_;
  TH1D * Chi2extallSSSW_;
  TH1D * DPbballSSSW_;
  TH1D * SumHED4SSSW_;
  TH1D * SumHPD4SSSW_;
  TH1D * SumHED6SSSW_;
  TH1D * SumHPD6SSSW_;
  TH1D * HEDSSSW_;
  TH1D * HPDSSSW_;
  TH1D * Et6SSSW_;
  TH1D * MwminSSSW_;
  TH1D * HbestcombSSSW_;
  TH1D * DrpairbestallSSSW_;
  TH1D * M3aSSSW_;
  TH1D * MwaSSSW_;
  TH1D * ScprodSSSW_;
  TH1D * ThdetaSSSW_;
  TH1D * M5SSSW_;
  TH1D * TTMS1SSSW_;
  TH1D * TTMS2SSSW_;
  TH1D * TTMS3SSSW_;

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

  // Parameters of Et correction functions
  double tpar[7];
  double upar[7];

  double rel_lik;
  double njsss;

  double Nsltt_hjj;

  double total[10][5];
  double totalpass[10][5];
  double grandtotaltt[5];
  double grandtotalh[10];
  double grandgrandtotal;
  double grandtotalttpass[5];
  double grandtotalhpass[10];
  double grandgrandtotalpass;

  // Matrix with H kinematics
  // ------------------------
  double H[1000];
  double Hnot[1000];
  double Hread[1000];
  double Hnotread[1000];

  // Matrix with T kinematics
  // ------------------------
  double T[1000];
  double Tnot[1000];
  double Tread[1000];
  double Tnotread[1000];

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

  // Histograms with tag matrix
  // --------------------------
  TH1D * PTagEt_; 
  TH1D * PTagEta_;
  TH1D * PTagNt_; 

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

  // Select if the pixelTrigger must be done
  bool doPixelTrigger_;

  // Trigger, offline and total events
  int l1Eff_;
  int l1PixelEff_;
  int l1Njets_;
  int l1MEtSig_;
  int l1Tags_[4];
  int pixelNjets_;
  int pixelMEtSig_;
  int pixelTags_[4];
};
