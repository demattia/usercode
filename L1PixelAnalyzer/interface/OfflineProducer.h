// -*- C++ -*-
//
// Package:    OfflineProducer
// Class:      OfflineProducer
// 
/**\class OfflineProducer OfflineProducer.cc AnalysisExamples/OfflineProducer/src/OfflineProducer.cc

 Description: <one line class summary>

 Implementation:
 This class shows how to access:
 - level 1 calorimetric quantities
 - offline corrected jets (calibration performed here)
 - offline corrected MET, depending on jets corrections
 - MC informations
 - B tagging

 Evaluates:
 - association of MC partons to offline jets
 - association of btags to offline jets

*/
//
// Original Author:  Marco De Mattia
//         Created:  Tue May  8 13:05:37 CEST 2007
// $Id: OfflineProducer.h,v 1.6 2008/05/28 15:04:49 demattia Exp $
//
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

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

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Calo MEt
#include "DataFormats/METReco/interface/CaloMET.h"

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"

// Multiple TH1F
#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"

//
// class declaration
//

class OfflineProducer : public edm::EDProducer {
 public:
  explicit OfflineProducer(const edm::ParameterSet&);
  ~OfflineProducer();

  /// Method to evaluate the invariant mass coming from the tracks used to evalute the tag
  double tagTracksInvariantMass( const reco::TrackRefVector & selectedTracks ) const;

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//   double PI_;

  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one OfflineProducer object is created
  static L1Trig L1Trigger;

  edm::ParameterSet conf_;

  std::string CaloJetAlgorithm, L2JetCorrectionService, L3JetCorrectionService, METCollection;
  std::string genParticleCandidates;
  std::string trackCountingHighEffJetTags;
  std::string trackCountingHighPurJetTags;
  std::string impactParameterTagInfos;

  edm::InputTag paramGlobalMuons_;

  edm::InputTag electronCandidates_;
  edm::InputTag electronHcalIsolation_;

  edm::InputTag tauTagInfo_;

  unsigned int numTkCut;
  std::string OutputEffFileName;

  // Trigger efficiency counters
  // Multijet
  int Eff_;
  int Eff_et1_;
  int Eff_et2_;
  int Eff_et3_;
  int Eff_et4_;

  // MEt+Jet
  int Eff_MEtJet_;

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

  // Collections
  // L1
  std::string cenJets_;
  std::string forJets_;
  std::string tauJets_;
  std::string l1MEt_;
  // Offline
  std::string offlineJets_;
  std::string offlineMEt_;
  std::string simpleTracks_;
  std::string MCParticles_;
  std::string globalMuons_;
  std::string simpleElectrons_;
  std::string simpleTaus_;
  // Summary
  std::string summary_;
  unsigned int eventType_;

  bool doL1Trig_;

  // TagTracksInvarianMass
  double chargedpi_mass_;

  // ----------member data ---------------------------
};
