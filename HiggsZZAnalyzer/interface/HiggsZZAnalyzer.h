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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

//#include "FWCore/Framework/interface/ESHandle.h"
//#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//#include "DataFormats/GeometryVector/interface/GlobalVector.h"
//#include "DataFormats/GeometryVector/interface/LocalVector.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
//#include "DataFormats/DetId/interface/DetId.h" 
//#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
//#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
//#include "DataFormats/GeometryVector/interface/LocalPoint.h"
//#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/TrackReco/interface/Track.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"


#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"

//
// class declaration
//
class HiggsZZAnalyzer : public edm::EDAnalyzer {

 public:
  explicit HiggsZZAnalyzer(const edm::ParameterSet&);
  ~HiggsZZAnalyzer();

  //
  // constants, enums and typedefs
  //

//
// static data member definitions
//

//
// constructors and destructor
//

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  unsigned nZ_;
  unsigned nZleptons_;
  unsigned nZjets_;
  unsigned nforwardjets_;
  unsigned njets_;

  int eventcounter_;
  int signalTopologySizeCounter_;
  int HtoZZeventcounter_;
  
  int muonCounter_;
  int goodMuonCounter_;
  int electronCounter_;
  int goodElectronCounter_;
  int goodIc5JetCounter_;
  int goodIc5centralJetCounter_;
  int goodIc5forwardJetCounter_;
  int Ic5forwardJetCounter_;
  std::vector<int> vecZtoParticleseventcounter_;

  int nbin_;
  double first_bin_et;
  double last_bin_et;
  double first_bin_eta;
  double last_bin_eta;
  double first_bin_phi;
  double last_bin_phi;

  // ----------member data ---------------------------
  edm::ParameterSet conf_;
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag combinedSVBJetTagsLabel_;
  double        leptonPtCut_;	          
  double        leptonIPsignificanceCut_;
  double        muonIso03SumPtCut_;     
  double        muonIso03emEtCut_;      
  double        centralJetEtCut_;
  double        centralJetEtaCut_;
  double        centralJetLeptonDeltaRCut_;
  double        centralJetEMfracCut_;
  double        forwardJetEtCut_;
  double        forwardJetEtaCut_;
  double        forwardJetLeptonDeltaRCut_;
  double        forwardDiJetMassCut_;
  double        pTbalanceCut_;
  double        ZdileptonMassminCut_;     
  double        ZdileptonMassmaxCut_;     
  double        ZdileptonDeltaRminCut_;
  double        ZdileptonDeltaRmaxCut_;
  double        ZdileptonShiftedEtaCut_;
  double        ZdijetMassminCut_;	
  double        ZdijetMassmaxCut_;	
  double        ZdijetDeltaRminCut_;	
  double        ZdijetDeltaRmaxCut_;	
  double        ZdijetShiftedEtaCut_;	
  double        ZZdeltaRminCut_;            
  double        ZZdeltaRmaxCut_;      
  std::string   MCParticleLabel_;      
  edm::InputTag simVtxLabel_;
  double        HMassCut_;

  TFile* OutputFile;
  // Use a dynamic construction,
  // or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.


  std::vector< std::pair<int,int> > vecTypeIndexPair;

  TH1D * eventsNumber_;
  TH1D * signalTopologySizeNumber_;
  TH1D * failedTopologyNumber_;
  std::vector<TH1D> vecZtoParticlesEventsNumber_;    
  TH1D * HtoZZEventsNumber_;

  TH1D * muonNumber_;
  TH1D * goodMuonNumber_;
  TH1D * electronNumber_;
  TH1D * goodElectronNumber_;
  TH1D * OfflinejetNumber_;
  TH1D * goodIc5JetNumber_;
  TH1D * goodIc5centralJetNumber_;
  TH1D * goodIc5forwardJetNumber_;
  TH1D * Ic5forwardJetNumber_;

  TH1D * OfflinejetUncorrEt_;
  TH1D * OfflinejetEt_;
  TH1D * OfflinejetEtScale_;
  TH2D * OfflinejetEtScaleVSUncorrEt_;
  TH2D * OfflinejetEtScaleVSEt_;
  TH2D * OfflinejetEtScaleVSEta_;
  TH2D * OfflinejetEtVSUncorrEt_;
  TH2D * OfflinejetEtVSEt_;
  TH2D * OfflinejetEtVSEta_;
  TH2D * OfflinejetUncorrEtVSUncorrEt_;
  TH2D * OfflinejetUncorrEtVSEt_;
  TH2D * OfflinejetUncorrEtVSEta_;
  TH1D * OfflinejetEMenergyfrac_;
  TH1D * OfflinejetHighEffDiscr_;
  TH1D * OfflinejetHighPurDiscr_;
  TH1D * OfflinejetMass_;
  TH1D * OfflinejetTKnumS2_;
  TH1D * OfflinejetTKsumPtS2_;
  TH1D * OfflinejetTKtagmassS2_;
  TH1D * OfflinejetPhi_;
  TH1D * OfflinejetEta_;
  TH1D * OfflinejetLeptonDeltaR_;
  TH1D * centralJetLeptonDeltaR_;

  TH1D * pTbalance_;
  TH1D * pTbalancePhi_;
  TH2D * pTbalanceVSMEt_;
  TH1D * pTbalanceMEtDeltaPhi_;

  TH1D * MuonNumber_;
  TH1D * MuonEta_;	 
  TH1D * MuonPhi_;	 
  TH1D * MuonPt_;	 
  TH1D * MuonPx_;	 
  TH1D * MuonPy_;	 
  TH1D * MuonPz_;	 
  TH1D * MuonEMcalo_;    
  TH1D * MuonHADcalo_;	 
  TH1D * MuonHOcalo_;	 
  TH1D * MuonIso03SumPt_;
  TH1D * MuonIso03emEt_;
  TH1D * MuonNormChi2_;	 
  TH1D * MuonHitsNum_;   
  TH1D * MuonImpactParamXY_;
  TH1D * MuonErrorImpactParamXY_;
  TH1D * MuonImpactParamSig_;

  TH1D * ElectronNumber_;
  TH1D * ElectronEt_;	    
  TH1D * ElectronEta_;	    
  TH1D * ElectronPhi_;	    
  TH1D * ElectronPt_;	    
  TH1D * ElectronPx_;	    
  TH1D * ElectronPy_;	    
  TH1D * ElectronPz_;	    
  TH1D * ElectronHadOverEm_;
  TH1D * ElectronIsoVal_;
  TH1D * ElectronImpactParamXY_;
  TH1D * ElectronErrorImpactParamXY_;
  TH1D * ElectronImpactParamSig_;

  TH1D * MEt_;
  TH1D * MEtSumEt_;
  TH1D * MEtSignificance_;
  TH1D * MEtPhi_;
  
  TH2D * dileptonMassVSdeltaR_;
  TH2D * dijetMassVSdeltaR_;   
  TH2D * dibjetMassVSdeltaR_;  
  TH1D * scaleddijetmass_;
  TH1D * dijetMassComparison_;
  TH2D * dijetMassVSscaleddibjetMass_;

  TH1D * recZllMass_;
  TH1D * recZjjMass_;
  TH1D * recZbbMass_;
  TH1D * recZllMassResolution_;
  TH1D * recZeeMassResolution_;
  TH1D * recZmmMassResolution_;
  TH1D * recZjjMassResolution_;
  TH1D * recZbbMassResolution_;
  TH1D * recZllDeltaR_;
  TH1D * recZjjDeltaR_;
  TH1D * recZbbDeltaR_;
  TH1D * recZllDeltaEta_;
  TH1D * recZjjDeltaEta_;
  TH1D * recZbbDeltaEta_;
  TH1D * recZllDeltaPhi_;
  TH1D * recZjjDeltaPhi_;
  TH1D * recZbbDeltaPhi_;
  TH1D * recZll1ShiftedEta_;
  TH1D * recZjj1ShiftedEta_;
  TH1D * recZbb1ShiftedEta_;
  TH1D * recZll2ShiftedEta_;
  TH1D * recZjj2ShiftedEta_;
  TH1D * recZbb2ShiftedEta_;
  TH1D * recZllShiftedEta_;
  TH1D * recZjjShiftedEta_;
  TH1D * recZbbShiftedEta_;

  TH1D * recZjjLeptonDeltaR_;
  TH1D * recZbbLeptonDeltaR_;

  TH1D * forwardJetEt_;
  TH1D * forwardJetLeptonDeltaR_;
  TH1D * forwardJetEta_;
  TH1D * forwardJet1Eta_;
  TH1D * forwardJet2Eta_;
  TH1D * forwardDiJetMass_;  
  TH1D * forwardDiJetDeltaPhi_;
  TH1D * forwardDiJetDeltaEta_;
  TH1D * forwardDiJetDeltaR_;

  TH1D * recHZZdeltaPhi_;
  TH1D * recHZZdeltaEta_;
  TH1D * recHZZdeltaR_;
  TH1D * recHlljj_;	     
  TH1D * recHlljj_forwardDiJetMass;	     
  TH1D * recHlljj_pTbalance;	     
  TH1D * recHlljj_ZdileptonMass;	     
  TH1D * recHlljj_leptonDeltaR;	     
  TH1D * recHlljj_leptonShiftedEta;	     
  TH1D * recHlljj_ZdijetMass;	     
  TH1D * recHlljj_ZjetDeltaR;	     
  TH1D * recHlljj_ZjetShiftedEta;	     
  TH2D * recHlljjVSforwardDiJetMass_;	     
  TH2D * recHlljjVSpTbalance_;	     
  TH2D * recHlljjVSZdileptonMass_;	     
  TH2D * recHlljjVSleptonDeltaR_;	     
  TH2D * recHlljjVSleptonShiftedEta_;	     
  TH2D * recHlljjVSZdijetMass_;	     
  TH2D * recHlljjVSZjetDeltaR_;	     
  TH2D * recHlljjVSZjetShiftedEta_;	     
  TH1D * recHlljj_AFTERforwardDiJetMass;	     
  TH1D * recHlljj_AFTERpTbalance;	     
  TH1D * recHlljj_AFTERZdileptonMass;	     
  TH1D * recHlljj_AFTERleptonDeltaR;	     
  TH1D * recHlljj_AFTERleptonShiftedEta;	     
  TH1D * recHlljj_AFTERZdijetMass;	     
  TH1D * recHlljj_AFTERZjetDeltaR;	     
  TH1D * recHlljj_AFTERZjetShiftedEta;	     
  TH1D * recHlljjMassResolution_;
  TH1D * recHllbbMassResolution_;

  TH1D * MHtoZZllqq_;	     
  TH1D * MHtoZZllbb_;	     
  TH1D * PtHParticles_;	     
  TH1D * PtHParticle1_;	     
  TH1D * PtHParticle2_;	     
  TH1D * DeltaPhiZZ_;
  TH1D * DeltaEtaZZ_;
  TH1D * DeltaRZZ_;

  TH2D * MZ1vsMZ2_;

  TH1D * Mbb_;
  TH1D * M4_;
  TH1D * M4irrBkg_;
  TH1D * MZZ_;

  TH2D * MZvsMbb_;

  std::vector<TH1D> vecMZtoParticles_;
  std::vector<TH1D> vecPtZParticles_;
  std::vector<TH1D> vecPtZParticle1_;
  std::vector<TH1D> vecPtZParticle2_;

  std::vector<TH1D> vecDeltaPhiZParticles_;
  std::vector<TH1D> vecDeltaEtaZParticles_;
  std::vector<TH1D> vecDeltaRZParticles_;

};
