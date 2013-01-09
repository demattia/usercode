#ifndef Onia2MuMu_h
#define Onia2MuMu_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h" 

#include "DataFormats/SiPixelDetId/interface/PixelModuleName.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include <memory>
#include <iostream>
#include <string>
#include <map>
#include <set>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TClonesArray.h>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/View.h"

using namespace reco;
using namespace edm;
using namespace std;

class Onia2MuMu : public edm::EDAnalyzer {
 public:
      explicit Onia2MuMu(const edm::ParameterSet&);
      ~Onia2MuMu();
 private:
      virtual void beginJob() ;
      virtual void beginRun(const Run & iRun, const EventSetup & iSetup);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void l1Report(const edm::Event &iEvent); 
      virtual void hltReport(const edm::Event &iEvent);
      virtual void fillGeneratorBlock(const edm::Event &iEvent);
      virtual void fillRecTracks(const edm::Event &iEvent);
      virtual void fillPhotons(const edm::Event &iEvent); 
      virtual void fillMuons(const edm::Event &iEvent);
      virtual void fillPATMuons(const edm::Event &iEvent);
      virtual void PAT_l1Report(const edm::Event &iEvent); 
      virtual void PAT_hltReport(const edm::Event &iEvent);
      virtual void fillBeamSpot(const edm::Event &iEvent);
      virtual void fillPrimaryVertex(const edm::Event &iEvent);
      virtual void fillPrimaryVertex2m(const edm::Event &iEvent);
      virtual void findOniaCategories(const edm::Event &iEvent);
      virtual void fillOniaMuMuTracks(TrackRef muon1, int m1, TrackRef muon2, int m2, TVector3 vperp2, int oniacato);

      double PhiInRange(const double& phi) const;
      template <class T, class U> double deltaR(const T & t, const U & u) const;
      double calcDeltaR(double eta1, double phi1,double eta2, double phi2);
      double GetTheta( TLorentzVector & a,  TLorentzVector & b) const;
      bool isAbHadron( int pdgID ) const;
      void printTrack(const reco::Track& muon) const;
      TLorentzVector lorentzMomentum(const reco::Track& muon) const;
      TLorentzVector lorentzMomentum(const pat::Muon& muon) const;
      TLorentzVector lorentzMomentum(const reco::PFCandidate& pfcl) const;
      TLorentzVector lorentzMomentumPi(const reco::Track& tr) const;
      TLorentzVector lorentzMomentumKa(const reco::Track& tr) const;
      // std::vector<unsigned int> muonStatHits(const reco::Track& tr);
      // std::vector<unsigned int> trackStatHits(const reco::Track& tr);
      std::vector<unsigned int> trackHits(const reco::Track& tr);
      TLorentzVector lorentzTriObj(const trigger::TriggerObject& muon) const;
      double invMass(const reco::Track& lhs, const reco::Track& rhs) const;

      string theOutFileName;               // Filename
      int runNb, eventNb, lumiBlock;       // run number, event number and lumi section
      int theOniaType;                     // 443 for Jpsi, 553 for Upsilon, etc
      int theDebugLevel;                   // 0 no prints, 1 some, 2 lots
      edm::InputTag thegenParticlesLabel;
      // edm::InputTag theStandAloneMuonsLabel;      // Muon information from standalone muons
      edm::InputTag theGlobalMuonsLabel;          // Muon information from global muons
      edm::InputTag theMuonsLabel;                // From this one can get both global and track info
      edm::InputTag theCaloMuonsLabel;             
      edm::InputTag theTrackLabel;                // Track information
      edm::InputTag thePhotonLabel;               // Photon information
      double thePhotonMinE;                       // Photon Mininum Energy
      edm::InputTag theBeamSpotLabel;            // Beam Spot
      edm::InputTag thePrimaryVertexLabel;        // Primary Vertex
      string thetriggerEventLabel;
      string theHLTriggerResults;                 // HLT trigger results
      string the8e29ProcName;
      string the1e31ProcName;
      edm::InputTag theL1GTReadoutRec;            // L1 trigger results
      edm::InputTag theL1MuonLabel;
      edm::InputTag thePATMuonsLabel;      

      bool theUseKinFit;                   // Yes or No to use of kinematic vertex fitting (T. Speer)
      bool theStoreGenFlag;                // Yes or No to store generator information
      bool theStoreHLTFlag;                // Yes or No to store HLT information
      bool theStoreL1Flag;                 // Yes or No to store L1 information
      bool theStoreTrkFlag;                // Yes or No to store track information (e.g. for hadr. act.)
      bool theStorePhotFlag;                // Yes or No to store photon information (e.g. for Chi_c)

      // bool theStoreSTAMuonFlag;
      bool theStoreGLBMuonFlag;
      bool theStoreTRKMuonFlag;
      bool theStoreCALMuonFlag;
      // bool theStoreAllMuonFlag;
      bool theStoreBeamSpotFlag; 
      bool theStorePriVtxFlag;             // Yes or No to store primary vertex
      bool theUsePrimaryNoMuons;           // Yes or No to remove muons from PV calculation
      bool theStoreOniaFlag;               // Yes or No to store Onium info
      bool theStoreChicFlag;               // Yes or No to store Onium + gamma combinations
      bool theStoreBpFlag;                 // Yes or No to store Onium + track combinations
      bool theStoreOniaRadiation;          // Yes or No to study Onium radiation
      bool theStoreWSOnia;                 // Yes or No to store wrong-sign mu-mu combinations
      bool theBeamSpotFlag;
      bool theminimumFlag;
      bool theAODFlag;
      bool theStorePATFlag;

      TFile* outFile;
      TTree *fTree;
      double oniaMass;
      double branch_ratio;      
      
      static const int NTRIGGERS = 5;
      static const int Max_track_size = 3000;
      static const int Max_QQ_size = 100;
      static const int Max_mu_size = 100;
      static const int Max_PriVtx_size = 20;
      static const unsigned int Max_trig_size = 10;

      int Mc_ProcessId;               // Process ID in PYTHIA (e.g. octet vs singlet)
      double Mc_EventScale;           // Pthat of process
      double Mc_EventWeight;          // Weight of event generated

      int Mc_QQ_size;                     // Number of Onia in event (usually 1)
      TClonesArray* Mc_QQ_4mom;       // Array of 4-momentum of Onium
      TClonesArray* Mc_QQ_3vec;       // Array of 3-d creation vertex of Onium
      TClonesArray* Mc_QQmoth_4mom;   // Array of 4-momentum of Onium mother
      TClonesArray* Mc_QQmoth_3vec;   // Array of 3-d creation vertex of Onium mother
      int Mc_QQmoth_id[20];           // particle id of Onium mother 
      int Mc_QQmupl_indx[20];         // Vector of index of mu+ from Onium in muon array
      int Mc_QQmumi_indx[20];         // Vector of index of mu- from Onium in muon array
      int Mc_mu_size;                     // Number of muons in event (usually 2 from 1 Onium)
      TClonesArray* Mc_mu_4mom;       // Array of 4-momentum of muons (all muons, could also be weak)
      TClonesArray* Mc_mu_3vec;       // Array of 3-d creation vertex of muons (all muons, could also be weak)
      int Mc_mu_id[20];               // ID (13 or -13) of muons (all muons, could also be weak)
      int Mc_mumoth_id[20];           // MotherID of muons (all muons, could also be weak)
      int Mc_chargedtrk_size;           // Amount of charged  particles in MC
      TClonesArray* Mc_chargedtrk_4mom; // Array of 4-vector of charged particles in MC
      int Mc_chargedtrk_charge[3000];   // Vector of charge of particles in MC
      
      int Reco_track_size;            // Number of reconstructed tracks
      TClonesArray* Reco_track_4mom;  // Array of 4-momentum of Reconstructed trackss
      TClonesArray* Reco_track_3vec;  // Array of 3-d creation vertex of Reconstructed trackss
      TClonesArray* Reco_track_CovM;  // Array of 5*5 covariance matrix
      double Reco_track_ptErr[3000];  // Vector of err on pt of global muons
      double Reco_track_phiErr[3000]; // Vector of err on phi of global muons
      double Reco_track_etaErr[3000]; // Vector of err on eta of global muons
      double Reco_track_d0[3000];     // Vector of d0 of tracks
      double Reco_track_d0err[3000];  // Vector of d0err of tracks
      double Reco_track_dz[3000];     // Vector of dz of tracks
      double Reco_track_dzerr[3000];  // Vector of dzerr of tracks
      int Reco_track_charge[3000]; // Vector of charge of tracks
      double Reco_track_chi2[3000];   // Vector of chi2 of tracks
      double Reco_track_ndof[3000];   // Vector of ndof of tracks
      int Reco_track_nhits[3000];  // Vector of nhits of tracks
   
      int Reco_gamma_size;            // Number of reconstructed gammas
      TClonesArray* Reco_gamma_4mom;  // Array of 4-momentum of Reconstructed gammas
      double Reco_gamma_phi[3000];  // phi 
      double Reco_gamma_eta[3000];  // eta

      int Reco_mu_glb_size;           // Number of reconstructed global muons
      TClonesArray* Reco_mu_glb_4mom; // Array of 4-momentum of Reconstructed global muons
      TClonesArray* Reco_mu_glb_track4mom; // Array of 4-momentum of Reconstructed global muons' track 
      TClonesArray* Reco_mu_glb_3vec; // Array of 3-d creation vertex of Reconstructed global muons
      // TClonesArray* Reco_mu_glb_CovM;
      double Reco_mu_glb_ptErr[200];   // Vector of err on pt of global muons
      double Reco_mu_glb_phiErr[200];  // Vector of err on phi of global muons
      double Reco_mu_glb_etaErr[200];  // Vector of err on eta of global muons
      double Reco_mu_glb_d0[200];      // Vector of d0 of global muons
      double Reco_mu_glb_d0err[200];   // Vector of d0err of global muons
      double Reco_mu_glb_dz[200];      // Vector of dz of global muons
      double Reco_mu_glb_dzerr[200];   // Vector of dzerr of global muons
      int Reco_mu_glb_charge[200];  // Vector of charge of global muons
      double Reco_mu_glb_normChi2[200];   // Vector of chi2/ndof of global muons
      // double Reco_mu_glb_ndof[200];   // Vector of ndof of global muons
      int Reco_mu_glb_nhitsCSC[200];    // Vector of number of valid hits of global muons
      int Reco_mu_glb_nhitsDT[200];    // Vector of number of valid hits of global muons
      int Reco_mu_glb_nhitstrack[200];    // Vector of number of valid hits of global muons
      double Reco_mu_glb_caloComp[200];    // Vector of calorimeter compatibilities
      double Reco_mu_glb_segmComp[200];    // Vector of muon segment compatibilities 
      double Reco_mu_glb_iso[200];    // Vector of isolations (NOW ONLY SUMPt OF TRACKS) 
      int Reco_mu_glb_nhitsStrip[200];  // Vectors of strip/pixel hits
      int Reco_mu_glb_nhitsPixB[200];
      int Reco_mu_glb_nhitsPixE[200];
      int Reco_mu_glb_nhitsPix1Hit[200];
      int Reco_mu_glb_nhitsPix1HitBE[200];


      int Reco_mu_trk_size;           // Number of reconstructed tracker muons
      TClonesArray* Reco_mu_trk_4mom; // Array of 4-momentum of Reconstructed tracker muons
      TClonesArray* Reco_mu_trk_3vec; // Array of 3-d creation vertex of Reconstructed tracker muons
      // TClonesArray* Reco_mu_trk_CovM;
      double Reco_mu_trk_ptErr[200];   // Vector of err on pt of tracker muons
      double Reco_mu_trk_phiErr[200];  // Vector of err on phi of tracker muons
      double Reco_mu_trk_etaErr[200];  // Vector of err on eta of tracker muons
      double Reco_mu_trk_d0[200];      // Vector of d0 of tracker muons
      double Reco_mu_trk_d0err[200];   // Vector of d0err of tracker muons
      double Reco_mu_trk_dz[200];      // Vector of dz of tracker muons
      double Reco_mu_trk_dzerr[200];   // Vector of dzerr of tracker muons
      int Reco_mu_trk_charge[200];  // Vector of charge of tracker muons
      double Reco_mu_trk_normChi2[200];   // Vector of chi2/ndof of tracker muons
      // double Reco_mu_trk_ndof[200];   // Vector of ndof of tracker muons
      // int Reco_mu_trk_nhitsCSC[200];    // Vector of number of valid hits of tracker muons
      // int Reco_mu_trk_nhitsDT[200];    // Vector of number of valid hits of tracker muons
      int Reco_mu_trk_nhitstrack[200];    // Vector of number of valid hits of tracker muons
      int Reco_mu_trk_nhitsStrip[200];  // Vectors of strip/pixel hits
      int Reco_mu_trk_nhitsPixB[200];
      int Reco_mu_trk_nhitsPixE[200];
      int Reco_mu_trk_nhitsPix1Hit[200];
      int Reco_mu_trk_nhitsPix1HitBE[200];

      int Reco_mu_trk_PIDmask[200];    // Bit word of TM selectors:
                                       // 1st bit : AllTrackerMuons = checks isTrackerMuon flag
                                       // 2nd bit : TrackerMuonArbitrated = resolve ambiguity of sharing segments
                                       // 3rd bit : TMLastStationLoose = penetration depth loose selector
                                       // 4th bit : TMLastStationTight = penetration depth tight selector 
                                       // 5th bit : TM2DCompatibilityLoose = likelihood based loose selector
                                       // 6th bit : TM2DCompatibilityTight = likelihood based tight selector
                                       // 7th bit : TMOneStationLoose = require one well matched segment
                                       // 8th bit : TMOneStationTight = require one well matched segment
                                       // 9th bit : TMLastStationOptimizedLowPtLoose = combination of TMLastStation and TMOneStation
                                       // 10th bit : TMLastStationOptimizedLowPtTight = combination of TMLastStation and TMOneStation     
      double Reco_mu_trk_caloComp[200];    // Vector of calorimeter compatibilities
      double Reco_mu_trk_segmComp[200];    // Vector of muon segment compatibilities 
      double Reco_mu_trk_iso[200];    // Vector of isolations (NOW ONLY SUMPt OF TRACKS) 

      int Reco_mu_cal_size;           // Number of reconstructed calo muons
      TClonesArray* Reco_mu_cal_4mom; // Array of 4-momentum of Reconstructed calo muons
      TClonesArray* Reco_mu_cal_3vec; // Array of 3-d creation vertex of Reconstructed calo muons
      // TClonesArray* Reco_mu_cal_CovM;
      double Reco_mu_cal_ptErr[200];   // Vector of err on pt of calo muons
      double Reco_mu_cal_phiErr[200];  // Vector of err on phi of calo muons
      double Reco_mu_cal_etaErr[200];  // Vector of err on eta of calo muons
      double Reco_mu_cal_d0[200];      // Vector of d0 of calo muons
      double Reco_mu_cal_d0err[200];   // Vector of d0err of calo muons
      double Reco_mu_cal_dz[200];      // Vector of dz of calo muons
      double Reco_mu_cal_dzerr[200];   // Vector of dzerr of calo muons
      int Reco_mu_cal_charge[200];  // Vector of charge of calo muons
      double Reco_mu_cal_normChi2[200];   // Vector of chi2/ndof of calo muons
      // double Reco_mu_cal_ndof[200];   // Vector of ndof of calo muons
      int Reco_mu_cal_nhitstrack[200];    // Vector of number of valid hits of calo muons
      int Reco_mu_cal_nhitsStrip[200];  // Vectors of strip/pixel hits
      int Reco_mu_cal_nhitsPixB[200];
      int Reco_mu_cal_nhitsPixE[200];
      int Reco_mu_cal_nhitsPix1Hit[200];
      int Reco_mu_cal_nhitsPix1HitBE[200];
      double Reco_mu_cal_caloComp[200];    // Vector of calorimeter compatibilities 

/////////////////// PAT muons
      
      int Pat_mu_glb_size;           // Number of reconstructed global muons
      TClonesArray* Pat_mu_glb_4mom; // Array of 4-momentum of Reconstructed global muons
      TClonesArray* Pat_mu_glb_3vec; // Array of 3-d creation vertex of Reconstructed global muons
      TClonesArray* Pat_mu_glb_CovM;
      double Pat_mu_glb_ptErr[200];   // Vector of err on pt of global muons
      double Pat_mu_glb_phiErr[200];  // Vector of err on phi of global muons
      double Pat_mu_glb_etaErr[200];  // Vector of err on eta of global muons
      double Pat_mu_glb_d0[200];      // Vector of d0 of global muons
      double Pat_mu_glb_d0err[200];   // Vector of d0err of global muons
      double Pat_mu_glb_dz[200];      // Vector of dz of global muons
      double Pat_mu_glb_dzerr[200];   // Vector of dzerr of global muons
      int Pat_mu_glb_charge[200];  // Vector of charge of global muons
      double Pat_mu_glb_chi2[200];   // Vector of chi2 of global muons
      double Pat_mu_glb_ndof[200];   // Vector of ndof of global muons
      int Pat_mu_glb_nhits[200];  // Vector of number of valid hits of global muons
      
      int Pat_mu_sta_size;           // Number of reconstructed stand alone muons
      TClonesArray* Pat_mu_sta_4mom; // Array of 4-momentum of Reconstructed stand alone muons
      TClonesArray* Pat_mu_sta_3vec; // Array of 3-d creation vertex of Reconstructed stand alone muons
      TClonesArray* Pat_mu_sta_CovM;
      double Pat_mu_sta_ptErr[200];   // Vector of err on pt of stand alone muons
      double Pat_mu_sta_phiErr[200];  // Vector of err on phi of stand alone muons
      double Pat_mu_sta_etaErr[200];  // Vector of err on eta of stand alone muons
      double Pat_mu_sta_d0[200];      // Vector of d0 of stand alone muons
      double Pat_mu_sta_d0err[200];   // Vector of d0err of stand alone muons
      double Pat_mu_sta_dz[200];      // Vector of dz of stand alone muons
      double Pat_mu_sta_dzerr[200];   // Vector of dzerr of stand alone muons
      int Pat_mu_sta_charge[200];  // Vector of charge of stand alone muons
      double Pat_mu_sta_chi2[200];   // Vector of chi2 of stand alone muons
      double Pat_mu_sta_ndof[200];   // Vector of ndof of stand alone muons
      int Pat_mu_sta_nhits[200];  // Vector of number of valid hits of stand alone muons
      
      int Pat_mu_trk_size;           // Number of reconstructed tracker muons
      TClonesArray* Pat_mu_trk_4mom; // Array of 4-momentum of Reconstructed tracker muons
      TClonesArray* Pat_mu_trk_3vec; // Array of 3-d creation vertex of Reconstructed tracker muons
      TClonesArray* Pat_mu_trk_CovM;
      double Pat_mu_trk_ptErr[200];   // Vector of err on pt of tracker muons
      double Pat_mu_trk_phiErr[200];  // Vector of err on phi of tracker muons
      double Pat_mu_trk_etaErr[200];  // Vector of err on eta of tracker muons
      double Pat_mu_trk_d0[200];      // Vector of d0 of tracker muons
      double Pat_mu_trk_d0err[200];   // Vector of d0err of tracker muons
      double Pat_mu_trk_dz[200];      // Vector of dz of tracker muons
      double Pat_mu_trk_dzerr[200];   // Vector of dzerr of tracker muons
      int Pat_mu_trk_charge[200];  // Vector of charge of tracker muons
      double Pat_mu_trk_chi2[200];   // Vector of chi2 of tracker muons
      double Pat_mu_trk_ndof[200];   // Vector of ndof of tracker muons
      int Pat_mu_trk_nhits[200];  // Vector of number of valid hits of tracker muons
      
      int Pat_mu_cal_size;           // Number of reconstructed calorimeter muons
      TClonesArray* Pat_mu_cal_4mom; // Array of 4-momentum of Reconstructed calorimeter muons
      TClonesArray* Pat_mu_cal_3vec; // Array of 3-d creation vertex of Reconstructed calorimeter muons
      TClonesArray* Pat_mu_cal_CovM;
      double Pat_mu_cal_ptErr[200];   // Vector of err on pt of calorimeter muons
      double Pat_mu_cal_phiErr[200];  // Vector of err on phi of calorimeter muons
      double Pat_mu_cal_etaErr[200];  // Vector of err on eta of calorimeter muons
      double Pat_mu_cal_d0[200];      // Vector of d0 of calorimeter muons
      double Pat_mu_cal_d0err[200];   // Vector of d0err of calorimeter muons
      double Pat_mu_cal_dz[200];      // Vector of dz of calorimeter muons
      double Pat_mu_cal_dzerr[200];   // Vector of dzerr of calorimeter muons
      int Pat_mu_cal_charge[200];  // Vector of charge of calorimeter muons
      double Pat_mu_cal_chi2[200];   // Vector of chi2 of calorimeter muons
      double Pat_mu_cal_ndof[200];   // Vector of ndof of calorimeter muons
      int Pat_mu_cal_nhits[200];  // Vector of number of valid hits of calorimeter muons             
                  
//////////////////////////////////////////      

      /* int Reco_mu_sta_size;           // Number of reconstructed standalone muons
      TClonesArray* Reco_mu_sta_4mom; // Array of 4-momentum of Reconstructed standalone muons
      TClonesArray* Reco_mu_sta_3vec; // Array of 3-d creation vertex of Reconstructed standalone muons
      TClonesArray* Reco_mu_sta_CovM;
      double Reco_mu_sta_ptErr[200];   // Vector of err on pt of standalone muons
      double Reco_mu_sta_phiErr[200];  // Vector of err on phi of standalone muons
      double Reco_mu_sta_etaErr[200];  // Vector of err on eta of standalone muons
      double Reco_mu_sta_d0[200];      // Vector of d0 of standalone muons
      double Reco_mu_sta_d0err[200];   // Vector of d0err of standalone muons
      double Reco_mu_sta_dz[200];      // Vector of dz of standalone muons
      double Reco_mu_sta_dzerr[200];   // Vector of dzerr of standalone muons
      int Reco_mu_sta_charge[200];   // Vector of charge of standalone muons
      double Reco_mu_sta_chi2[200];   // Vector of chi2 of standalone muons
      double Reco_mu_sta_ndof[200];   // Vector of ndof of standalone muons
      int Reco_mu_sta_nhits[200];  // Vector of number of valid hits of standalone muons */

      int Reco_QQ_size;           // Number of reconstructed Onia 
      int Reco_QQ_type[3000];     // Onia category:
                                  // 0 = golden (2 global muons)
                                  // 1 = silver (1 global - 1 tracker muon)
                                  // 2 = bronze (2 tracker muons)
                                  // 3 = copper (1 global - 1 calo muon)
                                  // 4 = iron (1 tracker - 1 calo muon)
                                  // 5 = crap (2 calo muons)
      int theOniaMaxCat;          // Maximum category above which to store cands
      bool theSkimOnOniaMaxCat;   // Skim on the maximum category above which to store cands
      int Reco_QQ_sign[3000];     // Mu Mu combinations sign:
                                  // 0 = +/- (signal)
                                  // 1 = +/+
                                  // 2 = -/-
      TClonesArray* Reco_QQ_4mom; // Array of 4-momentum of Reconstructed onia
      int Reco_QQ_mupl[3000];       // Index of muon plus in onia 
      int Reco_QQ_mumi[3000];       // Index of muon minus in onia
      int Reco_QQ_mulpt[3000];      // Index of lower-pT muon in onia 
      int Reco_QQ_muhpt[3000];      // Index of higher-pT minus in onia
      double Reco_QQ_DeltaR[3000];  // DeltaR of the two muons
      double Reco_QQ_cosTheta[3000];// Polarization angle 
      double Reco_QQ_s[3000];       // S : the sum of the to muons impact parameter significance 
      bool Reco_QQ_VtxIsVal[3000];  // Vertex is valid or not  
      TClonesArray* Reco_QQ_Vtx;  // 3-d Vertex 
      double Reco_QQ_VxxE[3000];     // errors of x
      double Reco_QQ_VyyE[3000];     // errors of y
      double Reco_QQ_VzzE[3000];     // errors of 
      double Reco_QQ_VyxE[3000];     // corr errors of x and y
      double Reco_QQ_VzxE[3000];     // corr errors of z and x
      double Reco_QQ_VzyE[3000];     // corr errors of z and y
      double Reco_QQ_lxy[3000];     // Decay length
      double Reco_QQ_lxyErr[3000];  // Decay length errors
      double Reco_QQ_normChi2[3000];// Normalized chi2 of vertex fitting 
      double Reco_QQ_probChi2[3000];// chi2 probability of vertex fitting 
      double Reco_QQ_cosAlpha[3000];// Alpha: the angle of lxy and onia moemtum
      double Reco_QQ_ctau[3000];    // ctau: flying time

      int Reco_Chic_size;           // Number of reconstructed chi_c
      TClonesArray* Reco_Chic_4mom; // Array of 4-momentum of Reconstructed chi_c
      int Reco_Chic_OniaDaug[3000];   // Index of onia in chi_c 
      int Reco_Chic_GammaDaug[3000];  // Index of gamma in chi_c
      double Reco_Chic_DeltaM[3000];  // M(onia-gamma) - M(onia)
      
      int Reco_Bp_size;           // Number of reconstructed Onia from global muons
      TClonesArray* Reco_Bp_4mom; // Array of 4-momentum of Reconstructed Bp
      int Reco_Bp_OniaDaug[3000];   // Index of onia in Bp 
      int Reco_Bp_KDaug[3000];      // Index of kaon in Bp
      bool Reco_Bp_VtxIsVal[3000];  // Vertex is valid or not  
      TClonesArray* Reco_Bp_Vtx;  // 3-d Vertex 
      double Reco_Bp_VxE[3000];     // errors of x
      double Reco_Bp_VyE[3000];     // errors of y
      double Reco_Bp_VzE[3000];     // errors of z
      double Reco_Bp_lxy[3000];     // Decay length
      double Reco_Bp_lxyErr[3000];  // Decay length errors
      double Reco_Bp_normChi2[3000];// Normalized chi2 of vertex fitting 
      double Reco_Bp_cosAlpha[3000];// Alpha: the angle of lxy and onia moemtum
      double Reco_Bp_ctau[3000];    // ctau: flying time

      double Reco_BeamSpot_x;
      double Reco_BeamSpot_y;
      double Reco_BeamSpot_z;
      double Reco_BeamSpot_xxE;
      double Reco_BeamSpot_yyE;
      double Reco_BeamSpot_zzE;
      double Reco_BeamSpot_yxE;
      double Reco_BeamSpot_zyE;
      double Reco_BeamSpot_zxE;

      int Reco_PriVtx_size;                // Number of reconstructed primary vertex
      TClonesArray* Reco_PriVtx_3vec;      // 3-d vector of primary vertex
      double Reco_PriVtx_xxE[100];           // X error
      double Reco_PriVtx_yyE[100];           // Y error
      double Reco_PriVtx_zzE[100];           // Z error
      double Reco_PriVtx_yxE[100];           // X error
      double Reco_PriVtx_zyE[100];           // Y error
      double Reco_PriVtx_zxE[100];           // Z error
      int Reco_PriVtx_trkSize[100];      // Number of tracks in this primary vertex
      double Reco_PriVtx_chi2[100];         // chi2 of primary vertex 
      double Reco_PriVtx_ndof[100];         // number of freedom degree 

      int Reco_PriVtx_1st_trkSize;        // Number of tracks in the first primary vertex
      int Reco_PriVtx_1st_trkindex[3000]; //index for the tracks in the first primary vertex

      int L1TBits_size;               // Number of L1 trigger bits
      bool L1TBits_accept[3000];       // L1 trigger bits
      bool L1TGlobal_Decision;        // L1 trigger global decision
    
      int L1_mu_size;
      TClonesArray* L1_mu_4mom;
      int L1_mu_charge[200];

      int HLTBits_size;               // Number of HLT trigger bits 
      bool HLTBits_wasrun[1000];       // Each HLT bits was run or not
      bool HLTBits_accept[1000];       // Each HLT bits fired or not
      bool HLTBits_error[1000];        // Each HLT bits run successfully or not 
      bool HLTGlobal_wasrun;          // The HLT was run or not
      bool HLTGlobal_Decision;        // Global HLT decision
      bool HLTGlobal_error;           // HLT path error or not
  
      int HLT1Mu3_L2_size;
      TClonesArray* HLT1Mu3_L2_4mom;
      int HLT1Mu3_L2_id[20];
      int HLT1Mu3_L3_size;
      TClonesArray* HLT1Mu3_L3_4mom;
      int HLT1Mu3_L3_id[20];

      int HLT1Mu5_L2_size;
      TClonesArray* HLT1Mu5_L2_4mom;
      int HLT1Mu5_L2_id[20];
      int HLT1Mu5_L3_size;
      TClonesArray* HLT1Mu5_L3_4mom;
      int HLT1Mu5_L3_id[20];

      int HLT1Mu9_L2_size;
      TClonesArray* HLT1Mu9_L2_4mom;
      int HLT1Mu9_L2_id[20];
      int HLT1Mu9_L3_size;
      TClonesArray* HLT1Mu9_L3_4mom;
      int HLT1Mu9_L3_id[20];

      int HLT1Mu11_L2_size;
      TClonesArray* HLT1Mu11_L2_4mom;
      int HLT1Mu11_L2_id[20];
      int HLT1Mu11_L3_size;
      TClonesArray* HLT1Mu11_L3_4mom;
      int HLT1Mu11_L3_id[20];

      int HLT2Mu0_L2_size;
      TClonesArray* HLT2Mu0_L2_4mom;
      int HLT2Mu0_L2_id[20];
      int HLT2Mu0_L3_size;
      TClonesArray* HLT2Mu0_L3_4mom;
      int HLT2Mu0_L3_id[20];

      int HLT2IsoMu3_L2_size;
      TClonesArray* HLT2IsoMu3_L2_4mom;
      int HLT2IsoMu3_L2_id[20];
      int HLT2IsoMu3_L3_size;
      TClonesArray* HLT2IsoMu3_L3_4mom;
      int HLT2IsoMu3_L3_id[20];

      int HLT2Mu3_L2_size;
      TClonesArray* HLT2Mu3_L2_4mom;
      int HLT2Mu3_L2_id[20];
      int HLT2Mu3_L3_size;
      TClonesArray* HLT2Mu3_L3_4mom;
      int HLT2Mu3_L3_id[20];

      int HLTJpsi2Mu_L2_size;
      TClonesArray* HLTJpsi2Mu_L2_4mom;
      int HLTJpsi2Mu_L2_id[20];
      int HLTJpsi2Mu_L3_size;
      TClonesArray* HLTJpsi2Mu_L3_4mom;
      int HLTJpsi2Mu_L3_id[20];

      int HLTUpsilon2Mu_L2_size;
      TClonesArray* HLTUpsilon2Mu_L2_4mom;
      int HLTUpsilon2Mu_L2_id[20];
      int HLTUpsilon2Mu_L3_size;
      TClonesArray* HLTUpsilon2Mu_L3_4mom;
      int HLTUpsilon2Mu_L3_id[20];

      InputTag hltModules[2][NTRIGGERS]; // in order: L2, L3
      std::vector<std::string> hltPaths;
      int hltBits[NTRIGGERS];

      HLTConfigProvider hltConfig;
      edm::ESHandle<TransientTrackBuilder> theB;

      // physics objects (RECO)
      edm::Handle<reco::PFCandidateCollection> pfAll;
      reco::PFCandidateCollection pfClusters;
      edm::Handle<TrackCollection> allTracks;
      reco::TrackCollection noMuonTracks;
      reco::MuonCollection theGlobalMuons;
      reco::MuonCollection theTrkMuons;
      reco::CaloMuonCollection theCaloMuons;  

      unsigned int fNevt;            // event number
      unsigned int maxCatToStoreChic;
      unsigned int maxCatToStoreBp;
};

#endif
