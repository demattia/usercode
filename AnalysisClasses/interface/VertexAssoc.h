// -*- C++ -*-
//
// Package:    VertexAssoc
// Class:      VertexAssoc
//
// Code extracted from Roberto.cc, 02/04/08 R.Casagrande
//
/**\class VertexAssoc VertexAssoc.cc AnalysisExamples/VertexAssoc/src/VertexAssoc.cc

Description: <one line class summary>

*/
//
// $Id: VertexAssoc.h
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

//=====ROBERTO===================
//-------------------------
//tracce e vertici
//-------------------------

// For the SimVertex
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

//For the SimTracks
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Data include files
// ------------------

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"

//=====ROBERTO==============================

// Associator for the jets
// -----------------------
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"


// Class declaration
// -----------------
class VertexAssoc : public edm::EDAnalyzer {
public:
  explicit VertexAssoc(const edm::ParameterSet&);
  ~VertexAssoc();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
 
  int eventcounter_;

  // Declare as static so that only one exists, even if more
  // than one VertexAssoc object is created
  // -------------------------------------------------------
  
  edm::ParameterSet conf_;
  TFile* OutputFile;

  // Declare here, since it does not have a default constructor
  // it will be initialized with an initialization list ( in
  // the VertexAssoc constructor ).
  //  HiVariables HiVar;

  // Use a dynamic construction, or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.
 
  
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag MCParticleLabel_;
  
  //======ROBERTO==============
  
  //------------------------
  //tracce e vertici
  //------------------------
  edm::InputTag simVtxLabel_;
  edm::InputTag simTkLabel_;
  std::string impactParameterTagInfos;

 


  //======ROBERTO==============
 
  bool QCD_;
  std::string OutputEffFileName;

  //=====ROBERTO============
   
 
  //Histograms for primary vertex and PU vertices
  TH1D * H_NPUVtx_;
 
  //Histograms for Z PV vertices
  TH1D * H_primaryVtxZ_;
  //Histograms for Z PU vertices  
  TH1D * H_PUVtxZ_;
  //Histograms for Z simTracks
  TH1D * H_Z_simTks_;
  //Histograms for Z reco tracks
  TH1D * H_TksZ_ ;
  TH1D * H_TksZ_Error_; 

  //Histogram Z jet 
  TH1D * H_TksZ_Weighted_Avg_;
  TH1D * H_TksZ_WAvg_Error_;
  
  //Histograms selected tracks  
  TH1D * H_SelTk_;
  
  //Histograms for 2D, 3D IP Significance
  TH1D * H_S2DTks_;
  TH1D * H_S3DTks_;
  TH2D * H_Tks_Z_S3D_;

  //Histograms for Njet  
  TH1D * H_NJet_withTracks_;
  TH1D * H_NJet_withNoTracks_;
  TH1D * H_NJet_fromPV_;
  TH1D * H_NJet_fromPU_;
  TH1D * H_NJet_NotAssigned_;
  
  //Histograms for Dz jet - vertices 
  TH1D * H_Vmin_;
  TH1D * H_Dz_WAvg_Jet_Vtx_tot_;
  TH1D * H_Dz_WAvg_Jet_Vtx_2Best_tot_; 
  TH2D * H_Dz_WAvg_Jet_Vtx_vs_2Best_;

  //Histograms for Dz vertices
  TH1D * H_Dz_Vtx_;
  
  //Histogram for Dz jet - vertices 
  TProfile * H_Dz_WAvg_Vtx_vs_Nj_;
  
  TProfile * H_Z_WAvg_Vtx_vs_Nj_;
  
  //Histogram (phases space) for Dz jet - PU vertices
  TProfile *  H_Profile_WAvg_DzVtx_;

  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
  TH1D * H_Z_WAvg_Jet_Ntracks_[13];
 
  //Histograms for eta, phi reco tracks
  TH1D *     H_TksEta_;
  TH1D *    H_TksPhi_;
  //Histograms for Pt, Chi2 reco tracks
  TH1D *   H_TksPt_;
  TH1D *    H_TksChi2_;

  //Histograms for Et, eta, phi of calojets
  TH1D * H_Et_Jets_;
  TH1D * H_eta_Jets_;
  TH1D * H_phi_Jets_;
  
  //Histograms Dphi, Deta, Delta2 jet - reco selected tracks
  TH1D * H_Dphi_Jet_Tk_;
  TH1D * H_Deta_Jet_Tk_;
  TH1D * H_D2_Jet_Tk_;

  //Histograms for eta, phi sim tracks
  TH1D * H_Phi_simTk_;
  TH1D * H_Eta_simTk_;

  //Histograms for Delta2 minimum simTrack - reco selected track
  TH1D * H_D2_Tk_sim_reco_; 
  TH1D * H_D2_Tk_sim_reco_2Best_;  
  TH2D * H_D2_Tk_vs_2Best_;
  
  //Histograms of Dz simTrack - reco selected track
  TH1D * H_Dz_Tk_sim_reco_;
  TH1D * H_Z_sim_jet_;
  TH1D * H_Dz_jet_sim_reco_;
  TH1D * H_Dz_jet_sim_reco_avg_ ;

  //Histogram of Dz selected sim tracks
  TH1D *  H_Dz_sim_selected_tk_;
  //TProfile max rate of same Z in selected sim tracks in function of the # of tracks 
  TProfile * H_Profile_Ntk_vs_sameZ_;
  
  //Histograms of rate of same selected sim tracks Z and 
  //rate of same selected sim tracks Z vs # of tracks in the jet
  TH1D * H_Frequenze_;
  TH2D * H_Ntk_vs_SameZ_;
  
  //Histogram of selected simulated tracks' Z
  TH1D *  H_Z_selected_sim_tracks_;
  
  //Histograms of Dz jet - vertex: |Z(true) - Z(vtx)| 
  TH1D * H_Dz_sim_Jet_Vtx_;
  TH1D * H_Dz_sim_Jet_Vtx_2Best_;
  TH2D * H_Dz_sim_Jet_Vtx_vs_2Best_;
  
  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
  TH1D * H_Dz_sim_Jet_Vtx_NTk_[13];
  
  //Histogram of |Z(vera) - Z(reco)| jet
  TH1D * H_Dz_jet_true_reco_;

  //Histograms of Dz vera jet - Z sim vertex associated
  TH1D *  H_Dz_sim_jet_simVtx_;         
  TH1D *  H_Dz_sim_jet_simVtx_2Best_;   
  TH2D *  H_Dz_sim_jet_simVtx_vs_2Best_;
  
  //Histograms of Z jet and error after the rejection algorithm (5 sigma)
  TH1D *  H_ZJet_RJAlg_;
  TH1D *  H_ZErrorJet_RJAlg_;
  TH1D *  H_DzAfterRejection_;
  
  //Histograms of dZ:  Z jet( weight average after rejection of i sigma ) - vertices 
  TH1D * H_dZJet_RJAlg_SFunc_[20];
  
  //Histogram of Ntk in function # sigma rejection
  TProfile * H_Ntk_vs_NS_;
 
  TH1D *  H_NJetRJAlgFromPV_;
  TH1D *  H_NJetRJAlgFromPU_;
  
  TProfile * H_ProfileDzVtxVsNjRJAlg_;
 
  //=====ROBERTO============

  float loose_;
  float medium_;
  float tight_;

};
