///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Package:    Onia2MuMu
// Class:      Onia2MuMu
// 
// Class: Onia2MuMu Onia2MuMu.cc Onia2MuMu/Onia2MuMu/src/Onia2MuMu.cc
//
// Description: Analyzer for Onia->MuMu events
//
//
//
// Original Author:   Aafke Kraan, Zongchang Yang
//          Created:  Mon Nov 19 10:24:55 CET 2007
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#include "HeavyFlavorAnalysis/Onia2MuMu/interface/Onia2MuMu.h"

////////////////////////////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////////////////////////////
Onia2MuMu::Onia2MuMu(const edm::ParameterSet& iConfig)
{
  theOutFileName             = iConfig.getParameter<string>("OutputFileName");
  theOniaType                = iConfig.getParameter<int>("OniaType");
  theOniaMaxCat              = iConfig.getParameter<int>("OniaMaxCat");
  theSkimOnOniaMaxCat        = iConfig.getParameter<bool>("skimOnOniaMaxCat");
  theDebugLevel              = iConfig.getParameter<int>("DebugLevel");
  thegenParticlesLabel       = iConfig.getParameter<edm::InputTag>("genParticlesLabel");  
  // theStandAloneMuonsLabel    = iConfig.getParameter<edm::InputTag>("StandAloneMuonsLabel");
  theGlobalMuonsLabel        = iConfig.getParameter<edm::InputTag>("GlobalMuonsLabel");
  theMuonsLabel              = iConfig.getParameter<edm::InputTag>("MuonsLabel");
  theCaloMuonsLabel          = iConfig.getParameter<edm::InputTag>("CaloMuonsLabel");
  theTrackLabel              = iConfig.getParameter<edm::InputTag>("TrackLabel");
  thePhotonLabel             = iConfig.getParameter<edm::InputTag>("PhotonLabel");
  thePhotonMinE              = iConfig.getParameter<double>("PhotonMinEnergy");
  theBeamSpotLabel           = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
  thePrimaryVertexLabel      = iConfig.getParameter<edm::InputTag>("PrimaryVertexLabel");
  thetriggerEventLabel       = iConfig.getParameter<string>("triggerEventLabel");
  theHLTriggerResults        = iConfig.getParameter<string>("triggerResultsLabel");
  the8e29ProcName            = iConfig.getParameter<string>("HLTprocessName8e29");
  the1e31ProcName            = iConfig.getParameter<string>("HLTprocessName1e31");
  theL1GTReadoutRec          = iConfig.getParameter<edm::InputTag>("L1GTReadoutRec");
  theL1MuonLabel             = iConfig.getParameter<edm::InputTag>("L1MuonLabel");
  thePATMuonsLabel           = iConfig.getParameter<edm::InputTag>("PATMuonsLabel");  
  theUsePrimaryNoMuons       = iConfig.getParameter<bool>("UsePrimaryNoMuons");
  theUseKinFit               = iConfig.getParameter<bool>("useKinFit");
  theStoreGenFlag            = iConfig.getParameter<bool>("StoreGenFlag");
  theStoreHLTFlag            = iConfig.getParameter<bool>("StoreHLTFlag");
  theStoreL1Flag             = iConfig.getParameter<bool>("StoreL1Flag");
  theStoreTrkFlag            = iConfig.getParameter<bool>("StoreTrkFlag");
  theStorePhotFlag           = iConfig.getParameter<bool>("StorePhotonFlag");
  // theStoreSTAMuonFlag        = iConfig.getParameter<bool>("StoreSTAMuonFlag");
  theStoreGLBMuonFlag        = iConfig.getParameter<bool>("StoreGLBMuonFlag");
  theStoreTRKMuonFlag        = iConfig.getParameter<bool>("StoreTRKMuonFlag");
  theStoreCALMuonFlag        = iConfig.getParameter<bool>("StoreCALMuonFlag");
  // theStoreAllMuonFlag        = iConfig.getParameter<bool>("StoreAllMuonFlag");
  theStoreBeamSpotFlag       = iConfig.getParameter<bool>("StoreBeamSpotFlag");
  theStorePriVtxFlag         = iConfig.getParameter<bool>("StorePriVtxFlag");
  theStoreOniaFlag           = iConfig.getParameter<bool>("StoreOniaFlag");
  theStoreChicFlag           = iConfig.getParameter<bool>("StoreChicFlag");
  theStoreBpFlag             = iConfig.getParameter<bool>("StoreBpFlag");
  theStoreOniaRadiation      = iConfig.getParameter<bool>("StoreOniaRadiation");
  theStoreWSOnia             = iConfig.getParameter<bool>("StoreWSOnia");
  theBeamSpotFlag            = iConfig.getParameter<bool>("UsingBeamSpot");
  theminimumFlag             = iConfig.getParameter<bool>("minimumFlag");
  theAODFlag                 = iConfig.getParameter<bool>("UsingAOD");
  theStorePATFlag            = iConfig.getParameter<bool>("StorePATFlag");
  fNevt=0;
}

////////////////////////////////////////////////////////////////////////
// destructor
////////////////////////////////////////////////////////////////////////
Onia2MuMu::~Onia2MuMu()
{
  //
}


////////////////////////////////////////////////////////////////////////
// called at beginning
////////////////////////////////////////////////////////////////////////
void Onia2MuMu::beginJob()
{

  outFile = new TFile(theOutFileName.c_str(), "RECREATE", "");
  outFile->cd();
  fTree = new TTree ("T1", "CMSSW Quarkonia tree");

  maxCatToStoreChic = 0;
  maxCatToStoreBp = 1;

  if ( theOniaType==443 )  {
    oniaMass=3.09688;
    branch_ratio = 0.0593;
  }
  // else if ( theOniaType==23 ) {
  //  oniaMass=91.1876;
  //  branch_ratio = 0.03366;
  // }
  else if ( theOniaType==553 ) {
    oniaMass=9.46030;
    branch_ratio = 0.0248;
  }
  else if ( theOniaType==100443 ) {
    oniaMass=3.68600;
    branch_ratio = 0.0083;
  }
  else if ( theOniaType==100553 ) {
    oniaMass=10.0233;
    branch_ratio = 0.014;
  }
  else if ( theOniaType==200553 ) {
    oniaMass=10.3552;
    branch_ratio = 0.0218;
  }
  else {
    cout<<"Please set the correct onia type: 443 (Jpsi), 100443(psi'), 553(Upsilon1S), 100553(Upsilon2S), 200553(Upsilon3S),"<<endl;
    // return;
  }

  Mc_QQ_4mom=new TClonesArray("TLorentzVector", 10000);
  Mc_QQ_3vec=new TClonesArray("TVector3", 10000);
  Mc_QQmoth_4mom=new TClonesArray("TLorentzVector", 10000);
  Mc_QQmoth_3vec=new TClonesArray("TVector3", 10000);
  Mc_mu_4mom=new TClonesArray("TLorentzVector", 10000);
  Mc_mu_3vec=new TClonesArray("TVector3", 10000);
  Mc_chargedtrk_4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_track_4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_track_3vec=new TClonesArray("TVector3", 10000);
  Reco_gamma_4mom=new TClonesArray("TLorentzVector", 10000);
  // Reco_gamma_3pos=new TClonesArray("TVector3", 10000);
  Reco_track_CovM = new TClonesArray("TMatrixD", 10000);
  Reco_mu_glb_4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_mu_glb_track4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_mu_glb_3vec=new TClonesArray("TVector3", 10000);
  Reco_mu_trk_4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_mu_trk_3vec=new TClonesArray("TVector3", 10000); 
  Reco_mu_cal_4mom=new TClonesArray("TLorentzVector", 10000);
  Reco_mu_cal_3vec=new TClonesArray("TVector3", 10000); 
  // Reco_mu_glb_CovM = new TClonesArray("TMatrixD", 10000);
  Pat_mu_glb_4mom=new TClonesArray("TLorentzVector", 10000);
  Pat_mu_glb_3vec=new TClonesArray("TVector3", 10000);
  Pat_mu_glb_CovM = new TClonesArray("TMatrixD", 10000);
  Pat_mu_sta_4mom=new TClonesArray("TLorentzVector", 10000);
  Pat_mu_sta_3vec=new TClonesArray("TVector3", 10000);
  Pat_mu_sta_CovM = new TClonesArray("TMatrixD", 10000);
  Pat_mu_trk_4mom=new TClonesArray("TLorentzVector", 10000);
  Pat_mu_trk_3vec=new TClonesArray("TVector3", 10000);
  Pat_mu_trk_CovM = new TClonesArray("TMatrixD", 10000);
  Pat_mu_cal_4mom=new TClonesArray("TLorentzVector", 10000);
  Pat_mu_cal_3vec=new TClonesArray("TVector3", 10000);
  Pat_mu_cal_CovM = new TClonesArray("TMatrixD", 10000);
  // Reco_mu_sta_4mom=new TClonesArray("TLorentzVector", 10000);
  // Reco_mu_sta_3vec=new TClonesArray("TVector3", 10000);
  // Reco_mu_sta_CovM = new TClonesArray("TMatrixD", 10000);
  Reco_QQ_4mom=new TClonesArray("TLorentzVector",10000);
  Reco_QQ_Vtx =new TClonesArray("TVector3",10000);
  Reco_Chic_4mom=new TClonesArray("TLorentzVector",10000);
  Reco_Bp_4mom=new TClonesArray("TLorentzVector",10000);
  Reco_Bp_Vtx =new TClonesArray("TVector3",10000);
  Reco_PriVtx_3vec =new TClonesArray("TVector3",10000);
  L1_mu_4mom = new TClonesArray("TLorentzVector",10000);
  // HLT_mu_L2_4mom = new TClonesArray("TLorentzVector",10000);
  // HLT_mu_L3_4mom = new TClonesArray("TLorentzVector",10000);
  // HLT1Mu3_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT1Mu3_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT1Mu5_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT1Mu5_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT1Mu9_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT1Mu9_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT1Mu11_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT1Mu11_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT2Mu0_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT2Mu0_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT2IsoMu3_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT2IsoMu3_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLT2Mu3_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLT2Mu3_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLTJpsi2Mu_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLTJpsi2Mu_L3_4mom=new TClonesArray("TLorentzVector",10000);
  // HLTUpsilon2Mu_L2_4mom=new TClonesArray("TLorentzVector",10000);
  HLTUpsilon2Mu_L3_4mom=new TClonesArray("TLorentzVector",10000);

  fTree->Branch("eventNb",             &eventNb,             "eventNb/I");
  fTree->Branch("runNb",               &runNb,               "runNb/I");
  fTree->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I"); 

  if(theStoreGenFlag){
    fTree->Branch("Mc_ProcessId",        &Mc_ProcessId,        "Mc_ProcessId/I");
    fTree->Branch("Mc_EventScale",       &Mc_EventScale,       "Mc_EventScale/D");
    fTree->Branch("Mc_EventWeight",      &Mc_EventWeight,      "Mc_EventWeight/D");
    fTree->Branch("Mc_QQ_size",          &Mc_QQ_size,          "Mc_QQ_size/I");
    fTree->Branch("Mc_QQ_4mom",          "TClonesArray",       &Mc_QQ_4mom, 32000, 0);
    fTree->Branch("Mc_QQ_3vec",          "TClonesArray",       &Mc_QQ_3vec, 32000, 0);
   
    if(!theminimumFlag) {
      fTree->Branch("Mc_QQmoth_4mom",      "TClonesArray",       &Mc_QQmoth_4mom, 32000, 0);
      fTree->Branch("Mc_QQmoth_3vec",      "TClonesArray",       &Mc_QQmoth_3vec, 32000, 0);
      fTree->Branch("Mc_QQmoth_id",         Mc_QQmoth_id,        "Mc_QQmoth_id[Mc_QQ_size]/I");
    }
    fTree->Branch("Mc_QQmupl_indx",       Mc_QQmupl_indx,      "Mc_QQmupl_indx[Mc_QQ_size]/I");
    fTree->Branch("Mc_QQmumi_indx",       Mc_QQmumi_indx,      "Mc_QQmumi_indx[Mc_QQ_size]/I");
    fTree->Branch("Mc_mu_size",          &Mc_mu_size,          "Mc_mu_size/I");
    fTree->Branch("Mc_mu_4mom",          "TClonesArray",       &Mc_mu_4mom, 32000, 0);
    fTree->Branch("Mc_mu_3vec",          "TClonesArray",       &Mc_mu_3vec, 32000, 0);
    fTree->Branch("Mc_mu_id",             Mc_mu_id,            "Mc_mu_id[Mc_mu_size]/I");
    fTree->Branch("Mc_mumoth_id",         Mc_mumoth_id,        "Mc_mumoth_id[Mc_mu_size]/I");
  }

  if(theStoreTrkFlag){
    fTree->Branch("Reco_track_size",     &Reco_track_size,     "Reco_track_size/I");
    fTree->Branch("Reco_track_4mom",     "TClonesArray",       &Reco_track_4mom, 32000, 0);
    fTree->Branch("Reco_track_3vec",     "TClonesArray",       &Reco_track_3vec, 32000, 0);
    if(!theminimumFlag) { 
      fTree->Branch("Reco_track_CovM",     "TClonesArray",       &Reco_track_CovM, 32000, 0);
      fTree->Branch("Reco_track_phiErr",    Reco_track_phiErr,   "Reco_track_phiErr[Reco_track_size]/D");
      fTree->Branch("Reco_track_etaErr",    Reco_track_etaErr,   "Reco_track_etaErr[Reco_track_size]/D");
      fTree->Branch("Reco_track_ptErr",     Reco_track_ptErr,    "Reco_track_ptErr[Reco_track_size]/D");
      fTree->Branch("Reco_track_d0",        Reco_track_d0,       "Reco_track_d0[Reco_track_size]/D");
      fTree->Branch("Reco_track_d0err",     Reco_track_d0err,    "Reco_track_d0err[Reco_track_size]/D");
      fTree->Branch("Reco_track_dz",        Reco_track_dz,       "Reco_track_dz[Reco_track_size]/D");
      fTree->Branch("Reco_track_dzerr",     Reco_track_dzerr,    "Reco_track_dzerr[Reco_track_size]/D");
    }
    fTree->Branch("Reco_track_charge",    Reco_track_charge,   "Reco_track_charge[Reco_track_size]/I");
    fTree->Branch("Reco_track_chi2",      Reco_track_chi2,     "Reco_track_chi2[Reco_track_size]/D");
    fTree->Branch("Reco_track_ndof",      Reco_track_ndof,     "Reco_track_ndof[Reco_track_size]/D");
    fTree->Branch("Reco_track_nhits",     Reco_track_nhits,    "Reco_track_nhits[Reco_track_size]/I");
  }

  if(theStorePhotFlag){
    fTree->Branch("Reco_gamma_size",     &Reco_gamma_size,     "Reco_gamma_size/I");
    fTree->Branch("Reco_gamma_4mom",     "TClonesArray",       &Reco_gamma_4mom, 32000, 0);
    fTree->Branch("Reco_gamma_phi",   Reco_gamma_phi,  "Reco_gamma_phi[Reco_gamma_size]/D");
    fTree->Branch("Reco_gamma_eta",   Reco_gamma_eta,  "Reco_gamma_eta[Reco_gamma_size]/D");  
  }

  if(theStoreGLBMuonFlag){
    fTree->Branch("Reco_mu_glb_size",    &Reco_mu_glb_size,    "Reco_mu_glb_size/I");
    fTree->Branch("Reco_mu_glb_4mom",    "TClonesArray",       &Reco_mu_glb_4mom, 32000, 0);
    fTree->Branch("Reco_mu_glb_track4mom",    "TClonesArray",       &Reco_mu_glb_track4mom, 32000, 0);
    fTree->Branch("Reco_mu_glb_3vec",    "TClonesArray",       &Reco_mu_glb_3vec, 32000, 0);
    if(!theminimumFlag) {
      // fTree->Branch("Reco_mu_glb_CovM",     "TClonesArray",       &Reco_mu_glb_CovM, 32000, 0);
      fTree->Branch("Reco_mu_glb_phiErr",   Reco_mu_glb_phiErr,  "Reco_mu_glb_phiErr[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_etaErr",   Reco_mu_glb_etaErr,  "Reco_mu_glb_etaErr[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_ptErr",    Reco_mu_glb_ptErr,   "Reco_mu_glb_ptErr[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_d0",       Reco_mu_glb_d0,      "Reco_mu_glb_d0[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_d0err",    Reco_mu_glb_d0err,   "Reco_mu_glb_d0err[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_dz",       Reco_mu_glb_dz,      "Reco_mu_glb_dz[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_dzerr",    Reco_mu_glb_dzerr,   "Reco_mu_glb_dzerr[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_normChi2",     Reco_mu_glb_normChi2,    "Reco_mu_glb_normChi2[Reco_mu_glb_size]/D");
    // fTree->Branch("Reco_mu_glb_ndof",     Reco_mu_glb_ndof,    "Reco_mu_glb_ndof[Reco_mu_glb_size]/D");
      fTree->Branch("Reco_mu_glb_nhitstrack",    Reco_mu_glb_nhitstrack,   "Reco_mu_glb_nhitstrack[Reco_mu_glb_size]/I");      
      fTree->Branch("Reco_mu_glb_nhitsStrip",    Reco_mu_glb_nhitsStrip,   "Reco_mu_glb_nhitsStrip[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsPixB",    Reco_mu_glb_nhitsPixB,   "Reco_mu_glb_nhitsPixB[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsPixE",    Reco_mu_glb_nhitsPixE,   "Reco_mu_glb_nhitsPixE[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsPix1Hit",    Reco_mu_glb_nhitsPix1Hit,   "Reco_mu_glb_nhitsPix1Hit[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsPix1HitBE",    Reco_mu_glb_nhitsPix1HitBE,   "Reco_mu_glb_nhitsPix1HitBE[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsDT",    Reco_mu_glb_nhitsDT,   "Reco_mu_glb_nhitsDT[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_nhitsCSC",    Reco_mu_glb_nhitsCSC,   "Reco_mu_glb_nhitsCSC[Reco_mu_glb_size]/I");
      fTree->Branch("Reco_mu_glb_caloComp",   Reco_mu_glb_caloComp,  "Reco_mu_glb_caloComp[Reco_mu_glb_size]/D"); 
      fTree->Branch("Reco_mu_glb_segmComp",   Reco_mu_glb_segmComp,  "Reco_mu_glb_segmComp[Reco_mu_glb_size]/D"); 
      fTree->Branch("Reco_mu_glb_iso",   Reco_mu_glb_iso,  "Reco_mu_glb_iso[Reco_mu_glb_size]/D");  
    }

    fTree->Branch("Reco_mu_glb_charge",   Reco_mu_glb_charge,  "Reco_mu_glb_charge[Reco_mu_glb_size]/I"); 
  }   

  if(theStoreTRKMuonFlag){
    fTree->Branch("Reco_mu_trk_size",    &Reco_mu_trk_size,    "Reco_mu_trk_size/I");
    fTree->Branch("Reco_mu_trk_4mom",    "TClonesArray",       &Reco_mu_trk_4mom, 32000, 0);
    fTree->Branch("Reco_mu_trk_3vec",    "TClonesArray",       &Reco_mu_trk_3vec, 32000, 0);
    if(!theminimumFlag) {
      // fTree->Branch("Reco_mu_trk_CovM",     "TClonesArray",       &Reco_mu_trk_CovM, 32000, 0);
      fTree->Branch("Reco_mu_trk_phiErr",   Reco_mu_trk_phiErr,  "Reco_mu_trk_phiErr[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_etaErr",   Reco_mu_trk_etaErr,  "Reco_mu_trk_etaErr[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_ptErr",    Reco_mu_trk_ptErr,   "Reco_mu_trk_ptErr[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_d0",       Reco_mu_trk_d0,      "Reco_mu_trk_d0[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_d0err",    Reco_mu_trk_d0err,   "Reco_mu_trk_d0err[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_dz",       Reco_mu_trk_dz,      "Reco_mu_trk_dz[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_dzerr",    Reco_mu_trk_dzerr,   "Reco_mu_trk_dzerr[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_normChi2",     Reco_mu_trk_normChi2,    "Reco_mu_trk_normChi2[Reco_mu_trk_size]/D");
    // fTree->Branch("Reco_mu_trk_ndof",     Reco_mu_trk_ndof,    "Reco_mu_trk_ndof[Reco_mu_trk_size]/D");
      fTree->Branch("Reco_mu_trk_nhitstrack",    Reco_mu_trk_nhitstrack,   "Reco_mu_trk_nhitstrack[Reco_mu_trk_size]/I");      
      fTree->Branch("Reco_mu_trk_nhitsStrip",    Reco_mu_trk_nhitsStrip,   "Reco_mu_trk_nhitsStrip[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_nhitsPixB",    Reco_mu_trk_nhitsPixB,   "Reco_mu_trk_nhitsPixB[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_nhitsPixE",    Reco_mu_trk_nhitsPixE,   "Reco_mu_trk_nhitsPixE[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_nhitsPix1Hit",    Reco_mu_trk_nhitsPix1Hit,   "Reco_mu_trk_nhitsPix1Hit[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_nhitsPix1HitBE",    Reco_mu_trk_nhitsPix1HitBE,   "Reco_mu_trk_nhitsPix1HitBE[Reco_mu_trk_size]/I");
      // fTree->Branch("Reco_mu_trk_nhitsDT",    Reco_mu_trk_nhitsDT,   "Reco_mu_trk_nhitsDT[Reco_mu_trk_size]/I");
      // fTree->Branch("Reco_mu_trk_nhitsCSC",    Reco_mu_trk_nhitsCSC,   "Reco_mu_trk_nhitsCSC[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_PIDmask",    Reco_mu_trk_PIDmask,   "Reco_mu_trk_PIDmask[Reco_mu_trk_size]/I");
      fTree->Branch("Reco_mu_trk_caloComp",   Reco_mu_trk_caloComp,  "Reco_mu_trk_caloComp[Reco_mu_trk_size]/D"); 
      fTree->Branch("Reco_mu_trk_segmComp",   Reco_mu_trk_segmComp,  "Reco_mu_trk_segmComp[Reco_mu_trk_size]/D"); 
      fTree->Branch("Reco_mu_trk_iso",   Reco_mu_trk_iso,  "Reco_mu_trk_iso[Reco_mu_trk_size]/D");  
    }

    fTree->Branch("Reco_mu_trk_charge",   Reco_mu_trk_charge,  "Reco_mu_trk_charge[Reco_mu_trk_size]/I"); 
  }   
  
  if(theStoreCALMuonFlag){
    fTree->Branch("Reco_mu_cal_size",    &Reco_mu_cal_size,    "Reco_mu_cal_size/I");
    fTree->Branch("Reco_mu_cal_4mom",    "TClonesArray",       &Reco_mu_cal_4mom, 32000, 0);
    fTree->Branch("Reco_mu_cal_3vec",    "TClonesArray",       &Reco_mu_cal_3vec, 32000, 0);
    if(!theminimumFlag) {
      // fTree->Branch("Reco_mu_cal_CovM",     "TClonesArray",       &Reco_mu_cal_CovM, 32000, 0);
      fTree->Branch("Reco_mu_cal_phiErr",   Reco_mu_cal_phiErr,  "Reco_mu_cal_phiErr[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_etaErr",   Reco_mu_cal_etaErr,  "Reco_mu_cal_etaErr[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_ptErr",    Reco_mu_cal_ptErr,   "Reco_mu_cal_ptErr[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_d0",       Reco_mu_cal_d0,      "Reco_mu_cal_d0[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_d0err",    Reco_mu_cal_d0err,   "Reco_mu_cal_d0err[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_dz",       Reco_mu_cal_dz,      "Reco_mu_cal_dz[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_dzerr",    Reco_mu_cal_dzerr,   "Reco_mu_cal_dzerr[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_normChi2",     Reco_mu_cal_normChi2,    "Reco_mu_cal_normChi2[Reco_mu_cal_size]/D");
    // fTree->Branch("Reco_mu_cal_ndof",     Reco_mu_cal_ndof,    "Reco_mu_cal_ndof[Reco_mu_cal_size]/D");
      fTree->Branch("Reco_mu_cal_nhitstrack",    Reco_mu_cal_nhitstrack,   "Reco_mu_cal_nhitstrack[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_nhitsStrip",    Reco_mu_cal_nhitsStrip,   "Reco_mu_cal_nhitsStrip[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_nhitsPixB",    Reco_mu_cal_nhitsPixB,   "Reco_mu_cal_nhitsPixB[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_nhitsPixE",    Reco_mu_cal_nhitsPixE,   "Reco_mu_cal_nhitsPixE[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_nhitsPix1Hit",    Reco_mu_cal_nhitsPix1Hit,   "Reco_mu_cal_nhitsPix1Hit[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_nhitsPix1HitBE",    Reco_mu_cal_nhitsPix1HitBE,   "Reco_mu_cal_nhitsPix1HitBE[Reco_mu_cal_size]/I");
      fTree->Branch("Reco_mu_cal_caloComp",   Reco_mu_cal_caloComp,  "Reco_mu_cal_caloComp[Reco_mu_cal_size]/D");   
    }

    fTree->Branch("Reco_mu_cal_charge",   Reco_mu_cal_charge,  "Reco_mu_cal_charge[Reco_mu_cal_size]/I"); 
  }   
  
  if(theStorePATFlag){
  
    fTree->Branch("Pat_mu_glb_size",    &Pat_mu_glb_size,    "Pat_mu_glb_size/I");
    fTree->Branch("Pat_mu_glb_4mom",    "TClonesArray",       &Pat_mu_glb_4mom, 32000, 0);
    fTree->Branch("Pat_mu_glb_3vec",    "TClonesArray",       &Pat_mu_glb_3vec, 32000, 0);
    if(!theminimumFlag) {
      fTree->Branch("Pat_mu_glb_CovM",     "TClonesArray",       &Pat_mu_glb_CovM, 32000, 0);
      fTree->Branch("Pat_mu_glb_phiErr",   Pat_mu_glb_phiErr,  "Pat_mu_glb_phiErr[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_etaErr",   Pat_mu_glb_etaErr,  "Pat_mu_glb_etaErr[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_ptErr",    Pat_mu_glb_ptErr,   "Pat_mu_glb_ptErr[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_d0",       Pat_mu_glb_d0,      "Pat_mu_glb_d0[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_d0err",    Pat_mu_glb_d0err,   "Pat_mu_glb_d0err[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_dz",       Pat_mu_glb_dz,      "Pat_mu_glb_dz[Pat_mu_glb_size]/D");
      fTree->Branch("Pat_mu_glb_dzerr",    Pat_mu_glb_dzerr,   "Pat_mu_glb_dzerr[Pat_mu_glb_size]/D");
    }
    fTree->Branch("Pat_mu_glb_charge",   Pat_mu_glb_charge,  "Pat_mu_glb_charge[Pat_mu_glb_size]/I");
    fTree->Branch("Pat_mu_glb_chi2",     Pat_mu_glb_chi2,    "Pat_mu_glb_chi2[Pat_mu_glb_size]/D");
    fTree->Branch("Pat_mu_glb_ndof",     Pat_mu_glb_ndof,    "Pat_mu_glb_ndof[Pat_mu_glb_size]/D");
    fTree->Branch("Pat_mu_glb_nhits",    Pat_mu_glb_nhits,   "Pat_mu_glb_nhits[Pat_mu_glb_size]/I");
    
    fTree->Branch("Pat_mu_sta_size",    &Pat_mu_sta_size,    "Pat_mu_sta_size/I");
    fTree->Branch("Pat_mu_sta_4mom",    "TClonesArray",       &Pat_mu_sta_4mom, 32000, 0);
    fTree->Branch("Pat_mu_sta_3vec",    "TClonesArray",       &Pat_mu_sta_3vec, 32000, 0);
    if(!theminimumFlag) {
      fTree->Branch("Pat_mu_sta_CovM",     "TClonesArray",       &Pat_mu_sta_CovM, 32000, 0);
      fTree->Branch("Pat_mu_sta_phiErr",   Pat_mu_sta_phiErr,  "Pat_mu_sta_phiErr[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_etaErr",   Pat_mu_sta_etaErr,  "Pat_mu_sta_etaErr[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_ptErr",    Pat_mu_sta_ptErr,   "Pat_mu_sta_ptErr[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_d0",       Pat_mu_sta_d0,      "Pat_mu_sta_d0[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_d0err",    Pat_mu_sta_d0err,   "Pat_mu_sta_d0err[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_dz",       Pat_mu_sta_dz,      "Pat_mu_sta_dz[Pat_mu_sta_size]/D");
      fTree->Branch("Pat_mu_sta_dzerr",    Pat_mu_sta_dzerr,   "Pat_mu_sta_dzerr[Pat_mu_sta_size]/D");
    }
    fTree->Branch("Pat_mu_sta_charge",   Pat_mu_sta_charge,  "Pat_mu_sta_charge[Pat_mu_sta_size]/I");
    fTree->Branch("Pat_mu_sta_chi2",     Pat_mu_sta_chi2,    "Pat_mu_sta_chi2[Pat_mu_sta_size]/D");
    fTree->Branch("Pat_mu_sta_ndof",     Pat_mu_sta_ndof,    "Pat_mu_sta_ndof[Pat_mu_sta_size]/D");
    fTree->Branch("Pat_mu_sta_nhits",    Pat_mu_sta_nhits,   "Pat_mu_sta_nhits[Pat_mu_sta_size]/I");
    
    fTree->Branch("Pat_mu_trk_size",    &Pat_mu_trk_size,    "Pat_mu_trk_size/I");
    fTree->Branch("Pat_mu_trk_4mom",    "TClonesArray",       &Pat_mu_trk_4mom, 32000, 0);
    fTree->Branch("Pat_mu_trk_3vec",    "TClonesArray",       &Pat_mu_trk_3vec, 32000, 0);
    if(!theminimumFlag) {
      fTree->Branch("Pat_mu_trk_CovM",     "TClonesArray",       &Pat_mu_trk_CovM, 32000, 0);
      fTree->Branch("Pat_mu_trk_phiErr",   Pat_mu_trk_phiErr,  "Pat_mu_trk_phiErr[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_etaErr",   Pat_mu_trk_etaErr,  "Pat_mu_trk_etaErr[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_ptErr",    Pat_mu_trk_ptErr,   "Pat_mu_trk_ptErr[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_d0",       Pat_mu_trk_d0,      "Pat_mu_trk_d0[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_d0err",    Pat_mu_trk_d0err,   "Pat_mu_trk_d0err[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_dz",       Pat_mu_trk_dz,      "Pat_mu_trk_dz[Pat_mu_trk_size]/D");
      fTree->Branch("Pat_mu_trk_dzerr",    Pat_mu_trk_dzerr,   "Pat_mu_trk_dzerr[Pat_mu_trk_size]/D");
    }
    fTree->Branch("Pat_mu_trk_charge",   Pat_mu_trk_charge,  "Pat_mu_trk_charge[Pat_mu_trk_size]/I");
    fTree->Branch("Pat_mu_trk_chi2",     Pat_mu_trk_chi2,    "Pat_mu_trk_chi2[Pat_mu_trk_size]/D");
    fTree->Branch("Pat_mu_trk_ndof",     Pat_mu_trk_ndof,    "Pat_mu_trk_ndof[Pat_mu_trk_size]/D");
    fTree->Branch("Pat_mu_trk_nhits",    Pat_mu_trk_nhits,   "Pat_mu_trk_nhits[Pat_mu_trk_size]/I");
    
    fTree->Branch("Pat_mu_cal_size",    &Pat_mu_cal_size,    "Pat_mu_cal_size/I");
    fTree->Branch("Pat_mu_cal_4mom",    "TClonesArray",       &Pat_mu_cal_4mom, 32000, 0);
    fTree->Branch("Pat_mu_cal_3vec",    "TClonesArray",       &Pat_mu_cal_3vec, 32000, 0);
    if(!theminimumFlag) {
      fTree->Branch("Pat_mu_cal_CovM",     "TClonesArray",       &Pat_mu_cal_CovM, 32000, 0);
      fTree->Branch("Pat_mu_cal_phiErr",   Pat_mu_cal_phiErr,  "Pat_mu_cal_phiErr[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_etaErr",   Pat_mu_cal_etaErr,  "Pat_mu_cal_etaErr[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_ptErr",    Pat_mu_cal_ptErr,   "Pat_mu_cal_ptErr[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_d0",       Pat_mu_cal_d0,      "Pat_mu_cal_d0[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_d0err",    Pat_mu_cal_d0err,   "Pat_mu_cal_d0err[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_dz",       Pat_mu_cal_dz,      "Pat_mu_cal_dz[Pat_mu_cal_size]/D");
      fTree->Branch("Pat_mu_cal_dzerr",    Pat_mu_cal_dzerr,   "Pat_mu_cal_dzerr[Pat_mu_cal_size]/D");
    }
    fTree->Branch("Pat_mu_cal_charge",   Pat_mu_cal_charge,  "Pat_mu_cal_charge[Pat_mu_cal_size]/I");
    fTree->Branch("Pat_mu_cal_chi2",     Pat_mu_cal_chi2,    "Pat_mu_cal_chi2[Pat_mu_cal_size]/D");
    fTree->Branch("Pat_mu_cal_ndof",     Pat_mu_cal_ndof,    "Pat_mu_cal_ndof[Pat_mu_cal_size]/D");
    fTree->Branch("Pat_mu_cal_nhits",    Pat_mu_cal_nhits,   "Pat_mu_cal_nhits[Pat_mu_cal_size]/I");         
  }   
 
  // if (theStoreSTAMuonFlag) { 
    /*  fTree->Branch("Reco_mu_sta_size",    &Reco_mu_sta_size,    "Reco_mu_sta_size/I");
    fTree->Branch("Reco_mu_sta_4mom",    "TClonesArray",       &Reco_mu_sta_4mom, 32000, 0);
    fTree->Branch("Reco_mu_sta_3vec",    "TClonesArray",       &Reco_mu_sta_3vec, 32000, 0);
    if(!theminimumFlag) {
      fTree->Branch("Reco_mu_sta_CovM",    "TClonesArray",       &Reco_mu_sta_CovM, 32000, 0);
      fTree->Branch("Reco_mu_sta_phiErr",   Reco_mu_sta_phiErr,  "Reco_mu_sta_phiErr[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_etaErr",   Reco_mu_sta_etaErr,  "Reco_mu_sta_etaErr[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_ptErr",    Reco_mu_sta_ptErr,   "Reco_mu_sta_ptErr[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_d0",       Reco_mu_sta_d0,      "Reco_mu_sta_d0[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_d0err",    Reco_mu_sta_d0err,   "Reco_mu_sta_d0err[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_dz",       Reco_mu_sta_dz,      "Reco_mu_sta_dz[Reco_mu_sta_size]/D");
      fTree->Branch("Reco_mu_sta_dzerr",    Reco_mu_sta_dzerr,   "Reco_mu_sta_dzerr[Reco_mu_sta_size]/D");
    }
    fTree->Branch("Reco_mu_sta_charge",   Reco_mu_sta_charge,  "Reco_mu_sta_charge[Reco_mu_sta_size]/I");
    fTree->Branch("Reco_mu_sta_chi2",     Reco_mu_sta_chi2,    "Reco_mu_sta_chi2[Reco_mu_sta_size]/D");
    fTree->Branch("Reco_mu_sta_ndof",     Reco_mu_sta_ndof,    "Reco_mu_sta_ndof[Reco_mu_sta_size]/D");
    fTree->Branch("Reco_mu_sta_nhits",    Reco_mu_sta_nhits,   "Reco_mu_sta_nhits[Reco_mu_sta_size]/I"); */
  // }   

  /* if (theStoreAllMuonFlag){
    fTree->Branch("Reco_mu_size",      &Reco_mu_size,     "Reco_mu_size/I");
    fTree->Branch("Reco_mu_Normsize",  &Reco_mu_Normsize, "Reco_mu_Normsize/I");
    fTree->Branch("Reco_mu_Calmsize",  &Reco_mu_Calmsize, "Reco_mu_Calmsize/I"); 
    fTree->Branch("Reco_mu_links_glb", Reco_mu_links_glb, "Reco_mu_links_glb[Reco_mu_size]/I");
    fTree->Branch("Reco_mu_links_sta", Reco_mu_links_sta, "Reco_mu_links_sta[Reco_mu_size]/I");
    fTree->Branch("Reco_mu_links_trk", Reco_mu_links_trk, "Reco_mu_links_trk[Reco_mu_size]/I");
    fTree->Branch("Reco_mu_is_sta",    Reco_mu_is_sta,    "Reco_mu_is_sta[Reco_mu_size]/B");
    fTree->Branch("Reco_mu_is_glb",    Reco_mu_is_glb,    "Reco_mu_is_glb[Reco_mu_size]/B");
    fTree->Branch("Reco_mu_is_trk",    Reco_mu_is_trk,    "Reco_mu_is_trk[Reco_mu_size]/B");
    fTree->Branch("Reco_mu_is_cal",    Reco_mu_is_cal,    "Reco_mu_is_cal[Reco_mu_size]/B");
    fTree->Branch("Reco_mu_caloComp",  Reco_mu_caloComp,  "Reco_mu_caloComp[Reco_mu_size]/D"); 
    } */

  if(theStoreOniaFlag){
    fTree->Branch("Reco_QQ_size",    &Reco_QQ_size,    "Reco_QQ_size/I");
    fTree->Branch("Reco_QQ_type",     Reco_QQ_type,    "Reco_QQ_type[Reco_QQ_size]/I");
    fTree->Branch("Reco_QQ_4mom",    "TClonesArray",       &Reco_QQ_4mom, 32000, 0);
    fTree->Branch("Reco_QQ_mupl",     Reco_QQ_mupl,    "Reco_QQ_mupl[Reco_QQ_size]/I");
    fTree->Branch("Reco_QQ_mumi",     Reco_QQ_mumi,    "Reco_QQ_mumi[Reco_QQ_size]/I");
    fTree->Branch("Reco_QQ_mulpt",    Reco_QQ_mulpt,   "Reco_QQ_mulpt[Reco_QQ_size]/I");
    fTree->Branch("Reco_QQ_muhpt",    Reco_QQ_muhpt,   "Reco_QQ_muhpt[Reco_QQ_size]/I");
    fTree->Branch("Reco_QQ_DeltaR",   Reco_QQ_DeltaR,  "Reco_QQ_DeltaR[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_cosTheta", Reco_QQ_cosTheta,"Reco_QQ_cosTheta[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_s",        Reco_QQ_s,       "Reco_QQ_s[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VtxIsVal", Reco_QQ_VtxIsVal,"Reco_QQ_VtxIsVal[Reco_QQ_size]/B");
    fTree->Branch("Reco_QQ_Vtx",     "TClonesArray",       &Reco_QQ_Vtx, 32000, 0);
    fTree->Branch("Reco_QQ_VxxE",     Reco_QQ_VxxE,    "Reco_QQ_VxxE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VyyE",     Reco_QQ_VyyE,    "Reco_QQ_VyyE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VzzE",     Reco_QQ_VzzE,    "Reco_QQ_VzzE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VyxE",     Reco_QQ_VyxE,    "Reco_QQ_VyxE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VzxE",     Reco_QQ_VzxE,    "Reco_QQ_VzxE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_VzyE",     Reco_QQ_VzyE,    "Reco_QQ_VzyE[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_lxy",      Reco_QQ_lxy,     "Reco_QQ_lxy[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_lxyErr",   Reco_QQ_lxyErr,  "Reco_QQ_lxyErr[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_normChi2", Reco_QQ_normChi2,"Reco_QQ_normChi2[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_probChi2", Reco_QQ_probChi2,"Reco_QQ_probChi2[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha,"Reco_QQ_cosAlpha[Reco_QQ_size]/D");
    fTree->Branch("Reco_QQ_ctau",     Reco_QQ_ctau,    "Reco_QQ_ctau[Reco_QQ_size]/D");
    if(theStoreWSOnia){
      fTree->Branch("Reco_QQ_sign",     Reco_QQ_sign,    "Reco_QQ_ctau[Reco_QQ_size]/I");
    }
    if(theStoreChicFlag){
      fTree->Branch("Reco_Chic_size",    &Reco_Chic_size,    "Reco_Chic_size/I");
      fTree->Branch("Reco_Chic_4mom",    "TClonesArray",       &Reco_Chic_4mom, 32000, 0);
      fTree->Branch("Reco_Chic_OniaDaug", Reco_Chic_OniaDaug,    "Reco_Chic_OniaDaug[Reco_Chic_size]/I");
      fTree->Branch("Reco_Chic_GammaDaug",Reco_Chic_GammaDaug,   "Reco_Chic_GammaDaug[Reco_Chic_size]/I");
      fTree->Branch("Reco_Chic_DeltaM",   Reco_Chic_DeltaM,  "Reco_Chic_DeltaM[Reco_Chic_size]/D");
    }
    if(theStoreBpFlag){
      fTree->Branch("Reco_Bp_size",    &Reco_Bp_size,    "Reco_Bp_size/I");
      fTree->Branch("Reco_Bp_4mom",    "TClonesArray",       &Reco_Bp_4mom, 32000, 0);
      fTree->Branch("Reco_Bp_OniaDaug", Reco_Bp_OniaDaug, "Reco_Bp_OniaDaug[Reco_Bp_size]/I");
      fTree->Branch("Reco_Bp_KDaug",    Reco_Bp_KDaug,    "Reco_Bp_KDaug[Reco_Bp_size]/I"); 
      fTree->Branch("Reco_Bp_VtxIsVal", Reco_Bp_VtxIsVal,"Reco_Bp_VtxIsVal[Reco_Bp_size]/B");
      fTree->Branch("Reco_Bp_Vtx",     "TClonesArray",       &Reco_Bp_Vtx, 32000, 0);
      fTree->Branch("Reco_Bp_VxE",     Reco_Bp_VxE,    "Reco_Bp_VxE[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_VyE",     Reco_Bp_VyE,    "Reco_Bp_VyE[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_VzE",     Reco_Bp_VzE,    "Reco_Bp_VzE[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_lxy",      Reco_Bp_lxy,     "Reco_Bp_lxy[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_lxyErr",   Reco_Bp_lxyErr,  "Reco_Bp_lxyErr[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_normChi2", Reco_Bp_normChi2,"Reco_Bp_normChi2[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_cosAlpha", Reco_Bp_cosAlpha,"Reco_Bp_cosAlpha[Reco_Bp_size]/D");
      fTree->Branch("Reco_Bp_ctau",     Reco_Bp_ctau,    "Reco_Bp_ctau[Reco_Bp_size]/D");
    }
  }
  
  if(theStoreOniaRadiation){
    fTree->Branch("Mc_chargedtrk_size",    &Mc_chargedtrk_size, "Mc_chargedtrk_size/I");
    fTree->Branch("Mc_chargedtrk_4mom",    "TClonesArray",       &Mc_chargedtrk_4mom, 32000, 0);
    fTree->Branch("Mc_chargedtrk_charge",   Mc_chargedtrk_charge,"Mc_chargedtrk_charge[Mc_chargedtrk_size]/I");
    fTree->Branch("Reco_PriVtx_1st_trkSize",  &Reco_PriVtx_1st_trkSize, "Reco_PriVtx_1st_trkSize/I");
    fTree->Branch("Reco_PriVtx_1st_trkindex", Reco_PriVtx_1st_trkindex, "Reco_PriVtx_1st_trkindex[Reco_PriVtx_1st_trkSize]/I");
  }  

  if(theStoreBeamSpotFlag){  
    fTree->Branch("Reco_BeamSpot_x",     &Reco_BeamSpot_x,  "Reco_BeamSpot_x/D");
    fTree->Branch("Reco_BeamSpot_y",     &Reco_BeamSpot_y,  "Reco_BeamSpot_y/D");
    fTree->Branch("Reco_BeamSpot_z",     &Reco_BeamSpot_z,  "Reco_BeamSpot_z/D");
    if(!theminimumFlag) {
      fTree->Branch("Reco_BeamSpot_xxE",    &Reco_BeamSpot_xxE, "Reco_BeamSpot_xxE/D");
      fTree->Branch("Reco_BeamSpot_yyE",    &Reco_BeamSpot_yyE, "Reco_BeamSpot_yyE/D");
      fTree->Branch("Reco_BeamSpot_zzE",    &Reco_BeamSpot_zzE, "Reco_BeamSpot_zzE/D");
      fTree->Branch("Reco_BeamSpot_yxE",    &Reco_BeamSpot_yxE, "Reco_BeamSpot_yxE/D");
      fTree->Branch("Reco_BeamSpot_zyE",    &Reco_BeamSpot_zyE, "Reco_BeamSpot_zyE/D");
      fTree->Branch("Reco_BeamSpot_zxE",    &Reco_BeamSpot_zxE, "Reco_BeamSpot_zxE/D");
    }
  }

  if(theStorePriVtxFlag){
    fTree->Branch("Reco_PriVtx_size",    &Reco_PriVtx_size,    "Reco_PriVtx_size/I"); 
    fTree->Branch("Reco_PriVtx_3vec",    "TClonesArray",       &Reco_PriVtx_3vec, 32000, 0); 
    if(!theminimumFlag) {
      fTree->Branch("Reco_PriVtx_xxE",       Reco_PriVtx_xxE,      "Reco_PriVtx_xxE[Reco_PriVtx_size]/D"); 
      fTree->Branch("Reco_PriVtx_yyE",       Reco_PriVtx_yyE,      "Reco_PriVtx_yyE[Reco_PriVtx_size]/D"); 
      fTree->Branch("Reco_PriVtx_zzE",       Reco_PriVtx_zzE,      "Reco_PriVtx_zzE[Reco_PriVtx_size]/D"); 
      fTree->Branch("Reco_PriVtx_yxE",       Reco_PriVtx_yxE,      "Reco_PriVtx_yxE[Reco_PriVtx_size]/D"); 
      fTree->Branch("Reco_PriVtx_zyE",       Reco_PriVtx_zyE,      "Reco_PriVtx_zyE[Reco_PriVtx_size]/D"); 
      fTree->Branch("Reco_PriVtx_zxE",       Reco_PriVtx_zxE,      "Reco_PriVtx_zxE[Reco_PriVtx_size]/D"); 
    }
    fTree->Branch("Reco_PriVtx_trkSize" , Reco_PriVtx_trkSize, "Reco_PriVtx_trkSize[Reco_PriVtx_size]/I");
    fTree->Branch("Reco_PriVtx_chi2",     Reco_PriVtx_chi2,    "Reco_PriVtx_chi2[Reco_PriVtx_size]/D"); 
    fTree->Branch("Reco_PriVtx_ndof",     Reco_PriVtx_ndof,    "Reco_PriVtx_ndof[Reco_PriVtx_size]/D"); 
  }

  if(theStoreL1Flag){
    fTree->Branch("L1TBits_size",        &L1TBits_size,        "L1TBits_size/I");
    fTree->Branch("L1TBits_accept",       L1TBits_accept,      "L1TBits_accept[L1TBits_size]/B");
    fTree->Branch("L1TGlobal_Decision",  &L1TGlobal_Decision,  "L1TGlobal_Decision/B");

    fTree->Branch("L1_mu_size",    &L1_mu_size,    "L1_mu_size/I");
    fTree->Branch("L1_mu_4mom",    "TClonesArray",       &L1_mu_4mom, 32000, 0);
    fTree->Branch("L1_mu_charge",   L1_mu_charge,  "L1_mu_charge[L1_mu_size]/I");

  }
  
  if(theStoreHLTFlag){
    fTree->Branch("HLTBits_size",        &HLTBits_size,        "HLTBits_size/I");
    fTree->Branch("HLTBits_wasrun",       HLTBits_wasrun,      "HLTBits_wasrun[HLTBits_size]/B");
    fTree->Branch("HLTBits_accept",       HLTBits_accept,      "HLTBits_accept[HLTBits_size]/B");
    fTree->Branch("HLTBits_error",        HLTBits_error,       "HLTBits_error[HLTBits_size]/B");
    fTree->Branch("HLTGlobal_wasrun",    &HLTGlobal_wasrun,    "HLTGlobal_wasrun/B");
    fTree->Branch("HLTGlobal_Decision",  &HLTGlobal_Decision,  "HLTGlobal_Decision/B");
    fTree->Branch("HLTGlobal_error",     &HLTGlobal_error,     "HLTGlobal_error/B");

    /* fTree->Branch("HLT1Mu3_L2_size",  &HLT1Mu3_L2_size,  "HLT1Mu3_L2_size/I");
    fTree->Branch("HLT1Mu3_L2_4mom",  "TClonesArray",    &HLT1Mu3_L2_4mom, 32000, 0);
    fTree->Branch("HLT1Mu3_L2_id",    HLT1Mu3_L2_id,     "HLT1Mu3_L2_id[HLT1Mu3_L2_size]/I");*/
    fTree->Branch("HLT1Mu3_L3_size",  &HLT1Mu3_L3_size,  "HLT1Mu3_L3_size/I");
    fTree->Branch("HLT1Mu3_L3_4mom",  "TClonesArray",    &HLT1Mu3_L3_4mom, 32000, 0);
    fTree->Branch("HLT1Mu3_L3_id",    HLT1Mu3_L3_id,     "HLT1Mu3_L3_id[HLT1Mu3_L3_size]/I");

    /* fTree->Branch("HLT1Mu5_L2_size",  &HLT1Mu5_L2_size,  "HLT1Mu5_L2_size/I");
    fTree->Branch("HLT1Mu5_L2_4mom",  "TClonesArray",    &HLT1Mu5_L2_4mom, 32000, 0);
    fTree->Branch("HLT1Mu5_L2_id",    HLT1Mu5_L2_id,     "HLT1Mu5_L2_id[HLT1Mu5_L2_size]/I"); */
    fTree->Branch("HLT1Mu5_L3_size",  &HLT1Mu5_L3_size,  "HLT1Mu5_L3_size/I");
    fTree->Branch("HLT1Mu5_L3_4mom",  "TClonesArray",    &HLT1Mu5_L3_4mom, 32000, 0);
    fTree->Branch("HLT1Mu5_L3_id",    HLT1Mu5_L3_id,     "HLT1Mu5_L3_id[HLT1Mu5_L3_size]/I");

    /* fTree->Branch("HLT1Mu9_L2_size",  &HLT1Mu9_L2_size,  "HLT1Mu9_L2_size/I");
    fTree->Branch("HLT1Mu9_L2_4mom",  "TClonesArray",    &HLT1Mu9_L2_4mom, 32000, 0);
    fTree->Branch("HLT1Mu9_L2_id",    HLT1Mu9_L2_id,     "HLT1Mu9_L2_id[HLT1Mu9_L2_size]/I"); */
    fTree->Branch("HLT1Mu9_L3_size",  &HLT1Mu9_L3_size,  "HLT1Mu9_L3_size/I");
    fTree->Branch("HLT1Mu9_L3_4mom",  "TClonesArray",    &HLT1Mu9_L3_4mom, 32000, 0);
    fTree->Branch("HLT1Mu9_L3_id",    HLT1Mu9_L3_id,     "HLT1Mu9_L3_id[HLT1Mu9_L3_size]/I");

    /* fTree->Branch("HLT1Mu11_L2_size",  &HLT1Mu11_L2_size,  "HLT1Mu11_L2_size/I");
    fTree->Branch("HLT1Mu11_L2_4mom",  "TClonesArray",     &HLT1Mu11_L2_4mom, 32000, 0);
    fTree->Branch("HLT1Mu11_L2_id",    HLT1Mu11_L2_id,     "HLT1Mu11_L2_id[HLT1Mu11_L2_size]/I"); */
    fTree->Branch("HLT1Mu11_L3_size",  &HLT1Mu11_L3_size,  "HLT1Mu11_L3_size/I");
    fTree->Branch("HLT1Mu11_L3_4mom",  "TClonesArray",     &HLT1Mu11_L3_4mom, 32000, 0);
    fTree->Branch("HLT1Mu11_L3_id",    HLT1Mu11_L3_id,     "HLT1Mu11_L3_id[HLT1Mu11_L3_size]/I");

    /* fTree->Branch("HLT2Mu0_L2_size",  &HLT2Mu0_L2_size,  "HLT2Mu0_L2_size/I");
    fTree->Branch("HLT2Mu0_L2_4mom",  "TClonesArray",       &HLT2Mu0_L2_4mom, 32000, 0);
    fTree->Branch("HLT2Mu0_L2_id",    HLT2Mu0_L2_id,     "HLT2Mu0_L2_id[HLT2Mu0_L2_size]/I"); */
    fTree->Branch("HLT2Mu0_L3_size",  &HLT2Mu0_L3_size,  "HLT2Mu0_L3_size/I");
    fTree->Branch("HLT2Mu0_L3_4mom",  "TClonesArray",       &HLT2Mu0_L3_4mom, 32000, 0);
    fTree->Branch("HLT2Mu0_L3_id",    HLT2Mu0_L3_id,     "HLT2Mu0_L3_id[HLT2Mu0_L3_size]/I");

    /* fTree->Branch("HLT2IsoMu3_L2_size",  &HLT2IsoMu3_L2_size,  "HLT2IsoMu3_L2_size/I");
    fTree->Branch("HLT2IsoMu3_L2_4mom",  "TClonesArray",     &HLT2IsoMu3_L2_4mom, 32000, 0);
    fTree->Branch("HLT2IsoMu3_L2_id",    HLT2IsoMu3_L2_id,     "HLT2IsoMu3_L2_id[HLT2IsoMu3_L2_size]/I"); */
    fTree->Branch("HLT2IsoMu3_L3_size",  &HLT2IsoMu3_L3_size,  "HLT2IsoMu3_L3_size/I");
    fTree->Branch("HLT2IsoMu3_L3_4mom",  "TClonesArray",     &HLT2IsoMu3_L3_4mom, 32000, 0);
    fTree->Branch("HLT2IsoMu3_L3_id",    HLT2IsoMu3_L3_id,     "HLT2IsoMu3_L3_id[HLT2IsoMu3_L3_size]/I");

    /* fTree->Branch("HLT2Mu3_L2_size",  &HLT2Mu3_L2_size,  "HLT2Mu3_L2_size/I");
    fTree->Branch("HLT2Mu3_L2_4mom",  "TClonesArray",       &HLT2Mu3_L2_4mom, 32000, 0);
    fTree->Branch("HLT2Mu3_L2_id",    HLT2Mu3_L2_id,     "HLT2Mu3_L2_id[HLT2Mu3_L2_size]/I"); */
    fTree->Branch("HLT2Mu3_L3_size",  &HLT2Mu3_L3_size,  "HLT2Mu3_L3_size/I");
    fTree->Branch("HLT2Mu3_L3_4mom",  "TClonesArray",       &HLT2Mu3_L3_4mom, 32000, 0);
    fTree->Branch("HLT2Mu3_L3_id",    HLT2Mu3_L3_id,     "HLT2Mu3_L3_id[HLT2Mu3_L3_size]/I");


    /* fTree->Branch("HLTJpsi2Mu_L2_size",  &HLTJpsi2Mu_L2_size,  "HLTJpsi2Mu_L2_size/I");
    fTree->Branch("HLTJpsi2Mu_L2_4mom",  "TClonesArray",       &HLTJpsi2Mu_L2_4mom, 32000, 0);
    fTree->Branch("HLTJpsi2Mu_L2_id",    HLTJpsi2Mu_L2_id,     "HLTJpsi2Mu_L2_id[HLTJpsi2Mu_L2_size]/I"); */
    fTree->Branch("HLTJpsi2Mu_L3_size",  &HLTJpsi2Mu_L3_size,  "HLTJpsi2Mu_L3_size/I");
    fTree->Branch("HLTJpsi2Mu_L3_4mom",  "TClonesArray",       &HLTJpsi2Mu_L3_4mom, 32000, 0);
    fTree->Branch("HLTJpsi2Mu_L3_id",    HLTJpsi2Mu_L3_id,     "HLTJpsi2Mu_L3_id[HLTJpsi2Mu_L3_size]/I");

    /* fTree->Branch("HLTUpsilon2Mu_L2_size",  &HLTUpsilon2Mu_L2_size,  "HLTUpsilon2Mu_L2_size/I");
    fTree->Branch("HLTUpsilon2Mu_L2_4mom",  "TClonesArray",       &HLTUpsilon2Mu_L2_4mom, 32000, 0);
    fTree->Branch("HLTUpsilon2Mu_L2_id",    HLTUpsilon2Mu_L2_id,     "HLTUpsilon2Mu_L2_id[HLTUpsilon2Mu_L2_size]/I"); */
    fTree->Branch("HLTUpsilon2Mu_L3_size",  &HLTUpsilon2Mu_L3_size,  "HLTUpsilon2Mu_L3_size/I");
    fTree->Branch("HLTUpsilon2Mu_L3_4mom",  "TClonesArray",       &HLTUpsilon2Mu_L3_4mom, 32000, 0);
    fTree->Branch("HLTUpsilon2Mu_L3_id",    HLTUpsilon2Mu_L3_id,     "HLTUpsilon2Mu_L3_id[HLTUpsilon2Mu_L3_size]/I");
 
    HLTBits_size = NTRIGGERS;

    // Level-2 FILTERS (module names)
    hltModules[0][0] = edm::InputTag("NotUsed","",the8e29ProcName);
    hltModules[0][1] = edm::InputTag("NotUsed","",the8e29ProcName);
    hltModules[0][2] = edm::InputTag("NotUsed","",the8e29ProcName);
    hltModules[0][3] = edm::InputTag("NotUsed","",the8e29ProcName);
    hltModules[0][4] = edm::InputTag("NotUsed","",the8e29ProcName);
    // Level-3 FILTERS (module names)
    hltModules[1][0] = edm::InputTag("hltSingleMu3L3Filtered3","",the8e29ProcName);
    hltModules[1][1] = edm::InputTag("hltSingleMu5L3Filtered5","",the8e29ProcName);
    hltModules[1][2] = edm::InputTag("hltSingleMu9L3Filtered9","",the8e29ProcName);
    hltModules[1][3] = edm::InputTag("hltDiMuonL3PreFiltered0","",the8e29ProcName);
    hltModules[1][4] = edm::InputTag("hltDiMuonL3PreFiltered","",the8e29ProcName);
  }
  
}

void Onia2MuMu::beginRun(const Run & iRun, const EventSetup & iSetup)
{
  if(theStoreHLTFlag){
    bool isChanged;
    if(hltConfig.init(iRun, iSetup, the8e29ProcName, isChanged)){
      if(isChanged){
        string HLTbitNames[NTRIGGERS] = {"HLT_Mu3", "HLT_Mu5", "HLT_Mu9", "HLT_DoubleMu0", "HLT_DoubleMu3"};
        // check if trigger name in config
        const unsigned int n(hltConfig.size());
        for (int ihlt = 0; ihlt < NTRIGGERS; ihlt++) {
          hltBits[ihlt] = 0;
          unsigned int triggerIndex( hltConfig.triggerIndex(HLTbitNames[ihlt]) );
          if (triggerIndex>=n) {
            cout << "OniaToMuMu::beginRun: "
                 << " TriggerName " << HLTbitNames[ihlt]
                 << " not available in config!" << endl;
          } else {
            hltBits[ihlt] = triggerIndex;
          }
        }
      }      
    }else{
      cout << "OniaToMuMu::beginRun:"
           << " HLT config extraction failure with process name HLT" << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////
// called for each event
////////////////////////////////////////////////////////////////////////
void Onia2MuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fNevt++;
  // get TTRHBuilder 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // BUILD USEFUL COLLECTIONS
  // -- 1 - get the collection of RecoTracks 
  iEvent.getByLabel(theTrackLabel, allTracks);  

  // -- 2 - get the collection of PF Photons
  iEvent.getByLabel(thePhotonLabel, pfAll);
  pfClusters.clear();
  for(PFCandidateCollection::const_iterator itPFC = pfAll->begin();
      itPFC != pfAll->end(); ++itPFC) {
    if (PFCandidate::ParticleType(itPFC->particleId()) == PFCandidate::gamma && itPFC->energy() > thePhotonMinE) {
      pfClusters.push_back(*itPFC);
    }
  }

  // -- 3 - Get NON-OVERLAPPING collections of muons
  Handle<reco::MuonCollection> allmuons;
  iEvent.getByLabel(theMuonsLabel,allmuons);

  edm::Handle<reco::CaloMuonCollection> allcalmuons;
  iEvent.getByLabel(theCaloMuonsLabel, allcalmuons);

  std::vector<int> theMuonTrkIndexes;

  for (MuonCollection::const_iterator nmuon = allmuons->begin();
       nmuon != allmuons->end(); ++nmuon) {
    if (nmuon->isGlobalMuon()) {
      theGlobalMuons.push_back(*nmuon);
      theMuonTrkIndexes.push_back(nmuon->innerTrack().index());
    }
    if (!(nmuon->isGlobalMuon()) && nmuon->isTrackerMuon()) {
      theTrkMuons.push_back(*nmuon);
      theMuonTrkIndexes.push_back(nmuon->innerTrack().index());
    }
  } 

  for (CaloMuonCollection::const_iterator cmuon = allcalmuons->begin();
       cmuon != allcalmuons->end(); ++cmuon) {
    bool storeThisOne = true; 
    for (int i = 0; i<(int)theMuonTrkIndexes.size(); ++i) {
      if ((int)cmuon->track().index() == theMuonTrkIndexes.at(i)) storeThisOne = false;
    }
    if (storeThisOne) theCaloMuons.push_back(*cmuon); 
  }   

  // -- 4 - Get track collection of NON-muons
  int k = 0;
  for (TrackCollection::const_iterator thetrk = allTracks->begin();
       thetrk != allTracks->end(); ++thetrk) {
    k++;
    bool storeThisOther = true;
    for (int j = 0; j<(int)theMuonTrkIndexes.size(); ++j) {
      if (k == theMuonTrkIndexes.at(j)) storeThisOther = false;
    }
    if (storeThisOther) noMuonTracks.push_back(*thetrk);
  } 

  if(theStoreL1Flag)       l1Report(iEvent);
  if(theStoreHLTFlag)      hltReport(iEvent);
  if(theStoreGenFlag)      fillGeneratorBlock(iEvent);
  if(theStoreTrkFlag)      fillRecTracks(iEvent);
  if(theStorePhotFlag)     fillPhotons(iEvent);
  if(theStoreGLBMuonFlag || theStoreTRKMuonFlag || theStoreCALMuonFlag)  fillMuons(iEvent);
  if(theStorePATFlag)      fillPATMuons(iEvent);
  if(theStorePATFlag)      PAT_l1Report(iEvent);
  if(theStorePATFlag)      PAT_hltReport(iEvent);
  if(theStoreBeamSpotFlag) fillBeamSpot(iEvent); 
  if(theStorePriVtxFlag&&!theUsePrimaryNoMuons)   fillPrimaryVertex(iEvent);
  if(theStorePriVtxFlag&&theUsePrimaryNoMuons)    fillPrimaryVertex2m(iEvent);
  if(theStoreOniaFlag)     findOniaCategories(iEvent);
   
  eventNb = iEvent.id().event();
  runNb = iEvent.id().run(); 
  lumiBlock = iEvent.luminosityBlock();
    
  if(!theSkimOnOniaMaxCat || Reco_QQ_size > 0) fTree->Fill(); 

  // CLEAR ALL VECTORS
  pfClusters.clear();
  noMuonTracks.clear();
  theGlobalMuons.clear();
  theTrkMuons.clear();
  theCaloMuons.clear();
  Mc_QQ_4mom->Clear();
  Mc_QQ_3vec->Clear();
  Mc_QQmoth_4mom->Clear();
  Mc_QQmoth_3vec->Clear();
  Mc_mu_4mom->Clear();
  Mc_mu_3vec->Clear();
  Reco_mu_glb_4mom->Clear();
  Reco_mu_glb_track4mom->Clear();
  Reco_mu_glb_3vec->Clear();
  Reco_mu_trk_4mom->Clear();
  Reco_mu_trk_3vec->Clear();
  Reco_mu_cal_4mom->Clear();
  Reco_mu_cal_3vec->Clear();
  Pat_mu_glb_4mom->Clear();
  Pat_mu_glb_3vec->Clear();
  Pat_mu_sta_4mom->Clear();
  Pat_mu_sta_3vec->Clear();
  Pat_mu_trk_4mom->Clear();
  Pat_mu_trk_3vec->Clear();
  Pat_mu_cal_4mom->Clear();
  Pat_mu_cal_3vec->Clear();
  // Reco_mu_sta_4mom->Clear();
  // Reco_mu_sta_3vec->Clear();
  Reco_QQ_4mom->Clear();
  Reco_QQ_Vtx->Clear();
  Reco_Chic_4mom->Clear();
  Reco_Bp_4mom->Clear();
  Reco_Bp_Vtx->Clear();
  Reco_PriVtx_3vec->Clear();
  Reco_track_4mom->Clear();
  Reco_track_3vec->Clear();
  Reco_gamma_4mom->Clear();
  Reco_track_CovM->Clear();
  // Reco_mu_glb_CovM->Clear();
  // Reco_mu_sta_CovM->Clear();
  L1_mu_4mom->Clear();
  // HLT1Mu3_L2_4mom->Clear();
  HLT1Mu3_L3_4mom->Clear();
  // HLT1Mu5_L2_4mom->Clear();
  HLT1Mu5_L3_4mom->Clear();
  // HLT1Mu9_L2_4mom->Clear();
  HLT1Mu9_L3_4mom->Clear();
  // HLT1Mu11_L2_4mom->Clear();
  HLT1Mu11_L3_4mom->Clear();
  // HLT2Mu0_L2_4mom->Clear();
  HLT2Mu0_L3_4mom->Clear();
  // HLT2IsoMu3_L2_4mom->Clear();
  HLT2IsoMu3_L3_4mom->Clear();
  // HLT2Mu3_L2_4mom->Clear();
  HLT2Mu3_L3_4mom->Clear();
  // HLTJpsi2Mu_L2_4mom->Clear();
  HLTJpsi2Mu_L3_4mom->Clear();
  HLTUpsilon2Mu_L3_4mom->Clear();  

}


///////////////////////////////////////////////////////////////
// called at end
///////////////////////////////////////////////////////////////
void Onia2MuMu::endJob()
{
  outFile->cd();
  fTree->Write();
  outFile->Close();
}

///////////////////////////////////////////////////////////////
// makes l1 report and fills ntuple
///////////////////////////////////////////////////////////////
void Onia2MuMu::l1Report(const edm::Event &iEvent) {
  if(theDebugLevel>0) cout << "l1Report called" << endl;
 
  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(theL1GTReadoutRec,L1GTRR);
  if (L1GTRR.isValid()) {
    L1TGlobal_Decision = L1GTRR->decision();
    L1TBits_size=L1GTRR->decisionWord().size(); 
    for (unsigned int i=0; i!=L1GTRR->decisionWord().size()&&i<Max_trig_size; i++) {
      L1TBits_accept[i]=L1GTRR->decisionWord()[i]; 
    }
  } 

  L1_mu_size = 0;
  Handle< l1extra::L1MuonParticleCollection > L1Muons;
  // L1 extra particles missing in some releases...
  if (iEvent.getByLabel(theL1MuonLabel, L1Muons)) {
    l1extra::L1MuonParticleCollection::const_iterator l1muon;
    for( l1muon = L1Muons->begin(); l1muon != L1Muons->end() && L1_mu_size<Max_mu_size; ++ l1muon ) {
      TLorentzVector a(0.0,0.0,0.0,0.0);
      a.SetPxPyPzE(l1muon->px(),l1muon->py(),l1muon->pz(),l1muon->energy());
      new((*L1_mu_4mom)[L1_mu_size])TLorentzVector(a);
      L1_mu_charge[L1_mu_size]=l1muon->charge();
      L1_mu_size++;
    }
  }
}

///////////////////////////////////////////////////////////////
// makes hlt report and fills ntuple
///////////////////////////////////////////////////////////////
void Onia2MuMu::hltReport(const edm::Event &iEvent) {
  using namespace trigger;
  if(theDebugLevel>0) cout << "hltReport called" << endl;
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(InputTag(theHLTriggerResults,"",the8e29ProcName), HLTR);
  if (HLTR.isValid()) {
    HLTGlobal_wasrun=HLTR->wasrun();
    HLTGlobal_Decision=HLTR->accept();
    HLTGlobal_error=HLTR->error();
    
    // HLTBits_size=HLTR->size();
    for (int i=0; i<HLTBits_size && i<(int)Max_trig_size; i++) {
      HLTBits_wasrun[i]=HLTR->wasrun(hltBits[i]);
      HLTBits_accept[i]=HLTR->accept(hltBits[i]);
      HLTBits_error[i]=HLTR->error(hltBits[i]);
    }
   
    // HLT1Mu3_L2_size=0;
    HLT1Mu3_L3_size=0;
    // HLT1Mu5_L2_size=0;
    HLT1Mu5_L3_size=0;
    // HLT1Mu9_L2_size=0;
    HLT1Mu9_L3_size=0;
    // HLT1Mu11_L2_size=0;
    HLT1Mu11_L3_size=0;
    // HLT2IsoMu3_L2_size=0;
    HLT2IsoMu3_L3_size=0;

    // HLT2Mu0_L2_size=0;
    HLT2Mu0_L3_size=0;
    // HLT2Mu3_L2_size=0;
    HLT2Mu3_L3_size=0;
    // HLTJpsi2Mu_L2_size=0;
    HLTJpsi2Mu_L3_size=0;
    // HLTUpsilon2Mu_L2_size=0;
    HLTUpsilon2Mu_L3_size=0;

    Handle<TriggerEvent> trgEvent;
    bool hltF = true;
    try {
      iEvent.getByLabel(InputTag(thetriggerEventLabel,"",the8e29ProcName), trgEvent);
    }
    catch (const cms::Exception& e) {
      hltF = false;
      cout<<"Error!! No TriggerEvent with label " << thetriggerEventLabel << endl;
    }
    if ( hltF ) {
      const TriggerObjectCollection& TOC(trgEvent->getObjects());
 
      for ( int lvl = 1; lvl<2; lvl++ ) { 
        for ( int ipath = 0; ipath < HLTBits_size; ipath++) {
          int muonsize = 0;
          const InputTag trigName = hltModules[lvl][ipath];
          size_type index = trgEvent->filterIndex(trigName);
          if ( index < trgEvent->sizeFilters() ) {
            const Keys& KEYS( trgEvent->filterKeys(index) );
            muonsize = KEYS.size();
            if(theDebugLevel>0) cout<<"Muon passing filters "<<trigName.label()<< " = " <<muonsize<<endl;
            int minNMuons = 1; if ( ipath>=3 ) minNMuons = 2;
            if (  muonsize < minNMuons && HLTR->accept(hltBits[ipath]) ) {  
              cout<<"Error!! Not enough HLT muons for "<<trigName.label()<<", but decision = "<<HLTR->accept(hltBits[ipath])<<endl;
            }
            for ( int hltm = 0; hltm < muonsize; hltm++ ) {
              size_type hltf = KEYS[hltm];
              const TriggerObject& TO(TOC[hltf]);
              TLorentzVector a = lorentzTriObj(TO);

              /* if ( lvl==0 ) {
                if ( ipath==0 ) {
                  new((*HLT1Mu3_L2_4mom)[HLT1Mu3_L2_size])TLorentzVector(a);
                  HLT1Mu3_L2_id[HLT1Mu3_L2_size]=TO.id();
		  HLT1Mu3_L2_size++;          
                }
                if ( ipath==1 ) {
                  new((*HLT1Mu5_L2_4mom)[HLT1Mu5_L2_size])TLorentzVector(a);
                  HLT1Mu5_L2_id[HLT1Mu5_L2_size]=TO.id();
		  HLT1Mu5_L2_size++;      
                }
                if ( ipath==2 ) {
		  new((*HLT1Mu9_L2_4mom)[HLT1Mu9_L2_size])TLorentzVector(a);
                  HLT1Mu9_L2_id[HLT1Mu9_L2_size]=TO.id();
		  HLT1Mu9_L2_size++;      
                }
                if ( ipath==3 ) {
		  new((*HLT1Mu11_L2_4mom)[HLT1Mu11_L2_size])TLorentzVector(a);
                  HLT1Mu11_L2_id[HLT1Mu11_L2_size]=TO.id();
		  HLT1Mu11_L2_size++;      
                }
                if ( ipath==4 ) {
                  new((*HLT2IsoMu3_L2_4mom)[HLT2IsoMu3_L2_size])TLorentzVector(a);
                  HLT2IsoMu3_L2_id[HLT2IsoMu3_L2_size]=TO.id();
		  HLT2IsoMu3_L2_size++;
                }
                if ( ipath==5 ) {
                  new((*HLT2Mu3_L2_4mom)[HLT2Mu3_L2_size])TLorentzVector(a);
                  HLT2Mu3_L2_id[HLT2Mu3_L2_size]=TO.id();
                  HLT2Mu3_L2_size++;
                }
                if ( ipath==6 ) {
                  new((*HLTJpsi2Mu_L2_4mom)[HLTJpsi2Mu_L2_size])TLorentzVector(a);
                  HLTJpsi2Mu_L2_id[HLTJpsi2Mu_L2_size]=TO.id();
                  HLTJpsi2Mu_L2_size++;
                }
                if ( ipath==7 ) {
                  new((*HLTUpsilon2Mu_L2_4mom)[HLTUpsilon2Mu_L2_size])TLorentzVector(a);
                  HLTUpsilon2Mu_L2_id[HLTUpsilon2Mu_L2_size]=TO.id();
                  HLTUpsilon2Mu_L2_size++;
                }
		} */

              if ( lvl==1 ) {
                if ( ipath==0 ) {
                  new((*HLT1Mu3_L3_4mom)[HLT1Mu3_L3_size])TLorentzVector(a);
                  HLT1Mu3_L3_id[HLT1Mu3_L3_size]=TO.id();
		  HLT1Mu3_L3_size++;          
                }
                if ( ipath==1 ) {
                  new((*HLT1Mu5_L3_4mom)[HLT1Mu5_L3_size])TLorentzVector(a);
                  HLT1Mu5_L3_id[HLT1Mu5_L3_size]=TO.id();
		  HLT1Mu5_L3_size++;      
                }
                if ( ipath==2 ) {
		  new((*HLT1Mu9_L3_4mom)[HLT1Mu9_L3_size])TLorentzVector(a);
                  HLT1Mu9_L3_id[HLT1Mu9_L3_size]=TO.id();
		  HLT1Mu9_L3_size++;      
                }
                /* if ( ipath==3 ) {
		  new((*HLT1Mu11_L3_4mom)[HLT1Mu11_L3_size])TLorentzVector(a);
                  HLT1Mu11_L3_id[HLT1Mu11_L3_size]=TO.id();
		  HLT1Mu11_L3_size++;      
                }
                if ( ipath==4 ) {
                  new((*HLT2IsoMu3_L3_4mom)[HLT2IsoMu3_L3_size])TLorentzVector(a);
                  HLT2IsoMu3_L3_id[HLT2IsoMu3_L3_size]=TO.id();
		  HLT2IsoMu3_L3_size++;
		  } */
                if ( ipath==3 ) {
                  new((*HLT2Mu0_L3_4mom)[HLT2Mu0_L3_size])TLorentzVector(a);
                  HLT2Mu0_L3_id[HLT2Mu0_L3_size]=TO.id();
                  HLT2Mu0_L3_size++;
                }
                if ( ipath==4 ) {
                  new((*HLT2Mu3_L3_4mom)[HLT2Mu3_L3_size])TLorentzVector(a);
                  HLT2Mu3_L3_id[HLT2Mu3_L3_size]=TO.id();
                  HLT2Mu3_L3_size++;
                }
                /* if ( ipath==6 ) {
                  new((*HLTJpsi2Mu_L3_4mom)[HLTJpsi2Mu_L3_size])TLorentzVector(a);
                  HLTJpsi2Mu_L3_id[HLTJpsi2Mu_L3_size]=TO.id();
                  HLTJpsi2Mu_L3_size++;
                }
                if ( ipath==7 ) {
                  new((*HLTUpsilon2Mu_L3_4mom)[HLTUpsilon2Mu_L3_size])TLorentzVector(a);
                  HLTUpsilon2Mu_L3_id[HLTUpsilon2Mu_L3_size]=TO.id();
                  HLTUpsilon2Mu_L3_size++;
		  }*/
              }
            } 
          } else {
	    if(theDebugLevel>0) cout<<"Filter index not found for filter name "<<trigName<<endl;
	  }
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////
// makes l1 report and fills ntuple
///////////////////////////////////////////////////////////////
void Onia2MuMu::PAT_l1Report(const edm::Event &iEvent) {
  cout << "L1 dummy function " << endl;
}

///////////////////////////////////////////////////////////////
// makes hlt report and fills ntuple
///////////////////////////////////////////////////////////////
void Onia2MuMu::PAT_hltReport(const edm::Event &iEvent) {
  cout << "HLT dummy function " << endl;
  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(thePATMuonsLabel,muonHandle);
  edm::View<pat::Muon> muons = *muonHandle;
    
  cout << "Number of muons " << muons.size() << endl;
  for(edm::View<pat::Muon>::const_iterator muoni = muons.begin();
      muoni!=muons.end(); 
      ++muoni) {
    cout << "new muon " << endl;
    // Single Muon Trigger
    const pat::TriggerObjectStandAloneCollection matches1=muoni->triggerObjectMatchesByFilter("hltSingleMu3L3Filtered3"); 
    if(matches1.empty()){
      cout << "this mu didn't pass the single mu trigger" << endl;
    }else{
      cout << "the associated single mu trigger muon has pt " << matches1[0].pt()  << endl;
    }
    // Double Muon Trigger
    const pat::TriggerObjectStandAloneCollection matches2=muoni->triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered"); 
    if(matches2.empty()){
      cout << "this mu didn't pass the double mu trigger" << endl;
    }else{
      cout << "the associated double mu trigger muon has pt " << matches2[0].pt()  << endl;
    }
  }

}


///////////////////////////////////////////////////////////////
// fills generator block
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillGeneratorBlock(const edm::Event &iEvent) {

  if(theDebugLevel>0) cout << "==>fillGeneratorBlocks, event: " << fNevt << endl;
 
  if ( theAODFlag ) {
    /*
    Handle< int > genProcessID;
    iEvent.getByLabel( "genEventProcID", genProcessID );
    Mc_ProcessId = *genProcessID;
    */
    //    Handle< double > genEventScale;
    //iEvent.getByLabel( "genEventScale", genEventScale );
    //Mc_EventScale = *genEventScale;
    /*
    Handle< double > genFilterEff;
    iEvent.getByLabel( "genEventRunInfo", "FilterEfficiency", genFilterEff);
    double filter_eff = *genFilterEff;
    
    Handle< double > genCrossSect;
    iEvent.getByLabel( "genEventRunInfo", "PreCalculatedCrossSection", genCrossSect);
    double cross_section = *genCrossSect;
    //Mc_EventWeight =cross_section * filter_eff*branch_ratio;
    Mc_EventWeight =cross_section * filter_eff;
    */
  } else {
    Handle< HepMCProduct > HepMCEvt;
    // if ( !iEvent.getByLabel( "evtgenproducer", HepMCEvt ) ) {
    //	iEvent.getByLabel( "source", HepMCEvt );
    // }
    iEvent.getByLabel( "generator", HepMCEvt );
    const HepMC::GenEvent* myGenEvent = HepMCEvt->GetEvent();
    Mc_ProcessId   = myGenEvent->signal_process_id();
    Mc_EventScale  = myGenEvent->event_scale();
    
    Handle< GenRunInfoProduct > gi;
    // iEvent.getRun().getByLabel( "source", gi);
    iEvent.getRun().getByLabel( "generator", gi);
    double auto_cross_section = gi->internalXSec().value(); // calculated at end of each RUN, in mb
    if(theDebugLevel>0) cout << "calculated cross-section" << auto_cross_section<<endl;
    double external_cross_section = gi->crossSection(); // is the precalculated one written in the cfg file -- units is pb
    double filter_eff = gi->filterEfficiency();
    //Mc_EventWeight = external_cross_section * filter_eff*branch_ratio ;  // in pb; in analyzer weight=this weight/Nr events analyzed
    Mc_EventWeight = external_cross_section * filter_eff;
  }

  Mc_QQ_size=0; 
  Mc_mu_size=0;
  Mc_chargedtrk_size=0;
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel( thegenParticlesLabel, genParticles );


  for( size_t i = 0; i < genParticles->size(); ++ i ) {
    const Candidate & p = (*genParticles)[ i ];
    int Mc_particleID=p.pdgId();
    ////// Store muon information (all muons)
    if (abs(Mc_particleID) == 13 && p.status()==1 && Mc_mu_size<Max_mu_size){
      TLorentzVector a(0.0,0.0,0.0,0.0);
      a.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
      new((*Mc_mu_4mom)[Mc_mu_size])TLorentzVector(a);
      TVector3 b(0.0,0.0,0.0);
      b.SetXYZ(p.vertex().x(),p.vertex().y(),p.vertex().z());
      new((*Mc_mu_3vec)[Mc_mu_size])TVector3(b);
      Mc_mu_id[Mc_mu_size]=p.pdgId();
      Mc_mumoth_id[Mc_mu_size]=(p.mother())->pdgId();
      if(theDebugLevel>1) cout << "Mc_mumoth_id" << Mc_mumoth_id[Mc_mu_size] << endl;
      if(theDebugLevel>1) cout<<Mc_mu_size<<" mc muon pt="<<p.pt()<<endl;

      Mc_mu_size++;
    }
  }
 
  for( size_t i = 0; i < genParticles->size(); ++ i ) {
    const Candidate & p = (*genParticles)[ i ];
    int Mc_particleID = p.pdgId();
    if (abs(Mc_particleID) == theOniaType && p.status()==2 && Mc_QQ_size<Max_QQ_size){
      TLorentzVector a(0.0,0.0,0.0,0.0);
      a.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
      new((*Mc_QQ_4mom)[Mc_QQ_size])TLorentzVector(a);
      TVector3 b(0.0,0.0,0.0);
      b.SetXYZ(p.vertex().x(),p.vertex().y(),p.vertex().z());
      new((*Mc_QQ_3vec)[Mc_QQ_size])TVector3(b);
    
      // mother (or B-mother) information     
      const Candidate & pmo=*(p.mother());
      const Candidate & pgmo=*(pmo.mother());

      unsigned int storeWhich = 0;
      if (!isAbHadron(pmo.pdgId()) && isAbHadron(pgmo.pdgId())) storeWhich = 1;

      TLorentzVector c(0.0,0.0,0.0,0.0);
      TVector3 d(0.0,0.0,0.0);

      if (storeWhich == 0) {
	Mc_QQmoth_id[Mc_QQ_size]=pmo.pdgId();
	c.SetPxPyPzE(pmo.px(),pmo.py(),pmo.pz(),pmo.energy());
	d.SetXYZ(pmo.vertex().x(),pmo.vertex().y(),pmo.vertex().z());
      } else {
	Mc_QQmoth_id[Mc_QQ_size]=pgmo.pdgId();
	c.SetPxPyPzE(pgmo.px(),pgmo.py(),pgmo.pz(),pgmo.energy());
	d.SetXYZ(pgmo.vertex().x(),pgmo.vertex().y(),pgmo.vertex().z());
      }
	
      new((*Mc_QQmoth_4mom)[Mc_QQ_size])TLorentzVector(c);
      new((*Mc_QQmoth_3vec)[Mc_QQ_size])TVector3(d);

      /////// Now loop over children and store corresponding line in muon vector
      int nchildrenOnia = p.numberOfDaughters();
      const Candidate & da1 = *(p.daughter( 0 ));
      const Candidate & da2 = *(p.daughter( 1 ));
       // If indeed jpsi decaying into two muons
      if(nchildrenOnia==2 && abs(da1.pdgId())==13 && abs(da2.pdgId())==13){
	TLorentzVector da1_4vec(p.daughter(0)->px(),
				p.daughter(0)->py(),
				p.daughter(0)->pz(),
				p.daughter(0)->energy());
	TLorentzVector da2_4vec(p.daughter(1)->px(),
				p.daughter(1)->py(),
				p.daughter(1)->pz(),
				p.daughter(1)->energy());
	for( int j=0; j<Mc_mu_size; j++){
	  if((p.daughter(0)->pdgId())== Mc_mu_id[j] && fabs(da1_4vec.Pt() - ((TLorentzVector*)Mc_mu_4mom->At(j))->Pt()) < 0.001){
	    if(Mc_mu_id[j]==13) Mc_QQmumi_indx[Mc_QQ_size]=j;
	    if(Mc_mu_id[j]==-13)Mc_QQmupl_indx[Mc_QQ_size]=j;
	  }
	  if((p.daughter(1)->pdgId())== Mc_mu_id[j] && fabs(da2_4vec.Pt() - ((TLorentzVector*)Mc_mu_4mom->At(j))->Pt()) < 0.001){
	    if(Mc_mu_id[j]==13) Mc_QQmumi_indx[Mc_QQ_size]=j;
	    if(Mc_mu_id[j]==-13)Mc_QQmupl_indx[Mc_QQ_size]=j;
	  } 
	}
      }
      
      Mc_QQ_size++;
    }
  } // end loop over genParticles

  
  if(theStoreOniaRadiation){ 
    for( size_t i = 0; i < genParticles->size(); ++ i ) {
      const Candidate & p = (*genParticles)[ i ];
      TLorentzVector a(0.0,0.0,0.0,0.0);
      a.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
      // Do not count the muons
      if(p.status()==1&&p.threeCharge()!=0&&Mc_chargedtrk_size<Max_track_size){
	TLorentzVector a;
	a.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
	new((*Mc_chargedtrk_4mom)[Mc_chargedtrk_size])TLorentzVector(a);
	Mc_chargedtrk_charge[Mc_chargedtrk_size]=p.threeCharge();
	Mc_chargedtrk_size++;
      }
    }
  } 
 
}

///////////////////////////////////////////////////////////////
// fills reconstructed tracks block
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillRecTracks(const edm::Event &iEvent) {
  
  if(theDebugLevel>0) cout << "==>fillRecTracks> Starting to fill reconstructed tracks, event: " << fNevt << endl;

  Reco_track_size=0;
  int count=0;
  for(TrackCollection::const_iterator itTrack = allTracks->begin();
      itTrack != allTracks->end()&&Reco_track_size<Max_track_size;
      ++itTrack) {
    //    count=count+1;
    if(theDebugLevel>2) cout << "track nr " <<  count << " Pt = " << itTrack->pt()<< "eta " << itTrack->eta() <<endl;
    count=count+1;
    TLorentzVector a=lorentzMomentum(*itTrack); 
    new((*Reco_track_4mom)[Reco_track_size])TLorentzVector(a); 
    TVector3 b(itTrack->vx(),itTrack->vy(),itTrack->vz());
    new((*Reco_track_3vec)[Reco_track_size])TVector3(b);
    Reco_track_ptErr[Reco_track_size]=itTrack->ptError();
    Reco_track_phiErr[Reco_track_size]=itTrack->phiError();
    Reco_track_etaErr[Reco_track_size]=itTrack->etaError();
    Reco_track_d0[Reco_track_size]=itTrack->d0();
    Reco_track_d0err[Reco_track_size]=itTrack->d0Error();
    Reco_track_dz[Reco_track_size]=itTrack->dz();
    Reco_track_dzerr[Reco_track_size]=itTrack->dzError();
    Reco_track_charge[Reco_track_size]=itTrack->charge();
    Reco_track_chi2[Reco_track_size]=itTrack->chi2();
    Reco_track_ndof[Reco_track_size]=itTrack->ndof();
    Reco_track_nhits[Reco_track_size]=itTrack->numberOfValidHits();
    
    TMatrixD cov(5,5); 
    for (int lan=0;lan<5;lan++ ) {
      for ( int len=0;len<5;len++ ) {
        cov(lan,len)=itTrack->covariance(lan,len);
      }
    }
    new((*Reco_track_CovM)[Reco_track_size])TMatrixD(cov); 
    Reco_track_size++;
  }
}

///////////////////////////////////////////////////////////////
// fills reconstructed photons
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillPhotons(const edm::Event &iEvent) {
  
  if(theDebugLevel>0) cout << "==>fillPhotons> Starting to fill reconstructed photons, event: " << fNevt << endl;

  Reco_gamma_size=0;
  int count=0;
  for(PFCandidateCollection::const_iterator itePFC = pfClusters.begin();
      itePFC != pfClusters.end()&&Reco_gamma_size<Max_track_size;
      ++itePFC) {
    PFCandidate itPFC = *itePFC;
    //    count=count+1;
    if(theDebugLevel>2) cout << "photon nr " <<  count << " E = " << itPFC.energy()<< " eta " << itPFC.eta() <<endl;
    count=count+1;
    TLorentzVector a=lorentzMomentum(itPFC); 
    new((*Reco_gamma_4mom)[Reco_gamma_size])TLorentzVector(a); 
    Reco_gamma_phi[Reco_gamma_size]=itPFC.phi();
    Reco_gamma_eta[Reco_gamma_size]=itPFC.eta();
    
    Reco_gamma_size++;
  }
}

//////////////////////////////////////////////////////////////
// fills muon block: 3 options for muons
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillMuons(const edm::Event &iEvent){
 
  if(theDebugLevel>0) cout << "fillMuons called " << endl;

  // Handle<reco::TrackCollection> stas;
  // iEvent.getByLabel(theStandAloneMuonsLabel, stas);

  // Handle<reco::TrackCollection> glbmuons;
  // iEvent.getByLabel(theGlobalMuonsLabel, glbmuons);

  /////////// StandAlone Muons
  /* if ( theStoreSTAMuonFlag ) { 
    Reco_mu_sta_size=0;    
    if(theDebugLevel>1) cout << "SIZE STA Muons " <<  stas->size() << endl;
    for (reco::TrackCollection::const_iterator muoni = stas->begin();
       muoni != stas->end()&&Reco_mu_sta_size<Max_mu_size;
       muoni++) {
      if(theDebugLevel>1) cout << "New STA muon " << endl;
      if(theDebugLevel>0) printTrack(*muoni);
      TLorentzVector a=lorentzMomentum(*muoni);
      if(theDebugLevel>1)cout << "StandAloneMuon PT " << a.Pt() << endl;
      new((*Reco_mu_sta_4mom)[Reco_mu_sta_size])TLorentzVector(a);
      TVector3 b(muoni->vx(),muoni->vy(),muoni->vz());
      new((*Reco_mu_sta_3vec)[Reco_mu_sta_size])TVector3(b);
      Reco_mu_sta_ptErr[Reco_mu_sta_size]=muoni->ptError();
      Reco_mu_sta_phiErr[Reco_mu_sta_size]=muoni->phiError();
      Reco_mu_sta_etaErr[Reco_mu_sta_size]=muoni->etaError();
      Reco_mu_sta_d0[Reco_mu_sta_size]=muoni->d0();
      Reco_mu_sta_d0err[Reco_mu_sta_size]=muoni->d0Error();
      Reco_mu_sta_dz[Reco_mu_sta_size]=muoni->dz();
      Reco_mu_sta_dzerr[Reco_mu_sta_size]=muoni->dzError();
      Reco_mu_sta_charge[Reco_mu_sta_size]=muoni->charge();
      Reco_mu_sta_chi2[Reco_mu_sta_size]=muoni->chi2();
      Reco_mu_sta_ndof[Reco_mu_sta_size]=muoni->ndof();
      Reco_mu_sta_nhits[Reco_mu_sta_size]=muoni->numberOfValidHits();

      TMatrixD cov2(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov2(lan,len)=muoni->covariance(lan,len);
        }
      }
      new((*Reco_mu_sta_CovM)[Reco_mu_sta_size])TMatrixD(cov2);

      Reco_mu_sta_size++;
    }
    }*/

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(thePrimaryVertexLabel, privtxs);
  VertexCollection::const_iterator privtx;

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel(theBeamSpotLabel ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;

  math::XYZPoint RefVtx;
  if ( theBeamSpotFlag ) {
    RefVtx = bs.position();
  }
  else {
    if ( privtxs->begin() != privtxs->end() ) {
      privtx=privtxs->begin();
      RefVtx = privtx->position();
    }
    else {
      RefVtx.SetXYZ(0, 0, 0);
    }
  }

  /////////// Global Muons 
  if ( theStoreGLBMuonFlag ) {
    Reco_mu_glb_size=0;
    if(theDebugLevel>1) cout << "SIZE Global Muons " <<  theGlobalMuons.size() << endl;
  
    for (reco::MuonCollection::const_iterator muoni = theGlobalMuons.begin();
         muoni != theGlobalMuons.end()&&Reco_mu_glb_size<Max_mu_size; 
         muoni++) {
      if(theDebugLevel>1) cout << "New GLB muon " << endl;
      TrackRef glbTrack = muoni->globalTrack();
      TrackRef innTrack = muoni->innerTrack();
      if(theDebugLevel>0) printTrack(*glbTrack);

      TLorentzVector a=lorentzMomentum(*glbTrack);
      if(theDebugLevel>1)cout << "GlobalMuon PT " << a.Pt() << endl;
      new((*Reco_mu_glb_4mom)[Reco_mu_glb_size])TLorentzVector(a); 

      TLorentzVector a2=lorentzMomentum(*innTrack);
      if(theDebugLevel>1)cout << "GlobalMuon track PT " << a2.Pt() << endl;
      new((*Reco_mu_glb_track4mom)[Reco_mu_glb_size])TLorentzVector(a2); 
 
      TVector3 b(glbTrack->vx(),glbTrack->vy(),glbTrack->vz());
      new((*Reco_mu_glb_3vec)[Reco_mu_glb_size])TVector3(b);

      Reco_mu_glb_ptErr[Reco_mu_glb_size]=glbTrack->ptError();
      Reco_mu_glb_phiErr[Reco_mu_glb_size]=glbTrack->phiError();
      Reco_mu_glb_etaErr[Reco_mu_glb_size]=glbTrack->etaError();
      Reco_mu_glb_d0[Reco_mu_glb_size]=-glbTrack->dxy(RefVtx);
      Reco_mu_glb_d0err[Reco_mu_glb_size]=glbTrack->d0Error();
      Reco_mu_glb_dz[Reco_mu_glb_size]=glbTrack->dz(RefVtx);
      Reco_mu_glb_dzerr[Reco_mu_glb_size]=glbTrack->dzError();
      Reco_mu_glb_charge[Reco_mu_glb_size]=glbTrack->charge();
      Reco_mu_glb_normChi2[Reco_mu_glb_size]=glbTrack->chi2()/glbTrack->ndof();
      // Reco_mu_glb_ndof[Reco_mu_glb_size]=glbTrack->ndof();
      Reco_mu_glb_nhitstrack[Reco_mu_glb_size]=innTrack->numberOfValidHits();
 
      // if ( !theAODFlag ) {
      std::vector<unsigned int> theHits = this->trackHits(*glbTrack);
      Reco_mu_glb_nhitsDT[Reco_mu_glb_size]=theHits.at(0);
      Reco_mu_glb_nhitsCSC[Reco_mu_glb_size]=theHits.at(1);
      Reco_mu_glb_nhitsStrip[Reco_mu_glb_size]=theHits.at(3);
      Reco_mu_glb_nhitsPixB[Reco_mu_glb_size]=theHits.at(4);
      Reco_mu_glb_nhitsPixE[Reco_mu_glb_size]=theHits.at(5);
      Reco_mu_glb_nhitsPix1Hit[Reco_mu_glb_size]=theHits.at(6);
      Reco_mu_glb_nhitsPix1HitBE[Reco_mu_glb_size]=theHits.at(7);
	// }

      Reco_mu_glb_caloComp[Reco_mu_glb_size]=muoni->caloCompatibility();
      Reco_mu_glb_segmComp[Reco_mu_glb_size]=muon::segmentCompatibility(*muoni);;
      Reco_mu_glb_iso[Reco_mu_glb_size]=muoni->isolationR03().sumPt;

      /* TMatrixD cov1(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov1(lan,len)=muoni->covariance(lan,len);
        }
      }
      new((*Reco_mu_glb_CovM)[Reco_mu_glb_size])TMatrixD(cov1); */

      Reco_mu_glb_size++;
    }
  } 

  /////////// Tracker Muons 
  if ( theStoreTRKMuonFlag ) {
    Reco_mu_trk_size=0;
    if(theDebugLevel>1) cout << "SIZE Trk Muons " <<  theTrkMuons.size() << endl;

    for (reco::MuonCollection::const_iterator trkmuoni = theTrkMuons.begin();
         trkmuoni != theTrkMuons.end()&&Reco_mu_trk_size<Max_mu_size; 
         trkmuoni++) {
      if(theDebugLevel>1) cout << "New TRK muon " << endl;
      // TrackRef glbTrack = trkmuoni->globalTrack();
      TrackRef innTrack = trkmuoni->innerTrack();
      if(theDebugLevel>0) printTrack(*innTrack);

      TLorentzVector a=lorentzMomentum(*innTrack);
      if(theDebugLevel>1)cout << "TrkMuon PT " << a.Pt() << endl;
      new((*Reco_mu_trk_4mom)[Reco_mu_trk_size])TLorentzVector(a); 

      TVector3 b(innTrack->vx(),innTrack->vy(),innTrack->vz());
      new((*Reco_mu_trk_3vec)[Reco_mu_trk_size])TVector3(b);

      Reco_mu_trk_ptErr[Reco_mu_trk_size]=innTrack->ptError();
      Reco_mu_trk_phiErr[Reco_mu_trk_size]=innTrack->phiError();
      Reco_mu_trk_etaErr[Reco_mu_trk_size]=innTrack->etaError();
      Reco_mu_trk_d0[Reco_mu_trk_size]=-innTrack->dxy(RefVtx);
      Reco_mu_trk_d0err[Reco_mu_trk_size]=innTrack->d0Error();
      Reco_mu_trk_dz[Reco_mu_trk_size]=innTrack->dz(RefVtx);
      Reco_mu_trk_dzerr[Reco_mu_trk_size]=innTrack->dzError();
      Reco_mu_trk_charge[Reco_mu_trk_size]=innTrack->charge();
      Reco_mu_trk_normChi2[Reco_mu_trk_size]=innTrack->chi2()/innTrack->ndof();
      // Reco_mu_trk_ndof[Reco_mu_trk_size]=innTrack->ndof();
      Reco_mu_trk_nhitstrack[Reco_mu_trk_size]=innTrack->numberOfValidHits();

      std::vector<unsigned int> theHits = this->trackHits(*innTrack);
      Reco_mu_trk_nhitsStrip[Reco_mu_trk_size]=theHits.at(3);
      Reco_mu_trk_nhitsPixB[Reco_mu_trk_size]=theHits.at(4);
      Reco_mu_trk_nhitsPixE[Reco_mu_trk_size]=theHits.at(5);
      Reco_mu_trk_nhitsPix1Hit[Reco_mu_trk_size]=theHits.at(6);
      Reco_mu_trk_nhitsPix1HitBE[Reco_mu_trk_size]=theHits.at(7);

      // std::vector<unsigned int> theMuonHits = this->muonStatHits(*glbTrack);
      // Reco_mu_trk_nhitsDT[Reco_mu_trk_size]=-1;
      // Reco_mu_trk_nhitsCSC[Reco_mu_trk_size]=-1;

      // STANDARD SELECTORS
      int myWord = 0;
      if (muon::isGoodMuon(*trkmuoni, muon::AllTrackerMuons))                   myWord += (int)pow(2.,0);
      if (muon::isGoodMuon(*trkmuoni, muon::TrackerMuonArbitrated))             myWord += (int)pow(2.,1);
      if (muon::isGoodMuon(*trkmuoni, muon::TMLastStationLoose))                myWord += (int)pow(2.,2);
      if (muon::isGoodMuon(*trkmuoni, muon::TMLastStationTight))                myWord += (int)pow(2.,3);
      if (muon::isGoodMuon(*trkmuoni, muon::TM2DCompatibilityLoose))            myWord += (int)pow(2.,4);
      if (muon::isGoodMuon(*trkmuoni, muon::TM2DCompatibilityTight))            myWord += (int)pow(2.,5);
      if (muon::isGoodMuon(*trkmuoni, muon::TMOneStationLoose))                 myWord += (int)pow(2.,6);
      if (muon::isGoodMuon(*trkmuoni, muon::TMOneStationTight))                 myWord += (int)pow(2.,7);
      if (muon::isGoodMuon(*trkmuoni, muon::TMLastStationOptimizedLowPtLoose))  myWord += (int)pow(2.,8);
      if (muon::isGoodMuon(*trkmuoni, muon::TMLastStationOptimizedLowPtTight))  myWord += (int)pow(2.,9);
      Reco_mu_trk_PIDmask[Reco_mu_trk_size]=myWord;

      Reco_mu_trk_caloComp[Reco_mu_trk_size]=trkmuoni->caloCompatibility();
      Reco_mu_trk_segmComp[Reco_mu_trk_size]=muon::segmentCompatibility(*trkmuoni);
      Reco_mu_trk_iso[Reco_mu_trk_size]=trkmuoni->isolationR03().sumPt;

      /* TMatrixD cov1(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov1(lan,len)=trkmuoni->covariance(lan,len);
        }
      }
      new((*Reco_mu_trk_CovM)[Reco_mu_trk_size])TMatrixD(cov1); */

      Reco_mu_trk_size++;
    }
  } 

  /////////// Calo Muons 
  if ( theStoreCALMuonFlag ) {
    Reco_mu_cal_size=0;
    if(theDebugLevel>1) cout << "SIZE Cal Muons " <<  theCaloMuons.size() << endl;
    for (reco::CaloMuonCollection::const_iterator calmuoni = theCaloMuons.begin();
         calmuoni != theCaloMuons.end()&&Reco_mu_cal_size<Max_mu_size; 
         calmuoni++) {
      if(theDebugLevel>1) cout << "New CAL muon " << endl;
  
      TrackRef innTrack = calmuoni->track();
      if(theDebugLevel>0) printTrack(*innTrack);

      TLorentzVector a=lorentzMomentum(*innTrack);
      if(theDebugLevel>1)cout << "CalMuon PT " << a.Pt() << endl;
      new((*Reco_mu_cal_4mom)[Reco_mu_cal_size])TLorentzVector(a); 

      TVector3 b(innTrack->vx(),innTrack->vy(),innTrack->vz());
      new((*Reco_mu_cal_3vec)[Reco_mu_cal_size])TVector3(b);

      Reco_mu_cal_ptErr[Reco_mu_cal_size]=innTrack->ptError();
      Reco_mu_cal_phiErr[Reco_mu_cal_size]=innTrack->phiError();
      Reco_mu_cal_etaErr[Reco_mu_cal_size]=innTrack->etaError();
      Reco_mu_cal_d0[Reco_mu_cal_size]=-innTrack->dxy(RefVtx);
      Reco_mu_cal_d0err[Reco_mu_cal_size]=innTrack->d0Error();
      Reco_mu_cal_dz[Reco_mu_cal_size]=innTrack->dz(RefVtx);
      Reco_mu_cal_dzerr[Reco_mu_cal_size]=innTrack->dzError();
      Reco_mu_cal_charge[Reco_mu_cal_size]=innTrack->charge();
      Reco_mu_cal_normChi2[Reco_mu_cal_size]=innTrack->chi2()/innTrack->ndof();
      // Reco_mu_cal_ndof[Reco_mu_cal_size]=innTrack->ndof();
      Reco_mu_cal_nhitstrack[Reco_mu_cal_size]=innTrack->numberOfValidHits();
 
      std::vector<unsigned int> theHits = this->trackHits(*innTrack);
      Reco_mu_cal_nhitsStrip[Reco_mu_cal_size]=theHits.at(3);
      Reco_mu_cal_nhitsPixB[Reco_mu_cal_size]=theHits.at(4);
      Reco_mu_cal_nhitsPixE[Reco_mu_cal_size]=theHits.at(5);
      Reco_mu_cal_nhitsPix1Hit[Reco_mu_cal_size]=theHits.at(6);
      Reco_mu_cal_nhitsPix1HitBE[Reco_mu_cal_size]=theHits.at(7);

      Reco_mu_cal_caloComp[Reco_mu_cal_size]=calmuoni->caloCompatibility();

      /* TMatrixD cov1(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov1(lan,len)=calmuoni->covariance(lan,len);
        }
      }
      new((*Reco_mu_cal_CovM)[Reco_mu_cal_size])TMatrixD(cov1); */

      Reco_mu_cal_size++;
    }
  } 

  /* Reco_mu_size=0;
  
  Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(theMuonsLabel,muons);

  if(theDebugLevel>1)cout << "SIZE muons " <<  muons->size() << endl;
  Reco_mu_Normsize=muons->size();
  for (reco::MuonCollection::const_iterator muoni=muons->begin() ;
       muoni !=muons->end()&&Reco_mu_size<Max_track_size ;
       muoni++ ){
    if(theDebugLevel>1) cout << "New MUON nr " << Reco_mu_size << endl;
    Reco_mu_is_sta[Reco_mu_size] = muoni->isStandAloneMuon();
    Reco_mu_is_glb[Reco_mu_size] = muoni->isGlobalMuon();
    Reco_mu_is_trk[Reco_mu_size] = muoni->isTrackerMuon(); 
    Reco_mu_is_cal[Reco_mu_size] = muoni->isCaloMuon();
    if(theDebugLevel>1) cout << "Reco_mu_is_sta[Reco_mu_size]" << Reco_mu_is_sta[Reco_mu_size] << endl;
    if(theDebugLevel>1) cout << "Reco_mu_is_glb[Reco_mu_size]" << Reco_mu_is_glb[Reco_mu_size] << endl;
    if(theDebugLevel>1) cout << "Reco_mu_is_trk[Reco_mu_size]" << Reco_mu_is_trk[Reco_mu_size] << endl;
    if(theDebugLevel>1) cout << "Reco_mu_is_cal[Reco_mu_size]" << Reco_mu_is_cal[Reco_mu_size] << endl;
    TrackRef glbTrack = muoni->globalTrack();
    Reco_mu_links_glb[Reco_mu_size]=glbTrack.index(); 
    if(theDebugLevel>1) cout << "Link to glb mu " << Reco_mu_links_glb[Reco_mu_size] << endl;
    TrackRef staTrack = muoni->outerTrack();
    Reco_mu_links_sta[Reco_mu_size]=staTrack.index();
    if(theDebugLevel>1) cout << "Link to sta mu " << Reco_mu_links_sta[Reco_mu_size] << endl;
    TrackRef trkTrack = muoni->innerTrack();
    Reco_mu_links_trk[Reco_mu_size]=trkTrack.index();
    if(theDebugLevel>1) cout << "Link to trk mu " << Reco_mu_links_trk[Reco_mu_size] << endl;
    Reco_mu_caloComp[Reco_mu_size]=muoni->caloCompatibility();
      if(theDebugLevel>1) cout << " calocompatibility " << muoni->caloCompatibility() << endl;
      Reco_mu_size++;
  }
  
  if(theDebugLevel>1) cout << "after normal Reco_mu_size " << Reco_mu_size << endl;
  
  edm::Handle<reco::CaloMuonCollection> calmuons;
  iEvent.getByLabel(theCaloMuonsLabel, calmuons);

  if(theDebugLevel>1)cout << "SIZE Calomuons " <<  calmuons->size() << endl;
  Reco_mu_Calmsize = calmuons->size();
  for (reco::CaloMuonCollection::const_iterator calmuoni=calmuons->begin() ;
       calmuoni !=calmuons->end()&&Reco_mu_size<Max_track_size ;
       calmuoni++ ){
    Reco_mu_is_sta[Reco_mu_size]= false;
    Reco_mu_is_glb[Reco_mu_size]= false;
    Reco_mu_is_trk[Reco_mu_size]= false;
    Reco_mu_is_cal[Reco_mu_size]= true;
    Reco_mu_links_glb[Reco_mu_size] = -10000;
    Reco_mu_links_sta[Reco_mu_size] = -10000;
    TrackRef trkTrack = calmuoni->track();
    Reco_mu_links_trk[Reco_mu_size]=trkTrack.index();
    if(theDebugLevel>1) cout << "Link to trk mu " << Reco_mu_links_trk[Reco_mu_size] << endl;    
    Reco_mu_caloComp[Reco_mu_size]=calmuoni->caloCompatibility(); 
    if(theDebugLevel>1) cout << " calocompatibility " << calmuoni->caloCompatibility() << endl;
    Reco_mu_size++;
  } 
 
  if(theDebugLevel>1) cout << "after calo Reco_mu_size " << Reco_mu_size << endl; */

}

////////////////
// PAT Muons 
///////////////
void Onia2MuMu::fillPATMuons(const edm::Event &iEvent){

  Pat_mu_glb_size=0;
  Pat_mu_sta_size=0;
  Pat_mu_trk_size=0;
  Pat_mu_cal_size=0;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(thePATMuonsLabel,muonHandle);
  edm::View<pat::Muon> muons = *muonHandle;
    
  for(edm::View<pat::Muon>::const_iterator muoni = muons.begin();
      muoni!=muons.end(); 
      ++muoni) {

    /// global muons      
    if (muoni->isGlobalMuon()) {      
      TVector3 b(muoni->vx(),muoni->vy(),muoni->vz());
      new((*Pat_mu_glb_3vec)[Pat_mu_glb_size])TVector3(b);
      Pat_mu_glb_charge[Pat_mu_glb_size]=muoni->charge();
      TLorentzVector a=lorentzMomentum(*muoni);
      if(theDebugLevel>1)cout << "PAT Muon PT " << a.Pt() << endl;
      new((*Pat_mu_glb_4mom)[Pat_mu_glb_size])TLorentzVector(a);
 
      TrackRef muonRefe = muoni->globalTrack();
      if ( muonRefe.isNull() ) { cout << "muon reference is null!"<<endl; } 
      Pat_mu_glb_ptErr[Pat_mu_glb_size]=muonRefe->ptError();
      Pat_mu_glb_phiErr[Pat_mu_glb_size]=muonRefe->phiError();
      Pat_mu_glb_etaErr[Pat_mu_glb_size]=muonRefe->etaError();
      Pat_mu_glb_d0[Pat_mu_glb_size]=muonRefe->d0();
      Pat_mu_glb_d0err[Pat_mu_glb_size]=muonRefe->d0Error();
      Pat_mu_glb_dz[Pat_mu_glb_size]=muonRefe->dz();
      Pat_mu_glb_dzerr[Pat_mu_glb_size]=muonRefe->dzError();
      Pat_mu_glb_chi2[Pat_mu_glb_size]=muonRefe->chi2();
      Pat_mu_glb_ndof[Pat_mu_glb_size]=muonRefe->ndof();
      Pat_mu_glb_nhits[Pat_mu_glb_size]=muonRefe->numberOfValidHits();
      TMatrixD cov1PAT(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov1PAT(lan,len)=muonRefe->covariance(lan,len);
        }
      }
      new((*Pat_mu_glb_CovM)[Pat_mu_glb_size])TMatrixD(cov1PAT);
      Pat_mu_glb_size++;     
    }
    /// stand alone muons      
    if (muoni->isStandAloneMuon()) {      
      TVector3 b(muoni->vx(),muoni->vy(),muoni->vz());
      new((*Pat_mu_sta_3vec)[Pat_mu_sta_size])TVector3(b);
      Pat_mu_sta_charge[Pat_mu_sta_size]=muoni->charge();
      TLorentzVector a=lorentzMomentum(*muoni);
      if(theDebugLevel>1)cout << "PAT Muon PT " << a.Pt() << endl;
      new((*Pat_mu_sta_4mom)[Pat_mu_sta_size])TLorentzVector(a); 
      TrackRef muonRefe = muoni->outerTrack();
      if ( muonRefe.isNull() ) { cout << "muon reference is null!"<<endl; } 
      Pat_mu_sta_ptErr[Pat_mu_sta_size]=muonRefe->ptError();
      Pat_mu_sta_phiErr[Pat_mu_sta_size]=muonRefe->phiError();
      Pat_mu_sta_etaErr[Pat_mu_sta_size]=muonRefe->etaError();
      Pat_mu_sta_d0[Pat_mu_sta_size]=muonRefe->d0();
      Pat_mu_sta_d0err[Pat_mu_sta_size]=muonRefe->d0Error();
      Pat_mu_sta_dz[Pat_mu_sta_size]=muonRefe->dz();
      Pat_mu_sta_dzerr[Pat_mu_sta_size]=muonRefe->dzError();
      Pat_mu_sta_chi2[Pat_mu_sta_size]=muonRefe->chi2();
      Pat_mu_sta_ndof[Pat_mu_sta_size]=muonRefe->ndof();
      Pat_mu_sta_nhits[Pat_mu_sta_size]=muonRefe->numberOfValidHits();
      TMatrixD cov2PAT(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov2PAT(lan,len)=muonRefe->covariance(lan,len);
        }
      }
      new((*Pat_mu_sta_CovM)[Pat_mu_sta_size])TMatrixD(cov2PAT);
      Pat_mu_sta_size++;     
    }
    /// tracker muons      
    if (muoni->isTrackerMuon()) {      
      TVector3 b(muoni->vx(),muoni->vy(),muoni->vz());
      new((*Pat_mu_trk_3vec)[Pat_mu_trk_size])TVector3(b);
      Pat_mu_trk_charge[Pat_mu_trk_size]=muoni->charge();
      TLorentzVector a=lorentzMomentum(*muoni);
      if(theDebugLevel>1)cout << "PAT Muon PT " << a.Pt() << endl;
      new((*Pat_mu_trk_4mom)[Pat_mu_trk_size])TLorentzVector(a); 
      TrackRef muonRefe = muoni->innerTrack();
      if ( muonRefe.isNull() ) { cout << "muon reference is null!"<<endl; } 
      Pat_mu_trk_ptErr[Pat_mu_trk_size]=muonRefe->ptError();
      Pat_mu_trk_phiErr[Pat_mu_trk_size]=muonRefe->phiError();
      Pat_mu_trk_etaErr[Pat_mu_trk_size]=muonRefe->etaError();
      Pat_mu_trk_d0[Pat_mu_trk_size]=muonRefe->d0();
      Pat_mu_trk_d0err[Pat_mu_trk_size]=muonRefe->d0Error();
      Pat_mu_trk_dz[Pat_mu_trk_size]=muonRefe->dz();
      Pat_mu_trk_dzerr[Pat_mu_trk_size]=muonRefe->dzError();
      Pat_mu_trk_chi2[Pat_mu_trk_size]=muonRefe->chi2();
      Pat_mu_trk_ndof[Pat_mu_trk_size]=muonRefe->ndof();
      Pat_mu_trk_nhits[Pat_mu_trk_size]=muonRefe->numberOfValidHits();
      TMatrixD cov3PAT(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov3PAT(lan,len)=muonRefe->covariance(lan,len);
        }
      }
      new((*Pat_mu_trk_CovM)[Pat_mu_trk_size])TMatrixD(cov3PAT);
      Pat_mu_trk_size++;     
    }
    /// calorimeter muons      
    if (muoni->isCaloMuon()) {      
      TVector3 b(muoni->vx(),muoni->vy(),muoni->vz());
      new((*Pat_mu_cal_3vec)[Pat_mu_cal_size])TVector3(b);
      Pat_mu_cal_charge[Pat_mu_cal_size]=muoni->charge();
      TLorentzVector a=lorentzMomentum(*muoni);
      if(theDebugLevel>1)cout << "PAT Muon PT " << a.Pt() << endl;
      new((*Pat_mu_cal_4mom)[Pat_mu_cal_size])TLorentzVector(a);
      TrackRef muonRefe = muoni->innerTrack();
      if ( muonRefe.isNull() ) { cout << "muon reference is null!"<<endl; }  
      Pat_mu_cal_ptErr[Pat_mu_cal_size]=muonRefe->ptError();
      Pat_mu_cal_phiErr[Pat_mu_cal_size]=muonRefe->phiError();
      Pat_mu_cal_etaErr[Pat_mu_cal_size]=muonRefe->etaError();
      Pat_mu_cal_d0[Pat_mu_cal_size]=muonRefe->d0();
      Pat_mu_cal_d0err[Pat_mu_cal_size]=muonRefe->d0Error();
      Pat_mu_cal_dz[Pat_mu_cal_size]=muonRefe->dz();
      Pat_mu_cal_dzerr[Pat_mu_cal_size]=muonRefe->dzError();
      Pat_mu_cal_chi2[Pat_mu_cal_size]=muonRefe->chi2();
      Pat_mu_cal_ndof[Pat_mu_cal_size]=muonRefe->ndof();
      Pat_mu_cal_nhits[Pat_mu_cal_size]=muonRefe->numberOfValidHits();
      TMatrixD cov4PAT(5,5);
      for (int lan=0;lan<5;lan++ ) {
        for ( int len=0;len<5;len++ ) {
          cov4PAT(lan,len)=muonRefe->covariance(lan,len);
        }
      }
      new((*Pat_mu_cal_CovM)[Pat_mu_cal_size])TMatrixD(cov4PAT);
      Pat_mu_cal_size++;     
    }

 
  }
} 

//////////////////////////////////////////////////////////////
// Get Primary vertex info
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillPrimaryVertex(const edm::Event &iEvent) {
  Reco_PriVtx_size=0;

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(thePrimaryVertexLabel, privtxs);

  for ( reco::VertexCollection::const_iterator vtx=privtxs->begin();
	vtx!=privtxs->end()&&Reco_PriVtx_size<Max_PriVtx_size; 
	++vtx) { 
    TVector3 vertex(0.0,0.0,0.0);
    vertex.SetXYZ(vtx->position().x(),vtx->position().y(),vtx->position().z());
    new((*Reco_PriVtx_3vec)[Reco_PriVtx_size])TVector3(vertex);
    Reco_PriVtx_xxE[Reco_PriVtx_size]=vtx->covariance(0,0);
    Reco_PriVtx_yyE[Reco_PriVtx_size]=vtx->covariance(1,1);
    Reco_PriVtx_zzE[Reco_PriVtx_size]=vtx->covariance(2,2);
    Reco_PriVtx_yxE[Reco_PriVtx_size]=vtx->covariance(1,0);
    Reco_PriVtx_zyE[Reco_PriVtx_size]=vtx->covariance(2,1);
    Reco_PriVtx_zxE[Reco_PriVtx_size]=vtx->covariance(2,0);
    Reco_PriVtx_trkSize[Reco_PriVtx_size]=vtx->tracksSize();
    Reco_PriVtx_chi2[Reco_PriVtx_size]=vtx->chi2();
    Reco_PriVtx_ndof[Reco_PriVtx_size]=vtx->ndof();
    if(theDebugLevel>1) cout<<Reco_PriVtx_size<<" Primary Vtx x="<<vtx->position().x()<<" y="<<vtx->position().y()<<endl;
    Reco_PriVtx_size++;
  }      

  if (theStoreOniaRadiation) { 
    reco::VertexCollection::const_iterator vtx=privtxs->begin();  
    Reco_PriVtx_1st_trkSize=0;
    for (Vertex::trackRef_iterator iter = vtx->tracks_begin(); 
         iter!=vtx->tracks_end();
         ++iter ) {
      Reco_PriVtx_1st_trkindex[Reco_PriVtx_1st_trkSize]=(*iter).key();
      Reco_PriVtx_1st_trkSize++;
    } 
  } 

}

//////////////////////////////////////////////////////////////
// Get Primary vertex info without muons tracks
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillPrimaryVertex2m(const edm::Event &iEvent) {
  Reco_PriVtx_size=0;
  bool ismuon=false;
  vector<TransientTrack> ttcoll;
  AdaptiveVertexFitter avf;
  GlobalError errv1;

  ttcoll.clear();
  vector<double> glmuontag;
  glmuontag.clear();
  for (reco::MuonCollection::const_iterator muoni = theGlobalMuons.begin();
         muoni != theGlobalMuons.end()&&Reco_mu_glb_size<Max_mu_size; 
         muoni++) {
    
    TrackRef innTrack = muoni->innerTrack();

    glmuontag.push_back(innTrack->p());
  }
  vector<double> tkmuontag;
  tkmuontag.clear();
  for (reco::MuonCollection::const_iterator trkmuoni = theTrkMuons.begin();
       trkmuoni != theTrkMuons.end()&&Reco_mu_trk_size<Max_mu_size; 
         trkmuoni++) {
    
      TrackRef innTrack = trkmuoni->innerTrack();

      tkmuontag.push_back(innTrack->p());
  }

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(thePrimaryVertexLabel, privtxs);

  for ( reco::VertexCollection::const_iterator vtx=privtxs->begin();
	vtx!=privtxs->end()&&Reco_PriVtx_size<Max_PriVtx_size; 
	++vtx) { 

    for (Vertex::trackRef_iterator t = vtx->tracks_begin();t !=vtx->tracks_end(); t++) {
      ismuon=false;
      double pref=(**t).p();

      for (unsigned int i=0;i<glmuontag.size();++i){if (pref==glmuontag[i]) ismuon=true;}
      for (unsigned int i=0;i<tkmuontag.size();++i){if (pref==tkmuontag[i]&&(glmuontag.size()<2)) ismuon=true;}
 
      if (!ismuon) {TransientTrack t1   = (*theB).build(&(**t));ttcoll.push_back(t1);}
 
    }
    TVector3 vertex(0.0,0.0,0.0);
    TransientVertex tv1;

    if (ttcoll.size()>=2){

      if (thePrimaryVertexLabel.label()=="offlinePrimaryVerticesWithBS"){
	reco::BeamSpot BSpot;
	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
	iEvent.getByLabel(theBeamSpotLabel ,recoBeamSpotHandle);
	if (recoBeamSpotHandle.isValid()){
	  BSpot = *recoBeamSpotHandle;
	  tv1 = avf.vertex(ttcoll,BSpot);
	}
	
      }
      else tv1 = avf.vertex(ttcoll);

      if (tv1.isValid()) {
	errv1 = tv1.positionError();

	vertex.SetXYZ(tv1.position().x(),tv1.position().y(),tv1.position().z());
	new((*Reco_PriVtx_3vec)[Reco_PriVtx_size])TVector3(vertex);
	Reco_PriVtx_xxE[Reco_PriVtx_size]=errv1.cxx();
	Reco_PriVtx_yyE[Reco_PriVtx_size]=errv1.cyy();
	Reco_PriVtx_zzE[Reco_PriVtx_size]=errv1.czz();
	Reco_PriVtx_yxE[Reco_PriVtx_size]=errv1.cyx();
	Reco_PriVtx_zyE[Reco_PriVtx_size]=errv1.czy();
	Reco_PriVtx_zxE[Reco_PriVtx_size]=errv1.czx();
	Reco_PriVtx_trkSize[Reco_PriVtx_size]=ttcoll.size();
	Reco_PriVtx_chi2[Reco_PriVtx_size]=tv1.totalChiSquared();
	Reco_PriVtx_ndof[Reco_PriVtx_size]=tv1.degreesOfFreedom();
	
	if(theDebugLevel>1) cout<<Reco_PriVtx_size<<" Primary Vtx x="<<tv1.position().x()<<" y="<<tv1.position().y()<<endl;
      }
      Reco_PriVtx_size++;
    }
  }      
	
  if (theStoreOniaRadiation) { 
    reco::VertexCollection::const_iterator vtx=privtxs->begin();  
    Reco_PriVtx_1st_trkSize=0;
    for (Vertex::trackRef_iterator iter = vtx->tracks_begin(); 
         iter!=vtx->tracks_end();
         ++iter ) {
      Reco_PriVtx_1st_trkindex[Reco_PriVtx_1st_trkSize]=(*iter).key();
      Reco_PriVtx_1st_trkSize++;
    } 
  } 

}




//////////////////////////////////////////////////////////////
// Fill BeamSpot info
///////////////////////////////////////////////////////////////
void Onia2MuMu::fillBeamSpot(const edm::Event &iEvent) {

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel(theBeamSpotLabel ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle; 

  Reco_BeamSpot_x=bs.x0();
  Reco_BeamSpot_y=bs.y0();
  Reco_BeamSpot_z=bs.z0();
  Reco_BeamSpot_xxE=bs.covariance(0,0);
  Reco_BeamSpot_yyE=bs.covariance(1,1);
  Reco_BeamSpot_zzE=bs.covariance(2,2);
  Reco_BeamSpot_yxE=bs.covariance(1,0);
  Reco_BeamSpot_zyE=bs.covariance(2,1);
  Reco_BeamSpot_zxE=bs.covariance(2,0);
  if(theDebugLevel>1) cout<<"Beam Spot x="<<bs.x0()<<" y="<<bs.y0()<<endl;


}

//////////////////////////////////////////////////////////////
// fills Onia candidate block
///////////////////////////////////////////////////////////////
void Onia2MuMu::findOniaCategories(const edm::Event &iEvent) {

  Reco_QQ_size=0;
  Reco_Chic_size=0;
  Reco_Bp_size=0;

  Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(thePrimaryVertexLabel, privtxs);
  VertexCollection::const_iterator privtx;

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel(theBeamSpotLabel ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;

  TVector3 vperp2;
  if ( theBeamSpotFlag ) {
    vperp2.SetXYZ(bs.x0(), bs.y0(), 0);
  }
  else {
    if ( privtxs->begin() != privtxs->end() ) {
      privtx=privtxs->begin();
      vperp2.SetXYZ(privtx->position().x(), privtx->position().y(), 0);
    }
    else {
      vperp2.SetXYZ(0, 0, 0);
    }
  }
  
  MuonCollection::const_iterator nmuon1;
  MuonCollection::const_iterator nmuon2;
  MuonCollection::const_iterator tmuon1;
  MuonCollection::const_iterator tmuon2;
  CaloMuonCollection::const_iterator cmuon1;
  CaloMuonCollection::const_iterator cmuon2;

  // int oniacato;
  TrackRef muon1;
  TrackRef muon2;
  int m1=0;
  int m2g=0;
  int m2t=0;
  int m2c=0;
  
  for ( nmuon1 = theGlobalMuons.begin(); nmuon1 != theGlobalMuons.end(); nmuon1++) {
    for( nmuon2 = nmuon1+1; nmuon2 != theGlobalMuons.end()&&Reco_QQ_size<3000; nmuon2++ ) {
      muon1 = nmuon1->innerTrack();    
      muon2 = nmuon2->innerTrack(); 
      m2g++;
      if (theOniaMaxCat >= 0) fillOniaMuMuTracks(muon1, m1, muon2, m1+m2g, vperp2, 0);
    }
    if ( theStoreTRKMuonFlag ) {
      for( tmuon2 = theTrkMuons.begin(); tmuon2 != theTrkMuons.end()&&Reco_QQ_size<3000; tmuon2++ ) {
	muon1 = nmuon1->innerTrack();    
	muon2 = tmuon2->innerTrack(); 
	 if (theOniaMaxCat >= 1) fillOniaMuMuTracks(muon1, m1, muon2, m2t, vperp2, 1);
	m2t++;
      }
    }
    if ( theStoreCALMuonFlag ) {
      for( cmuon2 = theCaloMuons.begin(); cmuon2 != theCaloMuons.end()&&Reco_QQ_size<3000; cmuon2++ ) {
	muon1 = nmuon1->innerTrack();    
	muon2 = cmuon2->track(); 
	 if (theOniaMaxCat >= 3) fillOniaMuMuTracks(muon1, m1, muon2, m2c, vperp2, 3);
	m2c++;
      }
    }
    m1++;
  }
  m1=0; m2t=0; m2c=0; 
  if ( theStoreTRKMuonFlag ) {
    for ( tmuon1 = theTrkMuons.begin(); tmuon1 != theTrkMuons.end(); tmuon1++) {
      for( tmuon2 = tmuon1+1; tmuon2 != theTrkMuons.end()&&Reco_QQ_size<3000; tmuon2++ ) {
	muon1 = tmuon1->innerTrack();    
	muon2 = tmuon2->innerTrack(); 
	m2t++; 
	 if (theOniaMaxCat >= 2) fillOniaMuMuTracks(muon1, m1, muon2, m2t+m1, vperp2, 2);
      }
      if ( theStoreCALMuonFlag ) {
	for( cmuon2 = theCaloMuons.begin(); cmuon2 != theCaloMuons.end()&&Reco_QQ_size<3000; cmuon2++ ) {
	  muon1 = tmuon1->innerTrack();    
	  muon2 = cmuon2->track(); 
	   if (theOniaMaxCat >= 4) fillOniaMuMuTracks(muon1, m1, muon2, m2c, vperp2, 4);
	  m2c++;
	}
      }
      m1++;
    }
  }
  m1=0; m2c=0;
  if ( theStoreCALMuonFlag ) {
    for ( cmuon1 = theCaloMuons.begin(); cmuon1 != theCaloMuons.end(); cmuon1++) {
      for( cmuon2 = cmuon1+1; cmuon2 != theCaloMuons.end()&&Reco_QQ_size<3000; cmuon2++ ) {
	muon1 = cmuon1->track();    
	muon2 = cmuon2->track(); 
	m2c++; 
	 if (theOniaMaxCat >= 5) fillOniaMuMuTracks(muon1, m1, muon2, m2c+m1, vperp2, 5);
      }
      m1++;
    }
  }
}
//////////////////////////////////////////////////////////////////////////
// checks if phi is in range
//////////////////////////////////////////////////////////////////////////
double Onia2MuMu::PhiInRange(const double& phi) const {
      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}


//////////////////////////////////////////////////////////////////////////
// calculates DeltaR
//////////////////////////////////////////////////////////////////////////
template <class T, class U>
double Onia2MuMu::deltaR(const T & t, const U & u) const {
      return sqrt(pow(t.Eta()-u.Eta(),2) +pow(PhiInRange(t.Phi()-u.Phi()),2));
}

////////////////////////////////////////////////////////////////////////
//calculates also DeltaR
////////////////////////////////////////////////////////////////////////
double Onia2MuMu::calcDeltaR(double eta1, double phi1,
                                          double eta2, double phi2){
  // Calculate delta phi and delta eta
  double delPhi = fabs(phi1 - phi2);
  if (delPhi >  TMath::Pi()) delPhi = 2.0*TMath::Pi() - delPhi;
  double delEta = eta1-eta2;
  // Calculate the delta R
  double delR = sqrt(delPhi*delPhi+delEta*delEta);
  return delR;
}

////////////////////////////////////////////////////////////////////////
//calculates theta
////////////////////////////////////////////////////////////////////////
double Onia2MuMu::GetTheta( TLorentzVector & a,  TLorentzVector & b) const {
  TLorentzVector c = a+b;
  TVector3 bv=c.BoostVector();
  a.Boost(-bv);
  b.Boost(-bv);
  double theta = c.Vect().Angle(a.Vect());
  return theta;
}

////////////////////////////////////////////////////////////////////////
// print track
////////////////////////////////////////////////////////////////////////
void Onia2MuMu::printTrack(const reco::Track& muon) const {
   cout << "phi=" << muon.phi();
   cout << " eta=" << muon.eta();
   cout << " p=" << muon.p();
   cout << " chg=" << muon.charge();
   cout << " chi2/hits "<<muon.chi2() <<"/"<<muon.found()<<endl;
}


///////////////////////////////////////////////////////////////////////////////
// Returns Lorentz-vector of muon
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Onia2MuMu::lorentzMomentum(const reco::Track& muon) const {

    double preco = muon.p();
    double pxreco = muon.px();
    double pyreco = muon.py();
    double pzreco = muon.pz();

    // energy = sqrt(p^2 +m^2)
    double ereco = sqrt(preco*preco + 0.011163613);
    //double ereco = sqrt(preco*preco + 0.105658*0.105658);

    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;

}

TLorentzVector Onia2MuMu::lorentzMomentum(const pat::Muon& muon) const {

    double preco = muon.p();
    double pxreco = muon.px();
    double pyreco = muon.py();
    double pzreco = muon.pz();

    // energy = sqrt(p^2 +m^2)
    double ereco = sqrt(preco*preco + 0.011163613);
    //double ereco = sqrt(preco*preco + 0.105658*0.105658);

    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;

}

TLorentzVector Onia2MuMu::lorentzTriObj(const trigger::TriggerObject& muon) const {

    double preco = muon.p();
    double pxreco = muon.px();
    double pyreco = muon.py();
    double pzreco = muon.pz();

    // energy = sqrt(p^2 +m^2)
    double ereco = sqrt(preco*preco + 0.011163691);
    //double ereco = sqrt(preco*preco + 0.105658*0.105658);

    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;

}
///////////////////////////////////////////////////////////////////////////////
// Returns Lorentz-vector of PF photon
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Onia2MuMu::lorentzMomentum(const reco::PFCandidate & pfcl) const {

    double ereco = pfcl.energy();
    double pxreco = pfcl.px();
    double pyreco = pfcl.py();
    double pzreco = pfcl.pz();
   
    // cout << "pfcl = " << pxreco << " "  << pyreco << " "  << pzreco << " "  << ereco << endl;
    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;

}
///////////////////////////////////////////////////////////////////////////////
// Returns Lorentz-vector of track with pion mass hypothesis
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Onia2MuMu::lorentzMomentumPi(const reco::Track & tr) const {
  
    double preco = tr.p();
    double pxreco = tr.px();
    double pyreco = tr.py();
    double pzreco = tr.pz();

    // energy = sqrt(p^2 +m^2)
    double ereco = sqrt(preco*preco + 0.0194798);
    //double ereco = sqrt(preco*preco + 0.105658*0.105658);

    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;
   

}
///////////////////////////////////////////////////////////////////////////////
// Returns Lorentz-vector of track with kaon mass hypothesis
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Onia2MuMu::lorentzMomentumKa(const reco::Track & tr) const {
  
    double preco = tr.p();
    double pxreco = tr.px();
    double pyreco = tr.py();
    double pzreco = tr.pz();

    // energy = sqrt(p^2 +m^2)
    double ereco = sqrt(preco*preco + 0.24372);
    //double ereco = sqrt(preco*preco + 0.105658*0.105658);

    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;
   

}

/////////////////////////////////////////////////////////////////////
// Fill any onia category
/////////////////////////////////////////////////////////////////////
void Onia2MuMu::fillOniaMuMuTracks(TrackRef muon1, int m1, TrackRef muon2, int m2, TVector3 vperp2, int oniacato) {
 
  if (oniacato<0  ) return;
  // if ( muon1->charge() == muon2->charge() ) continue;
  Reco_QQ_sign[Reco_QQ_size]=0;
  if ( muon1->charge() == muon2->charge() ) {
    if (theStoreWSOnia) {
      if (muon1->charge() == 1) {Reco_QQ_sign[Reco_QQ_size]=1;}
      else {Reco_QQ_sign[Reco_QQ_size]=-1;}
    }
    else return;
  }
  
  TLorentzVector mu1=lorentzMomentum(*muon1);
  TLorentzVector mu2=lorentzMomentum(*muon2);
  TLorentzVector onia = mu1 + mu2;
  Reco_QQ_type[Reco_QQ_size]=oniacato;
  if (!theUseKinFit) new((*Reco_QQ_4mom)[Reco_QQ_size])TLorentzVector(onia);

  Reco_QQ_DeltaR[Reco_QQ_size]=deltaR(mu1, mu2);
  Reco_QQ_s[Reco_QQ_size] = pow((muon1->d0()/muon1->d0Error()),2)+pow((muon2->d0()/muon2->d0Error()),2);
  // Reco_QQ_absD0Diff[Reco_QQ_size] = fabs(muon1->d0() - muon2->d0());
  // Reco_QQ_absDzDiff[Reco_QQ_size] = fabs(muon1->dz() - muon2->dz());

  if ( muon1->charge() == 1 ) {
    Reco_QQ_mupl[Reco_QQ_size]=m1;
    Reco_QQ_mumi[Reco_QQ_size]=m2;
    Reco_QQ_cosTheta[Reco_QQ_size]=cos(GetTheta(mu1, mu2));
  }
  else {
    Reco_QQ_mupl[Reco_QQ_size]=m2;
    Reco_QQ_mumi[Reco_QQ_size]=m1;
    Reco_QQ_cosTheta[Reco_QQ_size]=cos(GetTheta(mu2, mu1));
  }
  if (oniacato == 1 || oniacato == 3 || oniacato == 4 || mu1.Perp() > mu2.Perp()) {   // different muon categories 
    Reco_QQ_muhpt[Reco_QQ_size]=m1;
    Reco_QQ_mulpt[Reco_QQ_size]=m2;
  } else {
    Reco_QQ_muhpt[Reco_QQ_size]=m2;
    Reco_QQ_mulpt[Reco_QQ_size]=m1;
  }
  
  TransientTrack ttkp1   = (*theB).build(&(*muon1));
  TransientTrack ttkp2   = (*theB).build(&(*muon2));

  if (theUseKinFit) {
    // NUOVO VERTICE
    KinematicParticleFactoryFromTransientTrack pFactory;

    // The mass of a muon and the insignificant mass sigma to avoid singularities in the covariance matrix.
    ParticleMass muon_mass = 0.1056583;
    float muon_sigma = 0.0000000001;

    float chi = 0.;
    float ndf = 0.;

    // making particles
    vector<RefCountedKinematicParticle> allParticles;
    allParticles.push_back(pFactory.particle (ttkp1,muon_mass,chi,ndf,muon_sigma));
    allParticles.push_back(pFactory.particle (ttkp2,muon_mass,chi,ndf,muon_sigma));

    // fit to the vertex
    KinematicParticleVertexFitter kvFitter;
    RefCountedKinematicTree myTree = kvFitter.fit(allParticles);
    // cout << "Global vertex fit done\n";

    if (myTree->isValid()) {
      myTree->movePointerToTheTop();
      RefCountedKinematicParticle newQQ = myTree->currentParticle();
      RefCountedKinematicVertex QQVertex = myTree->currentDecayVertex();
      
      Reco_QQ_VtxIsVal[Reco_QQ_size]=true;
      TVector3 vtx(0.0,0.0,0.0);
      GlobalPoint v = QQVertex->vertexState().position();
      GlobalError err = QQVertex->vertexState().error();
      vtx.SetXYZ(v.x(),v.y(),v.z());
      new((*Reco_QQ_Vtx)[Reco_QQ_size])TVector3(vtx);
      
      Reco_QQ_VxxE[Reco_QQ_size]=err.cxx();
      Reco_QQ_VyyE[Reco_QQ_size]=err.cyy();
      Reco_QQ_VzzE[Reco_QQ_size]=err.czz();
      Reco_QQ_VyxE[Reco_QQ_size]=err.cyx();
      Reco_QQ_VzxE[Reco_QQ_size]=err.czx();
      Reco_QQ_VzyE[Reco_QQ_size]=err.czy();
      Reco_QQ_lxy[Reco_QQ_size]= v.perp();
      Reco_QQ_lxyErr[Reco_QQ_size]= sqrt(err.rerr(v));
      Reco_QQ_normChi2[Reco_QQ_size]= newQQ->chiSquared()/newQQ->degreesOfFreedom();
      Reco_QQ_probChi2[Reco_QQ_size]= TMath::Prob(newQQ->chiSquared(), (int)newQQ->degreesOfFreedom());
      TLorentzVector oniaPostFit(0.0,0.0,0.0,0.0);
      oniaPostFit.SetPxPyPzE(newQQ->currentState().globalMomentum().x(),
			     newQQ->currentState().globalMomentum().y(),
			     newQQ->currentState().globalMomentum().z(),
			     sqrt(pow(newQQ->currentState().mass(),2) +
		     pow(newQQ->currentState().globalMomentum().mag(),2)));
      
      new((*Reco_QQ_4mom)[Reco_QQ_size])TLorentzVector(oniaPostFit);
      TVector3 pperp(oniaPostFit.Px(), oniaPostFit.Py(), 0);
      TVector3 vperp1(v.x(), v.y(), 0);
      TVector3 vperp = vperp1 - vperp2;
      double cosAlpha = vperp.Dot(pperp)/(vperp.Perp()*pperp.Perp());
      double ctau = vperp.Perp()*cosAlpha*oniaMass/onia.Perp();
      Reco_QQ_cosAlpha[Reco_QQ_size]= cosAlpha;
      Reco_QQ_ctau[Reco_QQ_size]= ctau;
    } else {
      TVector3 vtx(-1,-1,-1);
      new((*Reco_QQ_Vtx)[Reco_QQ_size])TVector3(vtx);
      Reco_QQ_VxxE[Reco_QQ_size]=-1;
      Reco_QQ_VyyE[Reco_QQ_size]=-1;
      Reco_QQ_VzzE[Reco_QQ_size]=-1;
      Reco_QQ_VyxE[Reco_QQ_size]=-1;
      Reco_QQ_VzxE[Reco_QQ_size]=-1;
      Reco_QQ_VzyE[Reco_QQ_size]=-1;
      Reco_QQ_lxy[Reco_QQ_size]= -1;
      Reco_QQ_lxyErr[Reco_QQ_size]= -1;
      Reco_QQ_normChi2[Reco_QQ_size]= -1;
      Reco_QQ_probChi2[Reco_QQ_size]= -1;
      new((*Reco_QQ_4mom)[Reco_QQ_size])TLorentzVector(onia);
      Reco_QQ_cosAlpha[Reco_QQ_size]= -2;
      Reco_QQ_ctau[Reco_QQ_size]= -100;
    }
  } else {
    vector<TransientTrack> t_tks;
    t_tks.push_back(ttkp1);
    t_tks.push_back(ttkp2);

    KalmanVertexFitter kvf;
    TransientVertex tv = kvf.vertex(t_tks);
    
    if (  tv.isValid() ) {
      Reco_QQ_VtxIsVal[Reco_QQ_size]=true;
      GlobalPoint v = tv.position();
      GlobalError err = tv.positionError();
      TVector3 vtx(0.0,0.0,0.0);
      vtx.SetXYZ(v.x(),v.y(),v.z());
      new((*Reco_QQ_Vtx)[Reco_QQ_size])TVector3(vtx);
      
      Reco_QQ_VxxE[Reco_QQ_size]=err.cxx();
      Reco_QQ_VyyE[Reco_QQ_size]=err.cyy();
      Reco_QQ_VzzE[Reco_QQ_size]=err.czz();
      Reco_QQ_VyxE[Reco_QQ_size]=err.cyx();
      Reco_QQ_VzxE[Reco_QQ_size]=err.czx();
      Reco_QQ_VzyE[Reco_QQ_size]=err.czy();
      Reco_QQ_lxy[Reco_QQ_size]= v.perp();
      Reco_QQ_lxyErr[Reco_QQ_size]= sqrt(err.rerr(v));
      Reco_QQ_normChi2[Reco_QQ_size]= tv.normalisedChiSquared();
      Reco_QQ_probChi2[Reco_QQ_size]= TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
      TVector3 pperp(onia.Px(), onia.Py(), 0);
      TVector3 vperp1(v.x(), v.y(), 0);
      TVector3 vperp = vperp1 - vperp2;
      double cosAlpha = vperp.Dot(pperp)/(vperp.Perp()*pperp.Perp());
      double ctau = vperp.Perp()*cosAlpha*oniaMass/onia.Perp();
      Reco_QQ_cosAlpha[Reco_QQ_size]= cosAlpha;
      Reco_QQ_ctau[Reco_QQ_size]= ctau;
    } else {
      Reco_QQ_VtxIsVal[Reco_QQ_size]=false;
      TVector3 vtx(-1,-1,-1);
      new((*Reco_QQ_Vtx)[Reco_QQ_size])TVector3(vtx);
      Reco_QQ_VxxE[Reco_QQ_size]=-1;
      Reco_QQ_VyyE[Reco_QQ_size]=-1;
      Reco_QQ_VzzE[Reco_QQ_size]=-1;
      Reco_QQ_VyxE[Reco_QQ_size]=-1;
      Reco_QQ_VzxE[Reco_QQ_size]=-1;
      Reco_QQ_VzyE[Reco_QQ_size]=-1;
      Reco_QQ_lxy[Reco_QQ_size]= -1;
      Reco_QQ_lxyErr[Reco_QQ_size]= -1;
      Reco_QQ_normChi2[Reco_QQ_size]= -1;
      Reco_QQ_probChi2[Reco_QQ_size]= -1;
      Reco_QQ_cosAlpha[Reco_QQ_size]= -2;
      Reco_QQ_ctau[Reco_QQ_size]= -100;
    }
  }
  // Here implement onia+gamma / onia+track combos
  if (theStoreChicFlag && Reco_QQ_sign[Reco_QQ_size] == 0 && oniacato <= (int)maxCatToStoreChic) {
 
    int countClus = 0;
    for(PFCandidateCollection::const_iterator itePFC = pfClusters.begin();
	itePFC != pfClusters.end()&&Reco_Chic_size<3000;
	++itePFC) {
      PFCandidate itPFC = *itePFC;
      TLorentzVector gamma = lorentzMomentum(itPFC);
      TLorentzVector chic = onia + gamma;
      double deltaM = chic.M() - onia.M();
      if (deltaM < 5.0) {
	new((*Reco_Chic_4mom)[Reco_Chic_size])TLorentzVector(chic);
	Reco_Chic_DeltaM[Reco_Chic_size] = deltaM;
	Reco_Chic_OniaDaug[Reco_Chic_size] = Reco_QQ_size;
	Reco_Chic_GammaDaug[Reco_Chic_size] = countClus;
	countClus++;
	Reco_Chic_size++;
      }
    }
  }
  if (theStoreBpFlag && Reco_QQ_sign[Reco_QQ_size]==0 && oniacato<=(int)maxCatToStoreBp) {
 
    int countTracks = 0;
    for(TrackCollection::const_iterator itTr = noMuonTracks.begin();
	itTr != noMuonTracks.end()&&Reco_Bp_size<3000;
	++itTr) {
      TLorentzVector kappa = lorentzMomentumKa(*itTr);
      TLorentzVector bplus = onia + kappa;
      double bplusM = bplus.M();
      if (bplusM < 10.0) {
	new((*Reco_Bp_4mom)[Reco_Bp_size])TLorentzVector(bplus);
        Reco_Bp_OniaDaug[Reco_Bp_size] = Reco_QQ_size;
        Reco_Bp_KDaug[Reco_Bp_size] = countTracks;
	countTracks++;
	Reco_Bp_size++;
      }
    }
  }
    
  Reco_QQ_size++;
  return;
}

/////////////////////////////////////////////////////////////////////
// Calculate invariant mass
/////////////////////////////////////////////////////////////////////
double Onia2MuMu::invMass(const reco::Track& lhs, const reco::Track& rhs) const {

  return (lorentzMomentum(lhs) + lorentzMomentum(rhs)).Mag();

}

/////////////////////////////////////////////////////////////////////
// Check if the PDG of a long-lived B-hadron
/////////////////////////////////////////////////////////////////////
bool Onia2MuMu::isAbHadron(int pdgID) const {

  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;

}

/////////////////////////////////////////////////////////////////////
// Number of hits per muon/track
/////////////////////////////////////////////////////////////////////

// OLD WAY (non funziona con gli aoddi')
/* std::vector<unsigned int> Onia2MuMu::muonStatHits(const reco::Track& tr) {

  std::vector<unsigned int> theMuonHits;
  unsigned int nRecHitDT(0), nRecHitCSC(0), nRecHitRPC(0);

  for(trackingRecHit_iterator recHit = tr.recHitsBegin(); recHit != tr.recHitsEnd(); ++recHit){
     DetId detIdHit = (*recHit)->geographicalId();
     if (detIdHit.det() == DetId::Muon ){
       if (detIdHit.subdetId() == MuonSubdetId::DT ) nRecHitDT++;
       else if (detIdHit.subdetId() == MuonSubdetId::CSC ) nRecHitCSC++;
       else if (detIdHit.subdetId() == MuonSubdetId::RPC ) nRecHitRPC++;
     }
   }

  theMuonHits.push_back(nRecHitDT); 
  theMuonHits.push_back(nRecHitCSC);
  theMuonHits.push_back(nRecHitRPC);
  return theMuonHits; 
}

std::vector<unsigned int> Onia2MuMu::trackStatHits(const reco::Track& tr) {
  std::vector<unsigned int> theTrackHits;
  unsigned int nRecHitPixelB(0), nRecHitPixelE(0), nRecHitStrip(0), nhits(0);
  unsigned int firstlayer(0), firstdisk(0);
  for(trackingRecHit_iterator recHit = tr.recHitsBegin(); recHit != tr.recHitsEnd(); ++recHit){
     DetId detIdHit = (*recHit)->geographicalId();
     if (detIdHit.det() == DetId::Tracker ){
       if ((detIdHit.subdetId() == StripSubdetector::TIB)||(detIdHit.subdetId() == StripSubdetector::TOB)||(detIdHit.subdetId() == StripSubdetector::TEC)||(detIdHit.subdetId() == StripSubdetector::TIB)) nRecHitStrip++;
	   else if (detIdHit.subdetId() == PixelSubdetector::PixelBarrel) {nhits++;nRecHitPixelB++;if (nhits==1){PixelBarrelName PixB(detIdHit);firstlayer=PixB.layerName();}}
	   else if (detIdHit.subdetId() == PixelSubdetector::PixelEndcap) {nhits++;nRecHitPixelE++;if (nhits==1){PixelEndcapName PixE(detIdHit);firstdisk =PixE.diskName();}}
     }

   }

  theTrackHits.push_back(nRecHitStrip); 
  theTrackHits.push_back(nRecHitPixelB);
  theTrackHits.push_back(nRecHitPixelE);
  if (firstlayer>0) {theTrackHits.push_back(firstlayer);theTrackHits.push_back(1);}
  if (firstdisk>0)  {theTrackHits.push_back(firstdisk);theTrackHits.push_back(2);}
  if (firstlayer==0&&firstdisk==0) {theTrackHits.push_back(firstlayer);theTrackHits.push_back(0);}

  return theTrackHits; 
  } */

// SMART WAY (Boris)
std::vector<unsigned int> Onia2MuMu::trackHits(const reco::Track& tr) {

  std::vector<unsigned int> theTrackHits;
  const reco::HitPattern& p = tr.hitPattern();

  theTrackHits.push_back(p.numberOfValidMuonDTHits()); 
  theTrackHits.push_back(p.numberOfValidMuonCSCHits());
  theTrackHits.push_back(p.numberOfValidMuonRPCHits());
  theTrackHits.push_back(p.numberOfValidStripHits()); 
  theTrackHits.push_back(p.numberOfValidPixelBarrelHits());
  theTrackHits.push_back(p.numberOfValidPixelEndcapHits());

  // do not loop over the hits of the track
  uint32_t firsthit = p.getHitPattern(0);
    
  // if the hit is valid and in pixel barrel... etc. etc.
  if (p.validHitFilter(firsthit) && p.pixelBarrelHitFilter(firsthit)) {
    theTrackHits.push_back(p.getLayer(firsthit));
    theTrackHits.push_back(1);
  } else if (p.validHitFilter(firsthit) && p.pixelEndcapHitFilter(firsthit)) {
    theTrackHits.push_back(p.getLayer(firsthit));
    theTrackHits.push_back(2);
  } else {
    theTrackHits.push_back(p.getLayer(firsthit));
    theTrackHits.push_back(0);
  }

  return theTrackHits; 
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(Onia2MuMu);
