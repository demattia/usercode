//////////////////////////////////////////////////////////////////////////////
//
// AlgValidator.cc
// Code extracted from ProbAssoc.cc, 29/04/08 R.Casagrande
//
//  
//
// --------------------------------------------------------------------------------
//
// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/AlgValidator.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"

//======ROBERTO========================
// For Tracks
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"
// For jet
#include "AnalysisExamples/AnalysisObjects/interface/SimpleCaloJet.h"

// For Vertices
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// For the offline jets and corrections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// For the b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

//For RecoVertex
#include "DataFormats/VertexReco/interface/Vertex.h"

//======ROBERTO========================

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>

//======ROBERTO========================

//Implemented functions
//-----------------------
#include "AnalysisExamples/AnalysisClasses/interface/SortingDef.h"

#include "AnalysisExamples/AnalysisClasses/interface/DzAssociator.h"

#include "AnalysisExamples/AnalysisClasses/interface/D2Associator.h"

#include "AnalysisExamples/AnalysisClasses/interface/WAverager.h"


// #include <TF1.h>

//Implemented Classes
//------------------------
#include "AnalysisExamples/AnalysisClasses/interface/RejectionAlg.h"

#include "AnalysisExamples/AnalysisClasses/interface/VerticesJetCounter.h"

#include "AnalysisExamples/AnalysisClasses/interface/ProbCal.h"

#include "AnalysisExamples/AnalysisClasses/interface/MinimumDzVtx.h"

#include "AnalysisExamples/AnalysisClasses/interface/ProbEval.h"



//======ROBERTO========================
// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------


// Constructors and destructor
// ---------------------------
AlgValidator::AlgValidator(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  
  //--------------------
  // Jet
  //--------------------
  
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),

  //=====ROBERTO================================
   
  //----------------------------------
  // Simulated Tracks and Vertices
  //----------------------------------
  
  //Simulated vertices
  simVtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx" ) ),
  
  //Simulated tracks
  simTkLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimTk" ) ),
  
  simPUVtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimPUVtx" ) ),
 
  //----------------------------------
  // Reco Tracks and Vertices
  //----------------------------------  

  //Reco tracks impact parameter
  impactParameterTagInfos( iConfig.getUntrackedParameter<std::string>( "impactParameterTagInfos" ) ),

  //Reco vertexs
  vtxSample_( iConfig.getUntrackedParameter<edm::InputTag>( "vtxSample" ) ),

  //   genJetsLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "genJets" ) ),

  SIM_( iConfig.getUntrackedParameter<bool> ( "SIM" ) ), 

  //=====ROBERTO=============================== 

  QCD_( iConfig.getUntrackedParameter<bool> ( "QCD" ) ),
  OutputEffFileName( iConfig.getUntrackedParameter<std::string>( "OutputEffFileName" ) )

{
   
  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;
  
  // File for output histograms
  // --------------------------
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd(); 

  // Histograms
  // ----------
  //=====ROBERTO================
  
  //Histograms for Njet   
  H_NJet_withTracks_ =  new TH1D("H_NJet_withTracks", "# of jet (Et>25) with assigned tracks " ,40,0.,40.);
  H_NJet_withNoTracks_ =  new TH1D("H_NJet_withNoTracks", "# of jet (Et>25) without assigned tracks " ,20,0.,20.);
  
  //Histograms for Dz jet - vertices 
  H_Vmin_ =  new TH1D("H_Vmin", "vertex # assigned to a jet (W Avg) " ,50,0.,50.);
  //Histogram Z jet 
  H_ZJetWAvg_ = new TH1D("H_ZJetWAvg", "Z jet (W Average Z of tracks)" ,40,-20.,20.);
  H_ZJetWAvgError_ = new TH1D("H_ZJetWAvgError", "Error Z jet (W Average Z of tracks)" ,200,0.,0.2);
  H_DzJetVtxWAvg_ =  new TH1D("H_DzJetVtxWAvg", "#DeltaZ jet - vertex (W Avg) " ,800,-4.,4.);
  H_DzJetVtxWAvg2Best_ =  new TH1D("H_DzJetVtxWAvg2Best", "#DeltaZ second best jet - vertex(W Avg) " ,80,-4.,4.);
  H_DzJetVtxWAvg_vs_2Best_ =  new TH2D("H_DzJetVtxWAvg_vs_2Best","#DeltaZ vs #DeltaZ second best jet - vertex(W Avg)" ,80,0.,4.,80,0.,4.);
   
  //Histograms of id vertex assigned to a jet after rejection algorithm
  H_idRJAlg_ = new TH1D("H_idRJAlg", "vertex # assigned to a jet (RJAlg) " ,50,0.,50.);
  //Histograms of Z jet and error after the rejection algorithm (13 sigma)
  H_ZJet_RJAlg_ = new TH1D("H_ZJet_RJAlg", "Z jet (W Average Z of tracks after rejection alg)" ,40,-20.,20.);
  H_ZErrorJet_RJAlg_ = new TH1D("H_ZErrorJet_RJAlg", "#sigmaZ jet (W Average Z of tracks after rejection alg)" ,200,0.,0.2);
  H_DzJetVtx_RJAlg_ =  new TH1D("H_DzJetVtx_RJAlg", "#DeltaZ jet after rejection - vertex " ,800,-4.,4.);
  H_Dz2JetVtx_RJAlg_ =  new TH1D("H_Dz2JetVtx_RJAlg", "#DeltaZ jet after rejection - vertex " ,80,-4.,4.);
  H_DzJetVtx_RJAlg_vs_2Best_ =  new TH2D("H_DzJetVtxRJAlg_vs_2Best","#DeltaZ vs #DeltaZ second best jet - vertex RJAlg" ,80,0.,4.,80,0.,4.);

  //Histograms of fraction between probability of best vertex association and second best 
  H_P1P2_= new TH1D("H_P1P2","P1 / P2",500,0., 5.);
  
  //Histograms of id vertex assigned to a jet after SAAlg
  H_idVtxSAAlg_ = new TH1D("H_idVtxSAAlg","idVtxSAAlg ", 50, 0, 50);
  //Histograms of Dz jet - vertex after SAAlg
  H_DzAfterSAAlg_ = new TH1D("H_DzAfterSAAlg", "#DeltaZ jet - vertex after SAAlg " ,800,-4.,4.);

  //Histograms of difference of log of probability of best vertex association and second best vs 
  //difference of log of probability considering the vertexes jet multeplicity
  H_pVsp_ = new TH2D("H_pVsp","log(P1)- log(P2) vs (log(P1) + log(P(m+1,n))) - (log(P2) + log(P(m, n+1)))" ,200, 0.,20.,400,-20.,20.);

  //Histograms of # of vertexes with a specified jet multeplicity
  H_nVtxMap_  = new TH1D("H_nVtxMap","nVtxMap FAAlg", 20, 0, 20);
  H_nVtxMapRJAlg_ = new TH1D("H_nVtxMapRJAlg","nVtxMapRJAlg", 20, 0, 20);
  H_nVtxMapSAAlg_ = new TH1D("H_nVtxMapSAAlg","nVtxMapSAAlg", 20, 0, 20);
  H_nVtxMapProbEval_ = new TH1D("H_nVtxMapProbEval","nVtxMapProbEval", 20, 0, 20);

  //Histograms of probability  of jet assigned to vertexes configuration
  H_ProbRJAlg_ =  new TH1D("H_ProbRJAlg","Probability of event after RJAlg ", 1500, 0., 0.15);
  H_ProbAfterFAAlg_ =  new TH1D("H_ProbAfterFAAlg","Probability of event after first association ", 1500, 0., 0.15);
  H_ProbAfterSAAlg_ =  new TH1D("H_ProbAfterSAAlg","Probability of event after SAAlg ", 1500, 0., 0.15);
 
  //Histograms of id vertex assigned to a jet after first association
  H_idVtxFAAlg_ = new TH1D("H_idVtxFAAlg","idVtxFAAlg ", 50, 0, 50);
 
  //Histograms of # Tracks for jet with #DeltaZ best 
  H_nTrack_RJAlg_ =  new TH1D("H_nTrack_RJAlg","# Tracks for jet with #DeltaZ best < 5#sigma", 50, 0, 50);

  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after RJAlg
  H_nJMaxVsNjTot_RJAlg_ =  new TH2D("H_nJMaxVsNjTot_RJAlg","# jet for vertex with max multeplicity  vs # jet tot after RJAlg" ,20,0,20,20,0,20);
  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after SAAlg
  H_nJMaxVsNjTot_SAAlg_ = new TH2D("H_nJMaxVsNjTot_SAAlg","# jet for vertex with max multeplicity  vs # jet tot after SAAlg" ,20,0,20,20,0,20);

  //Histograms of # Tracks for jet with #DeltaZ best vs z jet (error)
  H_nTkRJAlgVSzJetError_ =  new TH2D("H_nTkRJAlgVSzJetError","# Tracks for jet after RJAlg vs z jet error" ,20,0,20,20,0.,2.);
  H_nTkRJAlgVSzJet_ = new TH2D("H_nTkRJAlgVSzJet","# Tracks for jet with after RJAlg best vs z jet" ,20,0.,20.,200,0.,20.);

  //Histograms of #jet with and without tracks (only cinematics cuts)
  H_nJet_ = new TH1D("H_nJet", "# of jet (Et>25)" ,40,0.,40.);
 
  //Histograms of #jet before, after RJAlg and SAAlg (with tracks)
  H_nJMax_SAAlg_  =  new TH1D("H_nJMax_SAAlg", "# of jet (Et>25) assigned after SAAlg" ,20,0.,20.);
  H_nJMax_RJAlg_  =  new TH1D("H_nJMax_RJAlg", "# of jet (Et>25) assigned after RJAlg" ,20,0.,20.);
  H_nJMax_BeforeAlg_  =  new TH1D("H_nJMax_BeforeAlg", "# of jet (Et>25) assigned BeforeAlg" ,20,0.,20.);
  H_nJMax_ProbEval_ =  new TH1D("H_nJMax_ProbEval", "# of jet (Et>25) assigned after ProbEval" ,20,0.,20.);

  //Histograms of Minimum Dz between vertexes
  H_VtxMinDz_ = new TH1D("H_VtxMinDz", "Minimum #DeltaZ between vertexes" ,150,0.,15.);

  H_ProfileVtxVsNjBeforeAlg_ = new TProfile("H_ProfileVtxVsNjBeforeAlg","#DeltaZ vtx - PV vs Nj Before Algorithm",80,-20.,20.);
  H_ProfileVtxVsNjRJAlg_ = new TProfile("H_ProfileVtxVsNjRJAlg","#DeltaZ vtx - PV vs Nj after RJAlgorithm",80,-20.,20.);
  H_ProfileVtxVsNjSAAlg_ = new TProfile("H_ProfileVtxVsNjSAAlg","#DeltaZ vtx - PV vs Nj after SAAlgorithm",80,-20.,20.);
  H_ProfileVtxVsNjProbEval_ = new TProfile("H_ProfileVtxVsNjProbEval","#DeltaZ vtx - PV vs Nj after ProbEval",80,-20.,20.);

  H_ZVtx_ = new TH1D("H_ZVtx", "Z Vertice" ,40,-20.,20.);
  H_ZVtxError_ = new TH1D("H_ZVtxError", "Error Z Vtx" ,200,0.,0.2);

  H_DzAfterProbEval_ = new TH1D("H_DzAfterProbEval", "#DeltaZ jet - vertex after ProbEval" ,800,-4.,4.);
  
  H_ZBeforevsAfterRJAlg_ = new TH2D("H_ZBeforevsAfterRJAlg","Z prima vs dopo rejezione tracce" ,400,-20,20,400,-20,20);

  H_ZPvVsZBeforeAlg_= new TH2D("H_ZPvVsZBeforeAlg","Z pv vs prima rejezione tracce" ,400,-20,20,400,-20,20);
  H_ZPvVsZRJAlg_= new TH2D("H_ZPvVsZRJAlg","Z pv vs dopo rejezione tracce" ,400,-20,20,400,-20,20);
  H_ZPvVsZProbEval_= new TH2D("H_ZPvVsZProbEval","Z pv vs dopo ProbEval" ,400,-20,20,400,-20,20);

  H_ZPvVsZRecoBeforeAlg_= new TH2D("H_ZPvVsZRecoBeforeAlg","Z pv vs z Reco prima rejezione tracce" ,400,-20,20,400,-20,20);
  H_ZPvVsZRecoRJAlg_= new TH2D("H_ZPvVsZRecoRJAlg","Z pv vs z Reco dopo rejezione tracce" ,400,-20,20,400,-20,20);
  H_ZPvVsZRecoProbEval_= new TH2D("H_ZPvVsZRecoProbEval","Z pv vs z Reco dopo ProbEval" ,400,-20,20,400,-20,20);

  H_ZjetSimVsReco_= new TH2D("H_ZjetSimVsReco","Z jet sim tk  vs z jet reco tk BeforeAlg",400,-20,20,400,-20,20);
  H_dzVtxSimRecoBeforeAlg_ = new TH1D("H_dzVtxSimRecoBeforeAlg","#DeltaZ  sim - reco vertex associated to the jet BeforeAlg",800,-4.,4.);
  H_dzVtxSimRecoBeforeAlgVsNtk_= new TH2D("H_dzVtxSimRecoBeforeAlgVsNtk","dZ sim vtx assigned to the jet - reco vtx  BeforeAlg vs # trk",40,-20,20,400,-20,20);
  H_ZjetSimVsReco_RJAlg_= new TH2D("H_ZjetSimVsReco_RJAlg","Z jet sim tk  vs z jet reco tk RJAlg",400,-20,20,400,-20,20);
  H_dzVtxSimRecoRJAlg_ = new TH1D("H_dzVtxSimRecoRJAlg","#DeltaZ  sim - reco vertex associated to the jet RjAlg",800,-4.,4.);
  H_dzVtxSimRecoRJAlgVsNtk_= new TH2D("H_dzVtxSimRecoRJAlgVsNtk","#DeltaZ sim vtx assigned to the jet - reco vtx  RJAlg vs # trk",40,-20,20,400,-20,20);
  H_dzVtxSimRecoProbEval_ = new TH1D("H_dzVtxSimRecoProbEval","#DeltaZ  sim - reco vertex associated to the jet ProbEval",800,-4.,4.);

  H_dzVtxHMSimRecoProbEval_ = new TH1D("H_dzVtxHMSimRecoProbEval","#DeltaZ  sim - reco vertex HM associated to the jet ProbEval",800,-4.,4.);
  H_dzVtxHMSimRecoBeforeAlg_= new TH1D("H_dzVtxHMSimRecoBeforeAlg","#DeltaZ  sim - reco vertex HM associated to the jet BeforeAlg",800,-4.,4.);
  H_dzVtxHMSimRecoRJAlg_ = new TH1D("H_dzVtxHMSimRecoRJAlg","#DeltaZ  sim - reco vertex HM associated to the jet RJAlg",800,-4.,4.);

  H_nJSimVsnJRecoVtxHMProbEval_= new TH2D("H_nJSimVsnJRecoVtxHMProbEval","# jet vtx HM sim vs reco ProbEval",30,0.,30.,30,0.,30.);
  H_nJSimVsnJRecoVtxHMBeforeAlg_= new TH2D("H_nJSimVsnJRecoVtxHMBeforeAlg","# jet vtx HM sim vs reco BeforeAlg", 30,0.,30.,30,0.,30.);
  H_nJSimVsnJRecoVtxHMRJAlg_= new TH2D("H_nJSimVsnJRecoVtxHMRJAlg","# jet vtx HM sim vs reco RJAlg",30,0.,30.,30,0.,30.);

  H_nJSimVsnJRecoVtxHMProbEvalWeighted_= new TH2D("H_nJSimVsnJRecoVtxHMProbEvalWeighted","# jet vtx HM sim vs reco ProbEvalWeighted",30,0.,30.,30,0.,30.);
  H_nJSimVsnJRecoVtxHMBeforeAlgWeighted_= new TH2D("H_nJSimVsnJRecoVtxHMBeforeAlgWeighted","# jet vtx HM sim vs reco BeforeAlgWeighted",30,0.,30.,30,0.,30.);
  H_nJSimVsnJRecoVtxHMRJAlgWeighted_= new TH2D("H_nJSimVsnJRecoVtxHMRJAlgWeighted","# jet vtx HM sim vs reco RJAlgWeighted",30,0.,30.,30,0.,30.);

  //--------------------------------------------------------------------------
  // Lettura da file .asc
  //--------------------------------------------------------------------------
  //matrice probabilita' singole
  //read from file
  //  ifstream b("MatrixProducerSing50.asc");
  ifstream b("MatrixProducerSing25.asc");
  string line1;
  while (b) {
    getline(b,line1);
    stringstream events;
    double number = 0.;
    if (line1 != "") {
      events << line1;
      events >> number;
      
      //      cout<<number<<endl;
      
      pMatrixSing25_.push_back(number);
      
    }
  }
  
  
  //matrice probabilita' combinate
  //read from file
  //  ifstream a("MatrixProducer50.asc");
  ifstream a("MatrixProducer25.asc"); 
  string line;
  while (a) {
    getline(a,line);
    stringstream events;
    double number = 0.;
    if (line != "") {
      events << line;
      events >> number;

      //      cout<<number<<endl;

      pMatrix25_.push_back(number);
      
    }
  }


  //=====ROBERTO================

  // Definitions for b-tagging
  // -------------------------
  loose_  = 2.3; 
  // high eff -> 70.49% b / 32.33% c / 8.64% uds / 10.43% g / 9.98% udsg // P.Schilling 23/10/07
  medium_ = 5.3; 
  // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
  tight_  = 4.8; 
  // high pur -> 31.94% b /  2.93% c / 0.10% uds /  0.11% g / 0.10% udsg // P.Schilling 23/10/07  

  // End of initializations
  // ----------------------
  std::cout << "Done with constructor" <<std::endl;

}


AlgValidator::~AlgValidator()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  // Write out the file
  // ------------------
  OutputFile->Write();

}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void AlgValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  //  double jetEtCut = 50.;
  // cut on Et jet
  double jetEtCut = 25.;
  // cut on sigma WAverage z jet
  int sigmaCut = 13;

  //Sigma medio vertici, calcolato come "deconvoluzione" dell'errore sul Dz Jet - Vertice dopo l'algoritmo di reiezione
  //il risultato del fit (-.05 , .05)  
  //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
  //    1  Constant     3.16510e+03   3.57364e+01   3.25501e-01   4.16937e-05
  //    2  Mean        -8.59255e-05   1.80908e-04   2.42854e-06  -2.99401e+00
  //    3  Sigma        2.19492e-02   2.10656e-04   2.76878e-05   4.84438e-01
  
  //la Sigma z media pesata del jet e' MEAN: 0.009336 RMS: 0.01597
  //per cui SigmaVtx: sqrt(pow(2.19492e-02,2) - pow(0.008,2)) = 0.02

  double SigmaVtx = 0.02;

  //======ROBERTO=======================
  //-------------------------------
  // Simulated vertices collection
  //-------------------------------

  // Take the SimVertex collection
  Handle<SimVertexContainer> SimVtx;
  std::vector<SimVertex> simVertexes;
  try{
    iEvent.getByLabel( simVtxLabel_, SimVtx );
  }
  catch(...){
    std::cerr << "Could not find the SimVertex collection" << endl;
  }
  simVertexes = *(SimVtx.product());
  
  //------------------------------
  // Pile-up vertices collection
  //------------------------------
  
  edm::Handle<HepMCProduct> theFamosPileUp;
  const HepMC::GenEvent* myGenPUEvent = 0;
  try{
    iEvent.getByLabel( simPUVtxLabel_ ,theFamosPileUp);
  } 
  catch (...) {
    std::cerr << "Could not find the pu collection" << std::endl;
    return;
  } 
  //--------------------------
  // Reco vertex collection
  //--------------------------

  Handle<reco::VertexCollection> recVtxs;
  std::vector<reco::Vertex> recoVertexes;
  try {
    iEvent.getByLabel( vtxSample_ , recVtxs);
  } 
  catch (...) {
    std::cerr << "Could not find the recovertex collection" << std::endl;
    return;
  } 
  recoVertexes = *(recVtxs.product());

  //------------------------------  
  // Simulated tracks collection 
  //------------------------------
  //simtrack per l'associazione tracce simulate vertici simulati 
  // Take the SimTracks collection
  
  Handle<SimTrackContainer> SimTk;
  std::vector<SimTrack> theSimTracks;
  try {
    iEvent.getByLabel( simTkLabel_, SimTk );
    //  cout << "SimTk size = " << SimTk->size() << endl;
  }
  catch (...) {
    std::cerr << "Could not find the SimTracks collection" << std::endl;
    return;
  } 
  theSimTracks = *(SimTk.product());

  //--------------------------
  // Reco tracks collection
  //--------------------------

  //Tracks Impact Parameter collection
  edm::Handle<reco::TrackIPTagInfoCollection> iPtagInfos;
  
  try {
    iEvent.getByLabel( impactParameterTagInfos, iPtagInfos );
  }
  catch (...) {
    std::cerr << "Could not find the b-tagging collection" << std::endl;
    return;
  } 

  //----------------------------
  //Offline calojets collection
  //----------------------------
  edm::Handle<OfflineJetCollection> caloJets;
  
  try {
    iEvent.getByLabel( offlineJetLabel_, caloJets );
  } 
  catch (...) {
    std::cerr << "Could not find the calojet collection" << std::endl;
    return;
  } 
  
  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }
  
  ////////////////////////////////// MC analysis ////////////////////////////

  //vettore di vertici, in cui il primo e' il vertice primario
  std::vector< math::XYZTLorentzVector > vec_Vertices;
  //metto via un vettore col vertice primario dell'evento
  //--------------------
  // PV + PU vertexes
  //--------------------

  myGenPUEvent = theFamosPileUp->GetEvent();

  HepMC::GenEvent::vertex_const_iterator vPUiter;
  HepMC::GenEvent::vertex_const_iterator vPUbegin = myGenPUEvent->vertices_begin();
  HepMC::GenEvent::vertex_const_iterator vPUend = myGenPUEvent->vertices_end();
  
    
  
  vec_Vertices.push_back( math::XYZTLorentzVector(SimVtx->begin()->position().x() ,
						  SimVtx->begin()->position().y() ,
						  SimVtx->begin()->position().z() ,
						  SimVtx->begin()->position().t() ) );
  
  // Loop on all pile-up events
  for ( vPUiter=vPUbegin; vPUiter!=vPUend; ++vPUiter ) { 
    
    // The origin vertex (turn it to cm's from GenEvent mm's)
    HepMC::GenVertex* vPU = *vPUiter;    
    //riempio un vettore con i vertici di pileup dell'evento
    vec_Vertices.push_back(math::XYZTLorentzVector(vPU->position().x()/10.,
						   vPU->position().y()/10.,
						   vPU->position().z()/10.,
						   vPU->position().t()/10.) );
  }

  // if(SIM_){
  //     std::vector<math::XYZTLorentzVector>::const_iterator vec_VerticesIt = vec_Vertices.begin();
  //     for(; vec_VerticesIt != vec_Vertices.end(); ++vec_VerticesIt){
  //       H_ZVtx_->Fill(vec_VerticesIt->z());
  //     }
    
  //--------------------------------
  //Minimum Dz Vertices
  //--------------------------------
  MinimumDzVtx<math::XYZTLorentzVector> simVtxMinDz(vec_Vertices);
  std::vector<double> vec_VtxMinDzSim(simVtxMinDz.evaluate());
    
  //     std::vector<double>::const_iterator vec_VtxMinDzSimIt =  vec_VtxMinDzSim.begin();
  //     for(; vec_VtxMinDzSimIt != vec_VtxMinDzSim.end(); ++vec_VtxMinDzSimIt ){
      
  //       //     H_VtxMinDz_->Fill(*vec_VtxMinDzSimIt);
      
  //     }
  // }  

  //--------------------------------------------------------- 

  std::vector<reco::Vertex>::const_iterator recoVertexesIt = recoVertexes.begin();
  for(; recoVertexesIt != recoVertexes.end(); ++recoVertexesIt){
    H_ZVtx_->Fill(recoVertexesIt->z());
    H_ZVtxError_->Fill(recoVertexesIt->zError());
  }

  //--------------------------------
  //Minimum Dz Vertices
  //--------------------------------
  MinimumDzVtx<reco::Vertex> recoVtxMinDz(recoVertexes);
  std::vector<double> vec_VtxMinDz(recoVtxMinDz.evaluate());
  
  std::vector<double>::const_iterator vec_VtxMinDzIt =  vec_VtxMinDz.begin();
  for(; vec_VtxMinDzIt != vec_VtxMinDz.end(); ++vec_VtxMinDzIt ){
      
    H_VtxMinDz_->Fill(*vec_VtxMinDzIt);
      
  }
 
 
  //--------------------------
  //Simulated tracks vertices
  //-------------------------- 
  //faccio un loop sul vettore di tracce simulate per estrarre le informazioni sulle Z della traccia
  //Come da esempio in AnalysisExamples/SimTrackerAnalysis/src/SimTrackSimVertexAnalyzer.cc
  
  //Vector of simTracks' Z, eta, phi  
  SimpleTrackCollection vec_SimTk;
  double simVertexes_size = 0.;
  simVertexes_size =  simVertexes.size();

  for (int isimvtx = 0; isimvtx < simVertexes_size;isimvtx++){  //, simvtx_it++){  
    for (std::vector<SimTrack>::iterator isimtk = theSimTracks.begin();
	 isimtk != theSimTracks.end(); ++isimtk){
      double Z_simTk = -1000.;
      double Phi_simTk = 0.;
      double Eta_simTk = 0.;
      double Pt_simTk = 0.;
      double zError_simTk = 0.000000001; //RMS = 0.0001
      //la classe simTrack non ha un metodo che restituisce eta e phi,
      //prendo eta e phi del momento, che e' un HepLorentzVector con 
      //questi metodi, e dovrebbe restituire una direzione della particella 
      //che ha generato la traccia,quindi una direzione per la traccia
      //poiche' sono tracce simulate, ricostruite con gli hits veri 
      
      Eta_simTk = (*isimtk).momentum().eta(); 
      Phi_simTk = (*isimtk).momentum().phi();

     
      if (isimtk->vertIndex() == isimvtx) {
	
	Z_simTk = simVertexes[isimvtx].position().z();

	//non considero le tracce senza vertice associato
	//L'errore sulle tracce lo pongo piccolo ripetto al valore dell'RMS del DeltaZ tra 
	//tutte le tracce simulate selezionate associate a un jet, per selezionare solo le 
	//tracce nel jet che provengono esattamente da uno dei vertici simulati
	vec_SimTk.push_back(SimpleTrack(Pt_simTk, Eta_simTk, Phi_simTk, Z_simTk, zError_simTk, 0.));

      }//end if track associated to this vertex
    }//end loop on sim tracks
  }//end loop on sim vertexes
    
  //--------------------------------
  //Reco tracks and Offline calojet (classe per simple calojet +z error e significanza traccie)
  //--------------------------------  
  std::map<int,int> map_VminSim; //map of id best Z vertex in association with dz minimum criterium 
  std::map<int,int> map_Vmin; //map of id best Z vertex in association with dz minimum criterium 
  std::map<int,int> map_idSimRJAlg; //map of id best Z vertex in association with dz minimum criterium after tracks rejection
  std::map<int,int> map_idRJAlg; //map of id best Z vertex in association with dz minimum criterium after tracks rejection

  //Et, eta, phi Offline caloJet
  double Et_Jet = 0.;
  double eta_Jet = 0.;
  double phi_Jet = 0.;

  SimpleCaloJetCollection vec_CaloJetSim; //Collection of simpleCalojet passing al cuts
  SimpleCaloJetCollection vec_CaloJet; //Collection of simpleCalojet passing al cuts

  int Njet = 0 ; //# jet
  int NjetBeforeAlg = 0 ; //# jet before cuts
  int NJet_withTracks = 0; //# jet with assigned tracks
  int NJet_withNoTracks = 0; //# jet without assigned tracks
  int Ntk = 0; //# tracks

  //loop on tracks (IptagInfo) and caloJets 
  anaobj::OfflineJetCollection::const_iterator caloJets_it = caloJets->begin();
  reco::TrackIPTagInfoCollection::const_iterator TkIpTagInfo_it = iPtagInfos->begin();
  for( ; TkIpTagInfo_it != iPtagInfos->end(); ++TkIpTagInfo_it ) {
    
    //vector of reco selected tracks' Pt, eta, phi, Z
    SimpleTrackCollection vec_RecoTk;

    //vector of sim selected tracks' Pt, eta, phi, Z
    SimpleTrackCollection vec_SimTk_selected;


    //Et, eta, phi Offline caloJet
    Et_Jet = caloJets_it->et();
    eta_Jet = caloJets_it->eta();
    phi_Jet =  caloJets_it->phi();

    // Take the selectedTracks vector
    // Return the vector of tracks for which the IP information is available Quality cuts are applied to reject fake tracks.
    const reco::TrackRefVector & vec_TkColl = TkIpTagInfo_it->selectedTracks();
    
    // Take the IP vector (ordered as the selectedTracks vector)
    const vector<reco::TrackIPTagInfo::TrackIPData> & vec_TkIP = TkIpTagInfo_it->impactParameterData();
    
    //loop on selected tracks collection and tracks IP  
    vector<reco::TrackIPTagInfo::TrackIPData>::const_iterator TkIP_it = vec_TkIP.begin();
    RefVector<reco::TrackCollection>::const_iterator TkColl_it = vec_TkColl.begin();
    for ( ; TkColl_it != vec_TkColl.end(); ++TkColl_it, ++TkIP_it ) {
      
      //2D Significance
      double Tk_IpS2D = 0.;
      Tk_IpS2D =TkIP_it->ip2d.significance();
      //      cout<<"significanza 2D IP traccia: "<<Tk_IpS2D<<endl;
 
      //Reco Track Z
      double Tk_dz = 0.;
      Tk_dz = (*TkColl_it)->dz();
           
      //Reco Track Z error
      double Tk_dz_error = 0.;
      Tk_dz_error= (*TkColl_it)->dzError();
         
      //Reco Track eta
      double Tk_eta = 0.;
      Tk_eta = (*TkColl_it)->eta();
     
      //Reco Track phi
      double Tk_phi = 0.;
      Tk_phi = (*TkColl_it)->phi(); 
    
      //Reco Track Pt
      double Tk_Pt = 0.;
      Tk_Pt = (*TkColl_it)->pt();
      // cout<<"Pt traccia: "<<Tk_Pt<<endl;
    
      //cuts in filling recoTrackCollection
      //-----------------------------------

      if(Et_Jet>jetEtCut && fabs(Tk_IpS2D) < 3.  && eta_Jet < 3 ){
	vec_RecoTk.push_back(SimpleTrack(Tk_Pt,
					 Tk_eta,
					 Tk_phi,
					 Tk_dz, 
					 Tk_dz_error,
					 Tk_IpS2D)
			     );
	
      }//end cut Et(jet)>25, IP significance < 3. , eta (jet) < 3.   
    }// end loop on selected tracks
  
    //----------------------------------
    //Taglio Et(jet)>25 , eta(jet) < 3
    //----------------------------------
    
    if(Et_Jet>jetEtCut && eta_Jet < 3){

      //----------------------------------------------
      //Simtracks - reco selected tracks association 
      //Delta2 minimum simTrack - reco selected track (ok)
      //----------------------------------------------
      //loop on reco selected tracks variabiles
      anaobj::SimpleTrackCollection::const_iterator vec_RecoTk_iter = vec_RecoTk.begin();
      for(; vec_RecoTk_iter !=  vec_RecoTk.end(); vec_RecoTk_iter++ ){
	    
	//richiamo la funzione d2Associator che riceve in input una collezione di SimpleTrack, 
	//il valore di eta e phi della traccia per cui cerchiamo l'associazione e restituisce 
	//un pair<pair<int ID, double VALUE>, pair<int ID, double VALUE> >, in cui la prima e' 
	//l'associazione migliore e la seconda la second best
	//------------------------------------------------------------------------------------
	    
	pair<pair<int, double>, pair<int, double> > d2SimReco( d2Associator( vec_SimTk, vec_RecoTk_iter->eta(), vec_RecoTk_iter->phi() ) );
	int idSimTk = d2SimReco.first.first ;
	    
	//fill vector collection of selected simulated tracks' Z 
	vec_SimTk_selected.push_back( vec_SimTk[idSimTk ] );
      }//end loop reco selected tracks 

      //----------------------------------------------------
      // Jet Z: Wieghted Average on tracks' Z, with tracks
      //        Z error
      //----------------------------------------------------
      
      //richiamo la funzione wAverager per calcolare la media pesata delle tracce ricostruite
      //richiede in input una simpletrackcollection e restituisce una coppia pair(contatore # 
      //tracce, pair (media, errore))
      //----------------------------------------------------------------------------------------
      pair< int, pair< double, double > > wAvgSimTk;
      pair< int, pair< double, double > > wAvgRecoTk;
     
      //per le tracce simulate l'errore sulla zl'ho posto 0.000000001
      wAvgSimTk = wAverager(vec_SimTk_selected);
      wAvgRecoTk = wAverager(vec_RecoTk);
	
      double wAvgSimTkValue = wAvgSimTk.second.first; //coordinata z media delle tracce
      //      double  wAvgSimTkError = wAvgSimTk.second.second; //errore sulla media pesata
      int NTks_S3_Sim = wAvgSimTk.first ; //contatore per il numero di tracce con significanza minore di 3

      double wAvgRecoTkValue = wAvgRecoTk.second.first; //coordinata z media delle tracce
      double  wAvgRecoTkError = wAvgRecoTk.second.second; //errore sulla media pesata
      int NTks_S3 = wAvgRecoTk.first ; //contatore per il numero di tracce con significanza minore di 3
    
      if (NTks_S3 != 0 && NTks_S3_Sim != 0 ){ 
	//Quanti vertici di jet con la stessa Z tra sim e reco	
	H_ZjetSimVsReco_->Fill( wAvgSimTkValue , wAvgRecoTkValue );
      }
      
      //if # of selected tracks with IP Significance < 3. != 0
      if (NTks_S3 != 0){ 	
	  
	//Histogram Z Jet
	H_ZJetWAvg_->Fill(wAvgRecoTkValue);
	H_ZJetWAvgError_->Fill(wAvgRecoTkError);

	NJet_withTracks++;
      } 
      else {
	NJet_withNoTracks++;
      }

      //---------------------------------
      // Delta Z minimo tra jet e vertici
      // con Z media pesata tracce
      //---------------------------------
      //if # of selected tracks with IP Significance < 3. != 0
      if ( NTks_S3 != 0 && NTks_S3_Sim != 0 ){
	
	//Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata e i vertici (PV e PU)	
	//la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
	//la funzione vuole in input un XYZTLorentz vector e un double
	//----------------------------------------------------------------------------------------------------------------------
	pair<pair<int, double>, pair<int, double> > dzSimFirstSecond;
	pair<pair<int, double>, pair<int, double> > dzFirstSecond;
		
	dzSimFirstSecond =  dzAssociator(vec_Vertices, wAvgSimTkValue) ;
	  
	dzFirstSecond = dzAssociator(recoVertexes, wAvgRecoTkValue) ;

	pair<int, double> * firstSim = &(dzSimFirstSecond.first);
	//	pair<int, double> * secondSim = &(dzSimFirstSecond.second);
	int VminSim = firstSim->first;
	//	double Dz_Jet_VtxSim = firstSim->second;
	//	double Dz_Jet_Vtx_2BestSim = secondSim->second;

	pair<int, double> * first = &(dzFirstSecond.first);
	pair<int, double> * second = &(dzFirstSecond.second);
	int Vmin = first->first;
	double Dz_Jet_Vtx = first->second;
	double Dz_Jet_Vtx_2Best = second->second;
	
	//Fill map with idVertex associated with minimum dz criterium with WAvg z jet
 	map_VminSim.insert(make_pair( NjetBeforeAlg, VminSim) );
	map_Vmin.insert(make_pair( NjetBeforeAlg, Vmin) );

	//distribuzione del Delta Z tra vertice simulato e ericostruito cui assegno il jet
	//se e' inferiore alla risoluzione dell'algoritmo di ricostruzione dei vertici ricostruiti
	//meglio non posso fare
	double dzVtxSimRecoBeforeAlg = vec_Vertices[VminSim].z() - recoVertexes[Vmin].z();

	H_dzVtxSimRecoBeforeAlg_->Fill(dzVtxSimRecoBeforeAlg);

	H_dzVtxSimRecoBeforeAlgVsNtk_->Fill( NTks_S3 ,dzVtxSimRecoBeforeAlg );
	
	NjetBeforeAlg++;

	//Histograms minimum Dz Jet - Vertices	
	H_Vmin_->Fill(Vmin);
	H_DzJetVtxWAvg_->Fill(Dz_Jet_Vtx);
	H_DzJetVtxWAvg2Best_->Fill(Dz_Jet_Vtx_2Best);
	H_DzJetVtxWAvg_vs_2Best_->Fill(fabs(Dz_Jet_Vtx),fabs(Dz_Jet_Vtx_2Best));
      
      }	//end if # of selected tracks with IP Significance < 3. != 0
	
	//----------------------------------------------
	//Algoritmo di reiezione tracce ricostruite nel jet
	//Z jet calcolata con la media pesata, scartando 
	//le tracce piu' distanti di 5 sigma dal valore 
	//della media pesata con algoritmo iterativo
	//----------------------------------------------
	//chiamo l'algoritmo di reiezione (e' una classe con il metodo eval(# sigma taglio)) sulla collezione di tracce ricostruite
	//agisce se il numero di trace nel jet e' > 2
      SimpleTrackCollection simTkClosest;
      SimpleTrackCollection recoTkClosest;
      
      RejectionAlg simTkRejection(vec_SimTk_selected);
      simTkClosest = (simTkRejection.eval(sigmaCut));
      
      
      RejectionAlg tkRejection(vec_RecoTk);
      recoTkClosest = (tkRejection.eval(sigmaCut));

      //chiamo la funzione che fa la media pesata wAvenger, che restituisce una pair(# tracce, pair(value, error))
      pair< int, pair< double, double > > zJetRejectionAlgSim( wAverager(simTkClosest));
      //z jet e suo errore
      double zJetSim = zJetRejectionAlgSim.second.first;
      double zJetSimError = zJetRejectionAlgSim.second.second;
      int NTksSimRJAlg_S3 = zJetRejectionAlgSim.first ;
      
      //chiamo la funzione che fa la media pesata wAvenger, che restituisce una pair(# tracce, pair(value, error))
      pair< int, pair< double, double > > zJetRejectionAlg( wAverager(recoTkClosest));
      //z jet e suo errore
      double zJet = zJetRejectionAlg.second.first;
      double zJetError = zJetRejectionAlg.second.second;
      int NTksRJAlg_S3 = zJetRejectionAlg.first ; //contatore per il numero di tracce con significanza minore di 3	
      //sono nello stesso numero di dei jet prima della reiezione, perche' quelli con 2 o meno tracce skippano l'algoritmo
      //e gli altri hanno tracce comunque, cambiera' il contenuto in tracce dei jet


      //quanti jet con la stessa z vertice tra sim e reco dopo la rejezione delle tracce 
      if (NTksRJAlg_S3 != 0 && NTksSimRJAlg_S3 != 0 ){ 
	//Quanti vertici di jet con la stessa Z tra sim e reco	
	H_ZjetSimVsReco_RJAlg_->Fill( zJetSim  , zJet );
      }

      if(NTksRJAlg_S3 != 0){
	//riempio gli istogrammi di value e error (Z del jet calcolata e suo errore) 
	//Histograms of Z jet and error after the rejection algorithm (5 sigma)
	H_ZJet_RJAlg_->Fill(zJet);
	H_ZErrorJet_RJAlg_->Fill(zJetError);
      }

      //---------------------------------
      // Delta Z minimo tra jet e vertici
      // con Z media pesata tracce
      //---------------------------------
	
      //Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata dopo la reiezione e i vertici (PV e PU)	
      //la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
      //la funzione vuole in input un XYZTLorentz vector e un double
      //----------------------------------------------------------------------------------------------------------------------
      //if # of selected tracks with IP Significance < 3. != 0
      if (NTksRJAlg_S3 != 0){	
	pair<pair<int, double>, pair<int, double> > dzAfterRejectionSim;
	pair<pair<int, double>, pair<int, double> > dzAfterRejection;

	dzAfterRejectionSim = dzAssociator(vec_Vertices,zJetSim ) ;

	dzAfterRejection = dzAssociator(recoVertexes,zJet ) ;

	int  idSimRJAlg = dzAfterRejectionSim.first.first;
	int  id2SimRJAlg = dzAfterRejectionSim.second.first;
	double dzSimRJAlg = dzAfterRejectionSim.first.second;
	double dz2SimRJAlg = dzAfterRejectionSim.second.second; //second best

	int  idRJAlg = dzAfterRejection.first.first;
	int  id2RJAlg = dzAfterRejection.second.first;
	double dzRJAlg = dzAfterRejection.first.second;
	double dz2RJAlg = dzAfterRejection.second.second; //second best

	
	//Histograms of id vertex assigned to a jet after rejection algorithm
	H_idRJAlg_->Fill(idRJAlg);
	//Histogram of minimum dz jet - vertices after rejection (5 sigma)
	H_DzJetVtx_RJAlg_->Fill(dzRJAlg);
	H_Dz2JetVtx_RJAlg_->Fill(dz2RJAlg);
	H_DzJetVtx_RJAlg_vs_2Best_->Fill(dzRJAlg, dz2RJAlg);

	//distribuzione del Delta Z tra vertice simulato e ericostruito cui assegno il jet
	//se e' inferiore alla risoluzione dell'algoritmo di ricostruzione dei vertici ricostruiti
	//meglio non posso fare
	
	double dzVtxSimRecoRJAlg = vec_Vertices[idSimRJAlg].z() - recoVertexes[idRJAlg].z();

	H_dzVtxSimRecoRJAlg_->Fill(dzVtxSimRecoRJAlg);

	H_dzVtxSimRecoRJAlgVsNtk_->Fill( NTksRJAlg_S3 ,dzVtxSimRecoRJAlg );


	//------------------
	// Calcolo P1/P2 Sim
	//------------------
	//Uso invece dell'errore sui vertici dato dalla collezione, l'errore sulla DeltaZ jet vertice dopo l'algoritmo di 
	//reiezione calcolato per il campione QCD170-230 = 0.03556

	//Ho ricavato l'errore sul vertice dalla "deconvoluzione" con l'errore medio sui jet, 
	//ora lo combino con l'errore per ogni jet per ottenere l'errore sulla Dz del jet

	double dzJetSimError = sqrt(pow(zJetSimError,2)+ pow(SigmaVtx ,2));
	double dz2JetSimError = sqrt(pow(zJetSimError,2)+ pow(SigmaVtx,2));

	double p1Sim = TMath::Erfc(fabs(dzSimRJAlg/dzJetSimError));
	double p2Sim = TMath::Erfc(fabs(dz2SimRJAlg/dz2JetSimError));
	
	if(p1Sim < 10e-30) p1Sim = 10e-30;
	if(p2Sim < 10e-30) p2Sim = 10e-30;
	
	map_idSimRJAlg.insert(make_pair( Njet, idSimRJAlg) );
	
	vec_CaloJetSim.push_back(SimpleCaloJet( Et_Jet,
						eta_Jet,
						phi_Jet,
						zJetSim,
						zJetSimError,
						Njet,
						idSimRJAlg,
						id2SimRJAlg,
						dzSimRJAlg,
						dz2SimRJAlg,
						p1Sim,
						p2Sim ) );

	//------------------
	// Calcolo P1/P2
	//------------------
	//Uso invece dell'errore sui vertici dato dalla collezione, l'errore sulla DeltaZ jet vertice dopo l'algoritmo di 
	//reiezione calcolato per il campione QCD170-230 = 0.03556

	//Ho ricavato l'errore sul vertice dalla "deconvoluzione" con l'errore medio sui jet, 
	//ora lo combino con l'errore per ogni jet per ottenere l'errore sulla Dz del jet

	double dzJetError = sqrt(pow(zJetError,2)+ pow(SigmaVtx ,2));
	double dz2JetError = sqrt(pow(zJetError,2)+ pow(SigmaVtx,2));

	double p1 = TMath::Erfc(fabs(dzRJAlg/dzJetError));
	double p2 = TMath::Erfc(fabs(dz2RJAlg/dz2JetError));
	
	if(p1 < 10e-30) p1 = 10e-30;
	if(p2 < 10e-30) p2 = 10e-30;

	
	map_idRJAlg.insert(make_pair( Njet, idRJAlg) );
	
	vec_CaloJet.push_back(SimpleCaloJet( Et_Jet,
					     eta_Jet,
					     phi_Jet,
					     zJet,
					     zJetError,
					     Njet,
					     idRJAlg,
					     id2RJAlg,
					     dzRJAlg,
					     dz2RJAlg,
					     p1,
					     p2 ) );
	

	//Histograms of fraction between probability of best vertex association and second best 
	H_P1P2_->Fill(p1/p2);
	//Histograms of # Tracks for jet with #DeltaZ best < 5#sigma
	H_nTrack_RJAlg_->Fill(NTksRJAlg_S3);
	//Histograms of # Tracks for jet with #DeltaZ best < 5#sigma vs z jet (error)
	H_nTkRJAlgVSzJetError_->Fill(NTksRJAlg_S3,zJetError);
	H_nTkRJAlgVSzJet_->Fill(NTksRJAlg_S3,zJet);
	
	Ntk += NTksRJAlg_S3;
	Njet++;
	

      }//end cut NTksRJAlg_S3 != 0
    }// end cut Et(jet>25) eta<3
    
    ++caloJets_it;
  }//end loop calojet

  cout<<"# jet= "<<Njet<<endl;
  cout<<"# tracks= "<<Ntk<<endl;
  
  //Histograms # of jet with (without) tracks in the sample  
  H_NJet_withTracks_->Fill(NJet_withTracks); //# jet nell'evento con tracce assegnate
  H_NJet_withNoTracks_->Fill(NJet_withNoTracks); //# jet nell'evento che non ha tracce assegnate
  H_nJet_->Fill(NJet_withTracks + NJet_withNoTracks); //# jet nell'evento 
    
  //if jetCollection after cuts !empty 
  if( !vec_CaloJet.empty() &&  !vec_CaloJetSim.empty()){

    //--------------------------------------
    // #jet before any algorithm 
    //--------------------------------------
    //------------------------------------------------------------------------------
    // # vertex with # of jet specified by key
    //------------------------------------------------------------------------------
    // VerticesJetCounter.countJet: count # jet to every vertex id specified by key
    // VerticesJetCounter.countVtx: count # vertex with # of jet specified by key
    //------------------------------------------------------------------------------ 
    VerticesJetCounter JetSimCountBeforeAlg;
    std::map<int,int> nJetSimMapBeforeAlg(JetSimCountBeforeAlg.countJet(map_VminSim));
    std::map<int,int> nVtxSimMapBeforeAlg(JetSimCountBeforeAlg.countVtx(nJetSimMapBeforeAlg));

    VerticesJetCounter JetCountBeforeAlg;
    std::map<int,int> nJetMapBeforeAlg(JetCountBeforeAlg.countJet(map_Vmin));
    std::map<int,int> nVtxMapBeforeAlg(JetCountBeforeAlg.countVtx(nJetMapBeforeAlg));

    if( !nVtxMapBeforeAlg.empty() )  H_nJMax_BeforeAlg_->Fill((*nVtxMapBeforeAlg.rbegin()).first + NJet_withNoTracks);

   
    //----------------------
    //PROBEVAL Sim
    //----------------------

    ProbEval ValueSimProb_;
    //mappa di associazione (idJet, idVtx)
    map<int,int> mapSimProbEval(ValueSimProb_.evalProb(vec_CaloJetSim));

    VerticesJetCounter JetSimCountProbEval;
    //Overloaded method
    std::map<int,int> nJetSimMapProbEval(JetSimCountProbEval.countJet(mapSimProbEval));
    //restituisce una mappa con # jet, # vertici con quel # jet
    std::map<int,int> nVtxSimMapProbEval(JetSimCountProbEval.countVtx(nJetSimMapProbEval));

    std::map<int,int>::const_iterator mapSimProbEvalIt = mapSimProbEval.begin();
    for(; mapSimProbEvalIt != mapSimProbEval.end() ; ++mapSimProbEvalIt ){
    
      //      int idJetSim =(*mapSimProbEvalIt).first;
      int idVtxSim =(*mapSimProbEvalIt).second;

      double zVtxSim;

      zVtxSim = vec_Vertices[idVtxSim].z();

      //      double zJetSim = vec_CaloJetSim[idJetSim].z();

    }

    //----------------------
    //PROBEVAL
    //----------------------

    ProbEval ValueProb_;
    //mappa di associazione (idJet, idVtx)
    map<int,int> mapProbEval(ValueProb_.evalProb(vec_CaloJet));

    VerticesJetCounter JetCountProbEval;
    //Overloaded method
    std::map<int,int> nJetMapProbEval(JetCountProbEval.countJet(mapProbEval));
    //restituisce una mappa con # jet, # vertici con quel # jet
    std::map<int,int> nVtxMapProbEval(JetCountProbEval.countVtx(nJetMapProbEval));
      
    std::map<int,int>::const_iterator nVtxMapProbEvalIt = nVtxMapProbEval.begin();
    for(; nVtxMapProbEvalIt!=  nVtxMapProbEval.end(); ++nVtxMapProbEvalIt ){
      //Histograms of # of vertexes with a specified jet multeplicity
      H_nVtxMapProbEval_->Fill(nVtxMapProbEvalIt->first,(double)nVtxMapProbEvalIt->second);
    }
    
  

    std::map<int,int>::const_iterator mapProbEvalIt = mapProbEval.begin();
    for(; mapProbEvalIt != mapProbEval.end() ; ++mapProbEvalIt ){
      int idJet =(*mapProbEvalIt).first;
      int idVtx =(*mapProbEvalIt).second;

      double zVtx;
      
      zVtx = recoVertexes[idVtx].z();

      double zJet = vec_CaloJet[idJet].z();
      
      //Histograms of Dz jet - vertex after SAAlg
      H_DzAfterProbEval_->Fill(zVtx - zJet);
      
      //distribuzione del Delta Z tra vertice simulato e ricostruito cui assegno il jet
      //se e' inferiore alla risoluzione dell'algoritmo di ricostruzione dei vertici ricostruiti
      //meglio non posso fare
      int idSimProbEval = mapSimProbEval[idJet];
      double dzVtxSimRecoProbEval = vec_Vertices[idSimProbEval].z() - recoVertexes[idVtx].z();
	
      H_dzVtxSimRecoProbEval_->Fill(dzVtxSimRecoProbEval);
    }

    if(!nVtxMapProbEval.empty())  H_nJMax_ProbEval_->Fill((*nVtxMapProbEval.rbegin()).first + NJet_withNoTracks);

  
    //--------------------------------------
    //RJAlg njet and prob event status
    //--------------------------------------
    
    //------------------------------------------------------------------------------
    // # vertex with # of jet specified by key (RJAlg)
    //------------------------------------------------------------------------------
    // VerticesJetCounter.countJet: count # jet to every vertex id specified by key
    // VerticesJetCounter.countVtx: count # vertex with # of jet specified by key
    //------------------------------------------------------------------------------
    VerticesJetCounter JetSimCountRJAlg;
    //Overloaded method
    std::map<int,int> nJetSimMapRJAlg(JetSimCountRJAlg.countJet(map_idSimRJAlg));
    //restituisce una mappa con # jet, # vertici con quel # jet
    std::map<int,int> nVtxSimMapRJAlg(JetSimCountRJAlg.countVtx(nJetSimMapRJAlg));

    VerticesJetCounter JetCountRJAlg;
    //Overloaded method
    std::map<int,int> nJetMapRJAlg(JetCountRJAlg.countJet(map_idRJAlg));
    //restituisce una mappa con # jet, # vertici con quel # jet
    std::map<int,int> nVtxMapRJAlg(JetCountRJAlg.countVtx(nJetMapRJAlg));

    std::map<int,int>::const_iterator nVtxMapRJAlgIt = nVtxMapRJAlg.begin();
    for(; nVtxMapRJAlgIt!=  nVtxMapRJAlg.end(); ++nVtxMapRJAlgIt ){
      if(nVtxMapRJAlgIt->first != 0){
	//Histograms of # of vertexes with a specified jet multeplicity
	H_nVtxMapRJAlg_->Fill(nVtxMapRJAlgIt->first,(double)nVtxMapRJAlgIt->second);
      }
    }
 
    //Histograms of # jet for vertex with max multeplicity  vs # jet tot after RJAlg
    if(!nVtxMapRJAlg.empty()) H_nJMaxVsNjTot_RJAlg_->Fill( (*nVtxMapRJAlg.rbegin()).first, Njet );
    
    if(!nVtxMapRJAlg.empty())  H_nJMax_RJAlg_->Fill((*nVtxMapRJAlg.rbegin()).first + NJet_withNoTracks);

    //-----------------------------
    //Dz vertices - primary vertex 
    //-----------------------------
    //cerco il vertice con max # jet assegnati	
    std::map<int , int> VtxMaxSimProbEval; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxSimBeforeAlg; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxSimRJAlg; //mappa con key il # di jet e value l'id del vertice

    //cerco il vertice con max # jet assegnati	
    std::map<int , int> VtxMaxProbEval; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxBeforeAlg; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxRJAlg; //mappa con key il # di jet e value l'id del vertice

    //riempimento mappe con key il # di jet e value l'id del vertice
    std::map<int,int>::const_iterator nJetSimMapProbEvalIt = nJetSimMapProbEval.begin();
    for(; nJetSimMapProbEvalIt != nJetSimMapProbEval.end() ; nJetSimMapProbEvalIt++){
      VtxMaxSimProbEval.insert( make_pair( nJetSimMapProbEvalIt->second, nJetSimMapProbEvalIt->first ) );
    }

    //riempimento mappe con key il # di jet e value l'id del vertice
    std::map<int,int>::const_iterator nJetMapProbEvalIt = nJetMapProbEval.begin();
    for(; nJetMapProbEvalIt != nJetMapProbEval.end() ; nJetMapProbEvalIt++){
      VtxMaxProbEval.insert( make_pair( nJetMapProbEvalIt->second, nJetMapProbEvalIt->first ) );
    }

    std::map<int,int>::const_iterator nJetSimMapBeforeAlgIt = nJetSimMapBeforeAlg.begin();
    for(; nJetSimMapBeforeAlgIt != nJetSimMapBeforeAlg.end() ; nJetSimMapBeforeAlgIt++){
      VtxMaxSimBeforeAlg.insert( make_pair( nJetSimMapBeforeAlgIt->second, nJetSimMapBeforeAlgIt->first ) );
    }

    std::map<int,int>::const_iterator nJetMapBeforeAlgIt = nJetMapBeforeAlg.begin();
    for(; nJetMapBeforeAlgIt != nJetMapBeforeAlg.end() ; nJetMapBeforeAlgIt++){
      VtxMaxBeforeAlg.insert( make_pair( nJetMapBeforeAlgIt->second, nJetMapBeforeAlgIt->first ) );
    }

    std::map<int,int>::const_iterator nJetSimMapRJAlgIt = nJetSimMapRJAlg.begin();
    for(; nJetSimMapRJAlgIt != nJetSimMapRJAlg.end() ; nJetSimMapRJAlgIt++){
      VtxMaxSimRJAlg.insert( make_pair( nJetSimMapRJAlgIt->second, nJetSimMapRJAlgIt->first ) );
    }

    std::map<int,int>::const_iterator nJetMapRJAlgIt = nJetMapRJAlg.begin();
    for(; nJetMapRJAlgIt != nJetMapRJAlg.end() ; nJetMapRJAlgIt++){
      VtxMaxRJAlg.insert( make_pair( nJetMapRJAlgIt->second, nJetMapRJAlgIt->first ) );
    }
    //-------------------------------
    //Calcolo DeltaZ
    //-------------------------------
    //    int idVtx = 0;

    //id vertice con max num jet e sua Z (la mappa e' ordinata in valori crescenti 
    //della key che per queste mappe e' scelta esser il # di jet)

    int  idPVSimProbEval = VtxMaxSimProbEval.rbegin()->second;
    int nJSimProbEval = VtxMaxSimProbEval.rbegin()->first;
    double zPVSimProbEval;

    int  idPVProbEval = VtxMaxProbEval.rbegin()->second;
    int nJProbEval = VtxMaxProbEval.rbegin()->first;
    double zPVProbEval;

    zPVSimProbEval = vec_Vertices[idPVSimProbEval].z(); 
    zPVProbEval = recoVertexes[idPVProbEval].z(); 

    int  idPVSimBeforeAlg = VtxMaxSimBeforeAlg.rbegin()->second;
    int nJSimBeforeAlg = VtxMaxSimBeforeAlg.rbegin()->first;
    double zPVSimBeforeAlg ;
    
    int  idPVBeforeAlg = VtxMaxBeforeAlg.rbegin()->second;
    int nJBeforeAlg = VtxMaxBeforeAlg.rbegin()->first;
    double zPVBeforeAlg ;
    
    zPVSimBeforeAlg = vec_Vertices[idPVSimBeforeAlg].z(); 
    zPVBeforeAlg = recoVertexes[idPVBeforeAlg].z();
 
    int  idPVSimRJAlg = VtxMaxSimRJAlg.rbegin()->second;
    int nJSimRJAlg = VtxMaxSimRJAlg.rbegin()->first;
    double zPVSimRJAlg;

    int  idPVRJAlg = VtxMaxRJAlg.rbegin()->second;
    int nJRJAlg = VtxMaxRJAlg.rbegin()->first;
    double zPVRJAlg;

    zPVSimRJAlg = vec_Vertices[idPVSimRJAlg].z(); 
    zPVRJAlg = recoVertexes[idPVRJAlg].z();

    H_ZBeforevsAfterRJAlg_->Fill(zPVBeforeAlg,zPVRJAlg);
      
    H_ZPvVsZBeforeAlg_->Fill(vec_Vertices.begin()->z(), zPVSimBeforeAlg);
    H_ZPvVsZRJAlg_->Fill(vec_Vertices.begin()->z(),zPVSimRJAlg);
    H_ZPvVsZProbEval_->Fill(vec_Vertices.begin()->z(),zPVSimProbEval);

    H_ZPvVsZRecoBeforeAlg_->Fill(vec_Vertices.begin()->z(), zPVBeforeAlg);
    H_ZPvVsZRecoRJAlg_->Fill(vec_Vertices.begin()->z(),zPVRJAlg);
    H_ZPvVsZRecoProbEval_->Fill(vec_Vertices.begin()->z(),zPVProbEval);

    H_dzVtxHMSimRecoProbEval_->Fill(zPVSimProbEval - zPVProbEval);
    H_dzVtxHMSimRecoBeforeAlg_->Fill(zPVSimBeforeAlg - zPVBeforeAlg);
    H_dzVtxHMSimRecoRJAlg_->Fill(zPVSimRJAlg - zPVRJAlg);

    H_nJSimVsnJRecoVtxHMProbEval_->Fill(nJSimProbEval, nJProbEval);
    H_nJSimVsnJRecoVtxHMBeforeAlg_->Fill(nJSimBeforeAlg, nJBeforeAlg);
    H_nJSimVsnJRecoVtxHMRJAlg_->Fill(nJSimRJAlg, nJRJAlg);

    //---------------------------------------------

    vector< int > vec_idJetRecoHMProbEval;//vettori di id jet assegnati al vertice con massima molteplicita'
    vector< int > vec_idJetSimHMProbEval;
   
    std::map<int,int>::const_iterator mapProbEvalIt1 = mapProbEval.begin();
    for(; mapProbEvalIt1 !=  mapProbEval.end(); ++mapProbEvalIt1 ){
      // se il vertice e' quello a massima molteplicita' mi salvo le id dei jet ad esso
      // associati per confrontarli con quelli assegnati al vertice con massima 
      // molteplicita' simulato
      if(mapProbEvalIt1->second == idPVProbEval){
	vec_idJetRecoHMProbEval.push_back(mapProbEvalIt1->first);
      }
    }
    cout<<"Reco HM #Jet"<<vec_idJetRecoHMProbEval.size()<<", "<<nJProbEval<<endl;
    
    std::map<int,int>::const_iterator mapSimProbEvalIt1 = mapSimProbEval.begin();
    for(; mapSimProbEvalIt1 !=  mapSimProbEval.end(); ++mapSimProbEvalIt1 ){
      if(mapSimProbEvalIt1->second == idPVSimProbEval){
	vec_idJetSimHMProbEval.push_back(mapSimProbEvalIt1->first);
      }
    }
    cout<<"Sim HM #Jet"<<vec_idJetSimHMProbEval.size()<<", "<<nJSimProbEval<<endl;

    //confronto tra le due collezioni di jet per contare quanti sono gli stessi 
    //jet assegnati al vertice con massima molteplicita'
    //---------------------------------------------
    int sameJetCounterProbEval = 0;

    std::vector< int >::const_iterator  vec_idJetSimHMProbEvalIt = vec_idJetSimHMProbEval.begin();
    for(; vec_idJetSimHMProbEvalIt != vec_idJetSimHMProbEval.end(); vec_idJetSimHMProbEvalIt++){
      std::vector< int >::const_iterator  vec_idJetRecoHMProbEvalIt = vec_idJetRecoHMProbEval.begin();
      for(; vec_idJetRecoHMProbEvalIt != vec_idJetRecoHMProbEval.end(); vec_idJetRecoHMProbEvalIt++){
	 //se l'id del jet ricostruito e' uguale all'id del simulato incremento il contatore
	 if( *vec_idJetRecoHMProbEvalIt == *vec_idJetSimHMProbEvalIt) sameJetCounterProbEval++;
       }
    }

    double WeightProbEval= 0.;
    if(vec_idJetRecoHMProbEval.size() > vec_idJetSimHMProbEval.size()){
      WeightProbEval = (double)sameJetCounterProbEval/(double)vec_idJetRecoHMProbEval.size();
    }else{
      WeightProbEval = (double)sameJetCounterProbEval/(double)vec_idJetSimHMProbEval.size();
    }
    cout<<"sameJetCounterProbEval: "<<sameJetCounterProbEval<<" , WeightProbEval: "<<WeightProbEval<<endl;

    H_nJSimVsnJRecoVtxHMProbEvalWeighted_->Fill(nJSimProbEval, nJProbEval, WeightProbEval);
  
    //------------------------------------------------------------

    vector< int > vec_idJetRecoHMRJAlg;//vettori di id jet assegnati al vertice con massima molteplicita'
    vector< int > vec_idJetSimHMRJAlg;
   
    std::map<int,int>::const_iterator mapRJAlgIt = map_idRJAlg.begin();
    for(; mapRJAlgIt!=  map_idRJAlg.end(); ++mapRJAlgIt ){
      // se il vertice e' quello a massima molteplicita' mi salvo le id dei jet ad esso
      // associati per confrontarli con quelli assegnati al vertice con massima 
      // molteplicita' simulato
      if(mapRJAlgIt->second == idPVRJAlg){
	vec_idJetRecoHMRJAlg.push_back(mapRJAlgIt->first);
      }
    }
    cout<<"Reco HM #Jet"<<vec_idJetRecoHMRJAlg.size()<<", "<<nJRJAlg<<endl;
    
    std::map<int,int>::const_iterator mapSimRJAlgIt = map_idSimRJAlg.begin();
    for(; mapSimRJAlgIt!=  map_idSimRJAlg.end(); ++mapSimRJAlgIt ){
      if(mapSimRJAlgIt->second == idPVSimRJAlg){
	vec_idJetSimHMRJAlg.push_back(mapSimRJAlgIt->first);
      }
    }
    cout<<"Sim HM #Jet"<<vec_idJetSimHMRJAlg.size()<<", "<<nJSimRJAlg<<endl;

    //confronto tra le due collezioni di jet per contare quanti sono gli stessi 
    //jet assegnati al vertice con massima molteplicita'
    //---------------------------------------------
    int sameJetCounterRJAlg = 0;
 
    std::vector< int >::const_iterator  vec_idJetSimHMRJAlgIt = vec_idJetSimHMRJAlg.begin();
    for(; vec_idJetSimHMRJAlgIt != vec_idJetSimHMRJAlg.end(); vec_idJetSimHMRJAlgIt++){
      std::vector< int >::const_iterator  vec_idJetRecoHMRJAlgIt = vec_idJetRecoHMRJAlg.begin();
       for(; vec_idJetRecoHMRJAlgIt != vec_idJetRecoHMRJAlg.end(); vec_idJetRecoHMRJAlgIt++){
	 //se l'id del jet ricostruito e' uguale all'id del simulato incremento il contatore
	 if( *vec_idJetRecoHMRJAlgIt == *vec_idJetSimHMRJAlgIt) sameJetCounterRJAlg++;
       }
    } 
    double WeightRJAlg = 0.;
    if(vec_idJetRecoHMRJAlg.size() > vec_idJetSimHMRJAlg.size()){
      WeightRJAlg = (double)sameJetCounterRJAlg/(double)vec_idJetRecoHMRJAlg.size();
    }else{
      WeightRJAlg = (double)sameJetCounterRJAlg/(double)vec_idJetSimHMRJAlg.size();
    }

    cout<<"sameJetCounterRJAlg: "<<sameJetCounterRJAlg<<" , WeightRJAlg: "<<WeightRJAlg<<endl;
    
    H_nJSimVsnJRecoVtxHMRJAlgWeighted_->Fill(nJSimRJAlg, nJRJAlg, WeightRJAlg);

 //------------------------------------------------------------

    vector< int > vec_idJetRecoHMBeforeAlg;//vettori di id jet assegnati al vertice con massima molteplicita'
    vector< int > vec_idJetSimHMBeforeAlg;
   
    std::map<int,int>::const_iterator mapBeforeAlgIt = map_Vmin.begin();
    for(; mapBeforeAlgIt!=  map_Vmin.end(); ++mapBeforeAlgIt ){
      // se il vertice e' quello a massima molteplicita' mi salvo le id dei jet ad esso
      // associati per confrontarli con quelli assegnati al vertice con massima 
      // molteplicita' simulato
      if(mapBeforeAlgIt->second == idPVBeforeAlg){
	vec_idJetRecoHMBeforeAlg.push_back(mapBeforeAlgIt->first);
      }
    }
    cout<<"Reco HM #Jet"<<vec_idJetRecoHMBeforeAlg.size()<<", "<<nJBeforeAlg<<endl;
    
    std::map<int,int>::const_iterator mapSimBeforeAlgIt = map_VminSim.begin();
    for(; mapSimBeforeAlgIt!=  map_VminSim.end(); ++mapSimBeforeAlgIt ){
      if(mapSimBeforeAlgIt->second == idPVSimBeforeAlg){
	vec_idJetSimHMBeforeAlg.push_back(mapSimBeforeAlgIt->first);
      }
    }
    cout<<"Sim HM #Jet"<<vec_idJetSimHMBeforeAlg.size()<<", "<<nJSimBeforeAlg<<endl;

    //confronto tra le due collezioni di jet per contare quanti sono gli stessi 
    //jet assegnati al vertice con massima molteplicita'
    //---------------------------------------------
    int sameJetCounterBeforeAlg = 0;

    std::vector< int >::const_iterator  vec_idJetSimHMBeforeAlgIt = vec_idJetSimHMBeforeAlg.begin();
    for(; vec_idJetSimHMBeforeAlgIt != vec_idJetSimHMBeforeAlg.end(); vec_idJetSimHMBeforeAlgIt++){
    std::vector< int >::const_iterator  vec_idJetRecoHMBeforeAlgIt = vec_idJetRecoHMBeforeAlg.begin();
       for(; vec_idJetRecoHMBeforeAlgIt != vec_idJetRecoHMBeforeAlg.end(); vec_idJetRecoHMBeforeAlgIt++){
	 //se l'id del jet ricostruito e' uguale all'id del simulato incremento il contatore
	 if( *vec_idJetRecoHMBeforeAlgIt == *vec_idJetSimHMBeforeAlgIt) sameJetCounterBeforeAlg++;
       }
    }

    double WeightBeforeAlg = 0.;
    if(vec_idJetRecoHMBeforeAlg.size() > vec_idJetSimHMBeforeAlg.size()){  
      WeightBeforeAlg = (double)sameJetCounterBeforeAlg/(double)vec_idJetRecoHMBeforeAlg.size() ;
    }else{
      WeightBeforeAlg = (double)sameJetCounterBeforeAlg/(double)vec_idJetSimHMBeforeAlg.size() ;
    }

    cout<<"sameJetCounterBeforeAlg: "<<sameJetCounterBeforeAlg<<" , WeightBeforeAlg: "<<WeightBeforeAlg<<endl;

    H_nJSimVsnJRecoVtxHMBeforeAlgWeighted_->Fill(nJSimBeforeAlg, nJBeforeAlg, WeightBeforeAlg);




  } //end if !caloJet.empty()  

  //=====ROBERTO=================================
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void AlgValidator::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void AlgValidator::endJob() {
  
  
  //=====ROBERTO=============

  OutputFile->cd();

  //Histograms for Njet  
  H_NJet_withTracks_->Write();
  H_NJet_withNoTracks_->Write();
  
  //Histograms for Dz jet - vertices 
  H_Vmin_->Write();
  //Histogram Z jet  
  H_ZJetWAvg_->Write();
  H_ZJetWAvgError_->Write();
  H_DzJetVtxWAvg_->Write();
  H_DzJetVtxWAvg2Best_ ->Write();
  H_DzJetVtxWAvg_vs_2Best_->Write();

  //Histograms of id vertex assigned to a jet after rejection algorithm
  H_idRJAlg_->Write();  
  //Histograms of Z jet and error after the rejection algorithm (5 sigma)
  H_ZJet_RJAlg_->Write();
  H_ZErrorJet_RJAlg_->Write();
  H_DzJetVtx_RJAlg_->Write();
  H_Dz2JetVtx_RJAlg_->Write();
  H_DzJetVtx_RJAlg_vs_2Best_->Write();

  //Histograms of fraction between probability of best vertex association and second best  
  H_P1P2_->Write();

  //Histograms of id vertex assigned to a jet after SAAlg
  H_idVtxSAAlg_->Write();
  
  //Histograms of difference of log of probability of best vertex association and second best vs 
  //difference of log of probability considering the vertexes jet multeplicity
  H_pVsp_->Write();
  //Histograms of # of vertexes with a specified jet multeplicity
  H_nVtxMap_->Write();
  H_nVtxMapRJAlg_->Write();
  H_nVtxMapSAAlg_->Write();
  H_nVtxMapProbEval_->Write();
  
  //Histograms of probability  of jet assigned to vertexes configuration
  H_ProbRJAlg_->Write();
  H_ProbAfterFAAlg_->Write();
  H_ProbAfterSAAlg_->Write();

  //Histograms of id vertex assigned to a jet after first association
  H_idVtxFAAlg_->Write();
  //Histograms of Dz jet - vertex after SAAlg
  H_DzAfterSAAlg_->Write();

  //Histograms of # Tracks for jet with #DeltaZ best
  H_nTrack_RJAlg_->Write();

  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after SAAlg
  H_nJMaxVsNjTot_SAAlg_->Write();
  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after RJAlg
  H_nJMaxVsNjTot_RJAlg_->Write();

  //Histograms of # Tracks for jet with #DeltaZ best < 5#sigma vs z jet (error)
  H_nTkRJAlgVSzJetError_->Write();
  H_nTkRJAlgVSzJet_->Write();

  //Histograms of #jet with and without tracks (only cinematics cuts)
  H_nJet_->Write();

  //Histograms of #jet before, after RJAlg and SAAlg (with tracks)
  H_nJMax_SAAlg_->Write();
  H_nJMax_RJAlg_->Write();
  H_nJMax_BeforeAlg_->Write();
  H_nJMax_ProbEval_->Write();
  
  //Histograms of Minimum Dz between vertexes
  H_VtxMinDz_->Write();

  H_ProfileVtxVsNjBeforeAlg_->Write();
  H_ProfileVtxVsNjRJAlg_->Write();
  H_ProfileVtxVsNjSAAlg_->Write();
  H_ProfileVtxVsNjProbEval_->Write();

  H_ZVtx_->Write();
  H_ZVtxError_->Write();

  H_DzAfterProbEval_->Write();

  H_ZBeforevsAfterRJAlg_->Write();

  H_ZPvVsZBeforeAlg_->Write();
  H_ZPvVsZRJAlg_->Write();
  H_ZPvVsZProbEval_->Write();

  H_ZPvVsZRecoBeforeAlg_->Write();
  H_ZPvVsZRecoRJAlg_->Write();
  H_ZPvVsZRecoProbEval_->Write();

  H_ZjetSimVsReco_->Write();		
  H_dzVtxSimRecoBeforeAlg_->Write();	
  H_dzVtxSimRecoBeforeAlgVsNtk_->Write();	
  H_ZjetSimVsReco_RJAlg_->Write();	
  H_dzVtxSimRecoRJAlg_->Write();		
  H_dzVtxSimRecoRJAlgVsNtk_->Write();	
  H_dzVtxSimRecoProbEval_->Write();      


  H_dzVtxHMSimRecoProbEval_->Write();
  H_dzVtxHMSimRecoBeforeAlg_->Write();
  H_dzVtxHMSimRecoRJAlg_->Write();

  H_nJSimVsnJRecoVtxHMProbEval_->Write();
  H_nJSimVsnJRecoVtxHMBeforeAlg_->Write();
  H_nJSimVsnJRecoVtxHMRJAlg_->Write();

  H_nJSimVsnJRecoVtxHMProbEvalWeighted_->Write();
  H_nJSimVsnJRecoVtxHMBeforeAlgWeighted_->Write();
  H_nJSimVsnJRecoVtxHMRJAlgWeighted_->Write();



  //=====ROBERTO=============
  
}


// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(AlgValidator);
