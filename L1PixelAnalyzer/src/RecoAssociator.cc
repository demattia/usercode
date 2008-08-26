//////////////////////////////////////////////////////////////////////////////
//
// RecoAssociator.cc
// Code extracted from ProbAssoc.cc, 29/04/08 R.Casagrande
//
//  
//
// --------------------------------------------------------------------------------
//
// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/RecoAssociator.h"

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

// //For GenJet
// #include "DataFormats/JetReco/interface/GenJet.h"

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

// #include "AnalysisExamples/AnalysisClasses/interface/SortedAssocAlg.h"

// #include "AnalysisExamples/AnalysisObjects/interface/AssociationSummary.h"

#include "AnalysisExamples/AnalysisClasses/interface/VerticesJetCounter.h"

#include "AnalysisExamples/AnalysisClasses/interface/ProbCal.h"

// #include "AnalysisExamples/AnalysisClasses/interface/FirstAssocAlgWithError.h"

#include "AnalysisExamples/AnalysisClasses/interface/MinimumDzVtx.h"

#include "AnalysisExamples/AnalysisClasses/interface/ProbEval.h"



//======ROBERTO========================
// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------


// Constructors and destructor
// ---------------------------
RecoAssociator::RecoAssociator(const edm::ParameterSet& iConfig) :
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

  //Histogram of Dz selected sim tracks
  H_Dz_sim_selected_tk_ = new TH1D("H_Dz_sim_selected_tk", "Delta Z sim selected tracks " ,2000,-1.,1.);

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


RecoAssociator::~RecoAssociator()
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
void RecoAssociator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  if(SIM_){
    try{
      iEvent.getByLabel( simVtxLabel_, SimVtx );
    }
    catch(...){
      std::cerr << "Could not find the SimVertex collection" << endl;
    }
    simVertexes = *(SimVtx.product());
  }
  
  //------------------------------
  // Pile-up vertices collection
  //------------------------------
  
  edm::Handle<HepMCProduct> theFamosPileUp;
  const HepMC::GenEvent* myGenPUEvent = 0;
  if(SIM_){
    try{
      iEvent.getByLabel( simPUVtxLabel_ ,theFamosPileUp);
    } 
    catch (...) {
      std::cerr << "Could not find the pu collection" << std::endl;
      return;
    } 
  }
  //--------------------------
  // Reco vertex collection
  //--------------------------

  Handle<reco::VertexCollection> recVtxs;
  std::vector<reco::Vertex> recoVertexes;
  if(!SIM_){
    try {
      iEvent.getByLabel( vtxSample_ , recVtxs);
    } 
    catch (...) {
      std::cerr << "Could not find the recovertex collection" << std::endl;
      return;
    } 
    recoVertexes = *(recVtxs.product());
  }
  //=======================================
  //------------------------------  
  // Simulated tracks collection 
  //------------------------------
  
  //simtrack per l'associazione tracce simulate vertici simulati 
  // Take the SimTracks collection
  
  Handle<SimTrackContainer> SimTk;
  std::vector<SimTrack> theSimTracks;
  
  if(QCD_){  
    try {
      iEvent.getByLabel( simTkLabel_, SimTk );
      //  cout << "SimTk size = " << SimTk->size() << endl;
    }
    catch (...) {
      std::cerr << "Could not find the SimTracks collection" << std::endl;
      return;
    } 
    theSimTracks = *(SimTk.product());
  }

  //==========================================
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

  //   edm::Handle<reco::GenJetCollection> genJets;

  //   try {
  //     iEvent.getByLabel( genJetsLabel_, genJets );
  //   } 
  //   catch (...) {
  //     std::cerr << "Could not find the genjet collection" << std::endl;
  //     return;
  //   } 

  //======ROBERTO=======================
  
  
  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }
  
  //======ROBERTO================================
  
  ////////////////////////////////// MC analysis ////////////////////////////
  //vettore di vertici, in cui il primo e' il vertice primario
  std::vector< math::XYZTLorentzVector > vec_Vertices;
  //metto via un vettore col vertice primario dell'evento

  if(SIM_){
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
  }

  if(SIM_){
    std::vector<math::XYZTLorentzVector>::const_iterator vec_VerticesIt = vec_Vertices.begin();
    for(; vec_VerticesIt != vec_Vertices.end(); ++vec_VerticesIt){
      H_ZVtx_->Fill(vec_VerticesIt->z());
    }
    
    //--------------------------------
    //Minimum Dz Vertices
    //--------------------------------
    MinimumDzVtx<math::XYZTLorentzVector> recoVtxMinDz(vec_Vertices);
    std::vector<double> vec_VtxMinDz(recoVtxMinDz.evaluate());
    
    std::vector<double>::const_iterator vec_VtxMinDzIt =  vec_VtxMinDz.begin();
    for(; vec_VtxMinDzIt != vec_VtxMinDz.end(); ++vec_VtxMinDzIt ){
      
      H_VtxMinDz_->Fill(*vec_VtxMinDzIt);
      
    }
  }  

  if(!SIM_){
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
  }  

  //=================================================================
  
 
  //--------------------------
  //Simulated tracks vertices
  //-------------------------- 
  //faccio un loop sul vettore di tracce simulate per estrarre le informazioni sulle Z della traccia
  //Come da esempio in AnalysisExamples/SimTrackerAnalysis/src/SimTrackSimVertexAnalyzer.cc
  
  //Vector of simTracks' Z, eta, phi  
  SimpleTrackCollection vec_SimTk;
  if(QCD_){ 
    double simVertexes_size = 0.;
    simVertexes_size =  simVertexes.size();

    for (int isimvtx = 0; isimvtx < simVertexes_size;isimvtx++){  //, simvtx_it++){  
      for (std::vector<SimTrack>::iterator isimtk = theSimTracks.begin();
	   isimtk != theSimTracks.end(); ++isimtk){
	double Z_simTk = -1000.;
	double Phi_simTk = 0.;
	double Eta_simTk = 0.;
	double Pt_simTk = 0.;
      
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
	  vec_SimTk.push_back(SimpleTrack(Pt_simTk, Eta_simTk, Phi_simTk, Z_simTk, 0., 0.));

	}//end if track associated to this vertex
      }//end loop on sim tracks
    }//end loop on sim vertexes
  }

  //=================================================================
    
    
  //--------------------------------
  //Reco tracks and Offline calojet (classe per simple calojet +z error e significanza traccie)
  //--------------------------------  

  std::map<int,int> map_Vmin; //map of id best Z vertex in association with dz minimum criterium 
  std::map<int,int> map_idRJAlg; //map of id best Z vertex in association with dz minimum criterium after tracks rejection

  //Et, eta, phi Offline caloJet
  double Et_Jet = 0.;
  double eta_Jet = 0.;
  double phi_Jet = 0.;

  SimpleCaloJetCollection vec_CaloJet; //Collection of simpleCalojet passing al cuts

  int Njet = 0 ; //# jet
  int NjetBeforeAlg = 0 ; //# jet before cuts
  int NJet_withTracks = 0; //# jet with assigned tracks
  int NJet_withNoTracks = 0; //# jet without assigned tracks
  int Ntk = 0; //# tracks

  map< int, pair< double , double > > probDzMap; //map of probability idJet is the key, value is a pair(Gaussian Prob best dz,Gaussian Prob second best dz)

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

      //=================================================================
      //----------------------------------------------
      //Simtracks - reco selected tracks association 
      //Delta2 minimum simTrack - reco selected track (ok)
      //----------------------------------------------
      if(QCD_){
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
      
	//--------------------------------------------
	//Delta Z tra le tracce simulate selezionate (ok)
	//--------------------------------------------
	//In un doppio loop su tutte le traccie simulate selezionate calcolo tutte le deltaZ partendo con l'iteratore 
	//del secondo loop dal secondo, per evitare la differenza tra gli stessi elementi (equivalente alla condizione it2>it1)
	anaobj::SimpleTrackCollection::const_iterator vec_SimTk_iter_1 = vec_SimTk_selected.begin();
	for(; vec_SimTk_iter_1 !=  vec_SimTk_selected.end(); vec_SimTk_iter_1++ ){
	  anaobj::SimpleTrackCollection::const_iterator vec_SimTk_iter_2 = (vec_SimTk_iter_1+1);
	  for(; vec_SimTk_iter_2 !=  vec_SimTk_selected.end(); vec_SimTk_iter_2++ ){
	    
	    double Dz_sim_selected_tk = 0.; //Dz selected sim tracks
	    //calcolo il Dz tra traccie
	    Dz_sim_selected_tk = vec_SimTk_iter_1->z() - vec_SimTk_iter_2->z();
	    
	    //Histogram of Dz selected sim tracks
	    H_Dz_sim_selected_tk_->Fill(Dz_sim_selected_tk);
	    
	  }//end second loop on Z selected sim tracks 
	}//end first loop on Z selected sim tracks 
	
      }//end if(QCD_)

      //=================================================================

      //----------------------------------------------------
      // Jet Z: Wieghted Average on tracks' Z, with tracks
      //        Z error
      //----------------------------------------------------
      
      //richiamo la funzione wAverager per calcolare la media pesata delle tracce ricostruite
      //richiede in input una simpletrackcollection e restituisce una coppia pair(contatore # 
      //tracce, pair (media, errore))
      //----------------------------------------------------------------------------------------
      pair< int, pair< double, double > > wAvgRecoTk;
      //per le tracce simulate l'errore sulla z e' 0. -> la media pesata degli errori potrebbe non funziona -> media semplice
      if(QCD_) wAvgRecoTk = wAverager(vec_SimTk_selected);

      if(!QCD_) wAvgRecoTk = wAverager(vec_RecoTk);

      double wAvgRecoTkValue = wAvgRecoTk.second.first; //coordinata z media delle tracce
      double  wAvgRecoTkError = wAvgRecoTk.second.second; //errore sulla media pesata
      int NTks_S3 = wAvgRecoTk.first ; //contatore per il numero di tracce con significanza minore di 3

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
      // Delta Z minimo tra jet e vertici (reco)
      // con Z media pesata tracce
      //---------------------------------
      //if # of selected tracks with IP Significance < 3. != 0
      if (NTks_S3 != 0){
	
	//Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata e i vertici (PV e PU)	
	//la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
	//la funzione vuole in input un XYZTLorentz vector e un double
	//----------------------------------------------------------------------------------------------------------------------
	pair<pair<int, double>, pair<int, double> > dzFirstSecond;
	  
	if(SIM_) dzFirstSecond =  dzAssociator(vec_Vertices, wAvgRecoTkValue) ;
	  
	if(!SIM_) dzFirstSecond = dzAssociator(recoVertexes, wAvgRecoTkValue) ;

	pair<int, double> * first = &(dzFirstSecond.first);
	pair<int, double> * second = &(dzFirstSecond.second);
	int Vmin = first->first;
	double Dz_Jet_Vtx = first->second;
	double Dz_Jet_Vtx_2Best = second->second;
	
	//Fill map with idVertex associated with minimum dz criterium with WAvg z jet 
	map_Vmin.insert(make_pair( NjetBeforeAlg, Vmin) );
	
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
      SimpleTrackCollection recoTkClosest;
      if(QCD_){
	RejectionAlg tkRejection(vec_SimTk_selected);
	recoTkClosest = (tkRejection.eval(sigmaCut));
      }

      if(!QCD_){
	RejectionAlg tkRejection(vec_RecoTk);
	recoTkClosest = (tkRejection.eval(sigmaCut));
      }

      //chiamo la funzione che fa la media pesata wAvenger, che restituisce una pair(# tracce, pair(value, error))
      pair< int, pair< double, double > > zJetRejectionAlg( wAverager(recoTkClosest));
      //z jet e suo errore
      double zJet = zJetRejectionAlg.second.first;
      double zJetError = zJetRejectionAlg.second.second;

      int NTksRJAlg_S3 = zJetRejectionAlg.first ; //contatore per il numero di tracce con significanza minore di 3	
      //sono nello stesso numero di dei jet prima della reiezione, perche' quelli con 2 o meno tracce skippano l'algoritmo
      //e gli altri hanno tracce comunque, cambiera' il contenuto in tracce dei jet
      if(NTksRJAlg_S3 != 0){
	//riempio gli istogrammi di value e error (Z del jet calcolata e suo errore) 
	//Histograms of Z jet and error after the rejection algorithm (5 sigma)
	H_ZJet_RJAlg_->Fill(zJet);
	H_ZErrorJet_RJAlg_->Fill(zJetError);
      }

      //---------------------------------
      // Delta Z minimo tra jet e vertici (reco)
      // con Z media pesata tracce
      //---------------------------------
	
      //Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata dopo la reiezione e i vertici (PV e PU)	
      //la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
      //la funzione vuole in input un XYZTLorentz vector e un double
      //----------------------------------------------------------------------------------------------------------------------
      //if # of selected tracks with IP Significance < 3. != 0
      if (NTksRJAlg_S3 != 0){	
	pair<pair<int, double>, pair<int, double> > dzAfterRejection;

	if(SIM_) dzAfterRejection = dzAssociator(vec_Vertices,zJetRejectionAlg.second.first ) ;

	if(!SIM_) dzAfterRejection = dzAssociator(recoVertexes,zJetRejectionAlg.second.first ) ;

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

	//mappa  con #jet,come key, coppia p1 e p2 come value
	probDzMap.insert(make_pair(Njet , make_pair( p1 , p2  ) ) );
	
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
  if( !vec_CaloJet.empty()){

    //--------------------------------------
    // #jet before any algorithm 
    //--------------------------------------
    //------------------------------------------------------------------------------
    // # vertex with # of jet specified by key
    //------------------------------------------------------------------------------
    // VerticesJetCounter.countJet: count # jet to every vertex id specified by key
    // VerticesJetCounter.countVtx: count # vertex with # of jet specified by key
    //------------------------------------------------------------------------------ 

    VerticesJetCounter JetCountBeforeAlg;
    std::map<int,int> nJetMapBeforeAlg(JetCountBeforeAlg.countJet(map_Vmin));
    std::map<int,int> nVtxMapBeforeAlg(JetCountBeforeAlg.countVtx(nJetMapBeforeAlg));

    if( !nVtxMapBeforeAlg.empty() )  H_nJMax_BeforeAlg_->Fill((*nVtxMapBeforeAlg.rbegin()).first + NJet_withNoTracks);


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

      if(SIM_)  zVtx = vec_Vertices[idVtx].z();

      if(!SIM_) zVtx = recoVertexes[idVtx].z();

      double zJet = vec_CaloJet[idJet].z();

      //Histograms of Dz jet - vertex after SAAlg
      H_DzAfterProbEval_->Fill(zVtx - zJet);
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
    
    //       //--------------------------------------------------------------------
    //       //Configuration probability after association (RJAlg)
    //       //--------------------------------------------------------------------
    //       // ProbCal.computeEventProb: combined probability of configuration of
    //       // # vertexes with specified jet multeplicity
    //       //--------------------------------------------------------------------

    //       //      ProbCal RJAssocEventProb( pMatrixSing25_ );
    //       ProbCal RJAssocEventProb;
    //       double ProbRJAlg(RJAssocEventProb.computeEventProb(nVtxMapRJAlg));
    //       //Histograms of probability  of jet assigned to vertexes configuration
    //       H_ProbRJAlg_->Fill(TMath::Exp(ProbRJAlg));

    //       //--------------------------------------------------------------------------
    //       // First Association Jet - Vertices after rejection alg
    //       //---------------------------------------------------------------------------
    //       // FirstAssocAlgWithError.associate(int): association done jet - vertexes 
    //       // if dz second best is > i(nt) * sigma so that probability P2 is 0
    //       //---------------------------------------------------------------------------
    //       AssociationSummaryCollection firstAssocMap;

    //       if(SIM_) {
    // 	FirstAssocAlgWithError<math::XYZTLorentzVector> firstAssociation(vec_CaloJet , vec_Vertices);
    // 	firstAssocMap = firstAssociation.associate(probDzMap);
    //       }
    //       if(!SIM_){
    // 	FirstAssocAlgWithError<reco::Vertex> firstAssociation(vec_CaloJet , recoVertexes);
    // 	firstAssocMap = firstAssociation.associate(probDzMap);
    //       }

    //       AssociationSummaryCollection::const_iterator firstAssocMapIt = firstAssocMap.begin();
    //       for(; firstAssocMapIt != firstAssocMap.end(); firstAssocMapIt++){
    // 	if (firstAssocMapIt->Associated()){
    // 	  //Histograms of id vertex assigned to a jet after first association
    // 	  H_idVtxFAAlg_->Fill(firstAssocMapIt->idVtx());
    // 	}
    //       }

    //       //------------------------------------------------------------------------------
    //       // # vertex with # of jet specified by key (FAAlg)
    //       //------------------------------------------------------------------------------
    //       // VerticesJetCounter.countJet: count # jet to every vertex id specified by key
    //       // VerticesJetCounter.countVtx: count # vertex with # of jet specified by key
    //       //------------------------------------------------------------------------------

    //       VerticesJetCounter JetCounting;
    //       std::map<int,int> nJetMap(JetCounting.countJet(firstAssocMap));

    //       std::map<int,int> nVtxMap(JetCounting.countVtx(nJetMap));

    //       std::map<int,int>::const_iterator nVtxMapIt = nVtxMap.begin();
    //       for(; nVtxMapIt!=  nVtxMap.end(); ++nVtxMapIt ){
    // 	//Histograms of # of vertexes with a specified jet multeplicity
    // 	H_nVtxMap_->Fill(nVtxMapIt->first,(double)nVtxMapIt->second);
    //       }

    //       //--------------------------------------------------------------------
    //       //Configuration probability after first association
    //       //--------------------------------------------------------------------
    //       // ProbCal.computeEventProb: combined probability of configuration of
    //       // # vertexes with specified jet multeplicity
    //       //--------------------------------------------------------------------
    //       ProbCal firstAssocEventProb;
    //       double ProbBefore(firstAssocEventProb.computeEventProb(nVtxMap));

    //       //Histograms of probability  of jet assigned to vertexes configuration 
    //       H_ProbAfterFAAlg_->Fill(TMath::Exp(ProbBefore));

    //       //-------------------------------------------------------------------
    //       //Iterative algorithm for association
    //       //-------------------------------------------------------------------
    //       // SortedAssocAlg:for jets not associated with FirstAssocAlg 
    //       // the combined probability of dz and (m+1, n) configuration 
    //       // vs dz secondbest and (m, n+1) are evaluated to make association 
    //       // jet are sorted for decrescent value of p1/p2 and every jet 
    //       //association modify (m,n) configuration 
    //       //-------------------------------------------------------------------
    //       SortedAssocAlg sortedAssocP1P2( pMatrix25_ );
    //       AssociationSummaryCollection assocMap( sortedAssocP1P2.associate(firstAssocMap, probDzMap) );
    //       //  cout<<"assocMapSize= "<<assocMap.size()<<endl; 
 
    //       AssociationSummaryCollection::const_iterator assocMapIt = assocMap.begin();
    //       for(; assocMapIt != assocMap.end(); assocMapIt++){
    // 	//Histograms of id vertex assigned to a jet after SAAlg
    // 	H_idVtxSAAlg_->Fill(assocMapIt->idVtx());
  
    // 	if(assocMapIt->SecondBest() != 1000.){ //1000. valore fake del second best che distingue i jet assegnati dal FirstAssocAlg, in quanto il second best incompatibile
    // 	  //Histograms of difference of log of probability of best vertex association and second best vs 
    // 	  //difference of log of probability considering the vertexes jet multeplicity
    // 	  H_pVsp_->Fill(TMath::Log(probDzMap[assocMapIt->idJet()].first ) - TMath::Log(probDzMap[assocMapIt->idJet()].second ) ,assocMapIt->Best() -  assocMapIt->SecondBest() );
    // 	}
    //      }

    //       //------------------------------------------------------------------------------
    //       // # vertex with # of jet specified by key (SAAlg)
    //       //------------------------------------------------------------------------------
    //       // VerticesJetCounter.countJet: count # jet to every vertex id specified by key
    //       // VerticesJetCounter.countVtx: count # vertex with # of jet specified by key
    //       //------------------------------------------------------------------------------

    //       VerticesJetCounter JetCountSAAlg;
    //       std::map<int,int> nJetMapSAAlg(JetCountSAAlg.countJet(assocMap));
    //       //  cout<<"nJetMapSAAlgSize= "<<nJetMapSAAlg.size()<<endl;
    //       std::map<int,int> nVtxMapSAAlg(JetCountSAAlg.countVtx(nJetMapSAAlg));
    //       //  cout<<"nVtxMapSAAlgSize= "<<nVtxMapSAAlg.size()<<endl;
  
    //       std::map<int,int>::const_iterator nVtxMapSAAlgIt = nVtxMapSAAlg.begin();
    //       for(; nVtxMapSAAlgIt!=  nVtxMapSAAlg.end(); ++nVtxMapSAAlgIt ){
    // 	//Histograms of # of vertexes with a specified jet multeplicity
    // 	H_nVtxMapSAAlg_->Fill(nVtxMapSAAlgIt->first,(double)nVtxMapSAAlgIt->second);
    //       }
    
    //       //Histograms of # jet for vertex with max multeplicity  vs # jet tot after SAAlg
    //       if(!nVtxMapSAAlg.empty()) H_nJMaxVsNjTot_SAAlg_->Fill( (*nVtxMapSAAlg.rbegin()).first, Njet );
    
    //       if(!nVtxMapSAAlg.empty()) H_nJMax_SAAlg_->Fill((*nVtxMapSAAlg.rbegin()).first + NJet_withNoTracks);
    
    //       //--------------------------------------------------------------------
    //       //Configuration probability after association (SAAlg)
    //       //--------------------------------------------------------------------
    //       // ProbCal.computeEventProb: combined probability of configuration of
    //       // # vertexes with specified jet multeplicity
    //       //--------------------------------------------------------------------
    //       ProbCal AssocEventProb;
    //       double ProbAfter(AssocEventProb.computeEventProb(nVtxMapSAAlg));
    //       //Histograms of probability  of jet assigned to vertexes configuration 
    //       H_ProbAfterSAAlg_->Fill(TMath::Exp(ProbAfter));

    //       //-------------------------------------
    //       // DZ after association (SAAlg) 
    //       //-------------------------------------

    //       //loop on associationMap
    //       AssociationSummaryCollection::const_iterator assocMapIt1 = assocMap.begin();
    //       for(; assocMapIt1 != assocMap.end() ; ++assocMapIt1 ){
    
    // 	int idJet = assocMapIt1->idJet();
    // 	int idVtx = assocMapIt1->idVtx();
    
    // 	double zVtx ;

    // 	if(SIM_)  zVtx = vec_Vertices[idVtx].z();
    // 	if(!SIM_) zVtx = recoVertexes[idVtx].z();

    // 	double zJet = vec_CaloJet[idJet].z();

    // 	//Histograms of Dz jet - vertex after SAAlg
    // 	H_DzAfterSAAlg_->Fill(zVtx - zJet);
    //       }
    //-----------------------------
    //Dz vertices - primary vertex 
    //-----------------------------
    //cerco il vertice con max # jet assegnati	
    std::map<int , int> VtxMaxProbEval; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxBeforeAlg; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxRJAlg; //mappa con key il # di jet e value l'id del vertice
    //       std::map<int , int> VtxMaxSAAlg; //mappa con key il # di jet e value l'id del vertice

    //riempimento mappe con key il # di jet e value l'id del vertice
    std::map<int,int>::const_iterator nJetMapProbEvalIt = nJetMapProbEval.begin();
    for(; nJetMapProbEvalIt != nJetMapProbEval.end() ; nJetMapProbEvalIt++){

      VtxMaxProbEval.insert( make_pair( nJetMapProbEvalIt->second, nJetMapProbEvalIt->first ) );
	
    }

    std::map<int,int>::const_iterator nJetMapBeforeAlgIt = nJetMapBeforeAlg.begin();
    for(; nJetMapBeforeAlgIt != nJetMapBeforeAlg.end() ; nJetMapBeforeAlgIt++){

      VtxMaxBeforeAlg.insert( make_pair( nJetMapBeforeAlgIt->second, nJetMapBeforeAlgIt->first ) );
	
    }

    std::map<int,int>::const_iterator nJetMapRJAlgIt = nJetMapRJAlg.begin();
    for(; nJetMapRJAlgIt != nJetMapRJAlg.end() ; nJetMapRJAlgIt++){

      VtxMaxRJAlg.insert( make_pair( nJetMapRJAlgIt->second, nJetMapRJAlgIt->first ) );
	
    }
 
    //       std::map<int,int>::const_iterator nJetMapSAAlgIt = nJetMapSAAlg.begin();
    //       for(; nJetMapSAAlgIt != nJetMapSAAlg.end() ; nJetMapSAAlgIt++){

    // 	VtxMaxSAAlg.insert( make_pair( nJetMapSAAlgIt->second, nJetMapSAAlgIt->first ) );
	
    //       }

    //calcolo il dz
    int idVtx = 0;

    //id vertice con max num jet e sua Z (la mappa e' ordinata in valori crescenti 
    //della key che per queste mappe e' scelta esser il # di jet)
    int  idPVProbEval = VtxMaxProbEval.rbegin()->second;
    double zPVProbEval;
    if(SIM_)  zPVProbEval = vec_Vertices[idPVProbEval].z(); 
    if(!SIM_) zPVProbEval = recoVertexes[idPVProbEval].z(); 

    int  idPVBeforeAlg = VtxMaxBeforeAlg.rbegin()->second;
    double zPVBeforeAlg ;
    if(SIM_)  zPVBeforeAlg = vec_Vertices[idPVBeforeAlg].z(); 
    if(!SIM_) zPVBeforeAlg = recoVertexes[idPVBeforeAlg].z(); 

    int  idPVRJAlg = VtxMaxRJAlg.rbegin()->second;
    double zPVRJAlg;
    if(SIM_)  zPVRJAlg = vec_Vertices[idPVRJAlg].z(); 
    if(!SIM_)  zPVRJAlg = recoVertexes[idPVRJAlg].z();
 
    //       int  idPVSAAlg = VtxMaxSAAlg.rbegin()->second;
    //       double zPVSAAlg;
    //       if(SIM_)  zPVSAAlg = vec_Vertices[idPVSAAlg].z();
    //       if(!SIM_) zPVSAAlg = recoVertexes[idPVSAAlg].z();

    H_ZBeforevsAfterRJAlg_->Fill(zPVBeforeAlg,zPVRJAlg);
      
    if(SIM_){
      H_ZPvVsZBeforeAlg_->Fill(vec_Vertices.begin()->z(), zPVBeforeAlg);
      H_ZPvVsZRJAlg_->Fill(vec_Vertices.begin()->z(),zPVRJAlg);
      H_ZPvVsZProbEval_->Fill(vec_Vertices.begin()->z(),zPVProbEval);
    }
      
    //loop on vertex collection
    if(SIM_){
      std::vector<math::XYZTLorentzVector>::const_iterator vec_VerticesIt = vec_Vertices.begin();
      for(; vec_VerticesIt!=   vec_Vertices.end(); ++vec_VerticesIt ,idVtx++){

	double zVtx = vec_VerticesIt->z();
	//Associazione in base all'algoritmo di massima probabilita'
	if(idVtx != idPVProbEval){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxProbEval = zVtx - zPVProbEval;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int    nJetProbEval = nJetMapProbEval[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjProbEval_->Fill(tempDzVtxProbEval,(double) nJetProbEval);
	}
	//Associazione in base al minimo dz prima di qualsiasi algoritmo
	if(idVtx != idPVBeforeAlg){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxBeforeAlg = zVtx - zPVBeforeAlg;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int  nJetBeforeAlg = nJetMapBeforeAlg[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjBeforeAlg_->Fill(tempDzVtxBeforeAlg,(double) nJetBeforeAlg);
	}
	//Associazione in base al minimo dz dopo la reiezione delle tracce
	if(idVtx != idPVRJAlg){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxRJAlg = zVtx - zPVRJAlg;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int    nJetRJAlg = nJetMapRJAlg[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjRJAlg_->Fill(tempDzVtxRJAlg,(double) nJetRJAlg);
	}
	// 	  //Associazione in base all'algoritmo di probabilita'
	// 	  if(idVtx != idPVSAAlg){
	// 	    //Delta Z vertice - Vertice con max # jet
	// 	    double tempDzVtxSAAlg = zVtx - zPVSAAlg;
	// 	    //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	// 	    //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	// 	    int  nJetSAAlg = nJetMapSAAlg[idVtx];
	// 	    //riempio l'istogramma
	// 	    H_ProfileVtxVsNjSAAlg_->Fill(tempDzVtxSAAlg,(double) nJetSAAlg);
	// 	  }

      }//end loop on vertex collection
    }//end if(SIM_)
    if(!SIM_){
      std::vector<reco::Vertex>::const_iterator recoVertexesIt = recoVertexes.begin();
      for(; recoVertexesIt!=   recoVertexes.end(); ++recoVertexesIt ,idVtx++){

	double zVtx = recoVertexesIt->z();
	//Associazione in base all'algoritmo di massima probabilita'
	if(idVtx != idPVProbEval){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxProbEval = zVtx - zPVProbEval;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int    nJetProbEval = nJetMapProbEval[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjProbEval_->Fill(tempDzVtxProbEval,(double) nJetProbEval);
	}
	//Associazione in base al minimo dz prima di qualsiasi algoritmo
	if(idVtx != idPVBeforeAlg){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxBeforeAlg = zVtx - zPVBeforeAlg;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int  nJetBeforeAlg = nJetMapBeforeAlg[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjBeforeAlg_->Fill(tempDzVtxBeforeAlg,(double) nJetBeforeAlg);
	}
	//Associazione in base al minimo dz dopo la reiezione delle tracce
	if(idVtx != idPVRJAlg){
	  //Delta Z vertice - Vertice con max # jet
	  double tempDzVtxRJAlg = zVtx - zPVRJAlg;
	  //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	  //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	  int    nJetRJAlg = nJetMapRJAlg[idVtx];
	  //riempio l'istogramma
	  H_ProfileVtxVsNjRJAlg_->Fill(tempDzVtxRJAlg,(double) nJetRJAlg);
	}
	// 	  //Associazione in base all'algoritmo di probabilita'
	// 	  if(idVtx != idPVSAAlg){
	// 	    //Delta Z vertice - Vertice con max # jet
	// 	    double tempDzVtxSAAlg = zVtx - zPVSAAlg;
	// 	    //# jet assegnati a quel vertice NB. l'operatore [] inserisce una nuova key con value 0 
	// 	    //se idVtx non e' una key della mappa gia' presente quindi le mappe vengono modificate
	// 	    int  nJetSAAlg = nJetMapSAAlg[idVtx];
	// 	    //riempio l'istogramma
	// 	    H_ProfileVtxVsNjSAAlg_->Fill(tempDzVtxSAAlg,(double) nJetSAAlg);
	// 	  }

      }//end loop on vertex collection
    }//end if(!SIM_)
  } //end if !caloJet.empty()  

  //=====ROBERTO=================================
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void RecoAssociator::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void RecoAssociator::endJob() {
  
  
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

  //Histogram of Dz selected sim tracks
  H_Dz_sim_selected_tk_->Write();

  //=====ROBERTO=============
  
}


// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(RecoAssociator);
