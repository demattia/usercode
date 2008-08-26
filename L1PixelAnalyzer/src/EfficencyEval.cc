//////////////////////////////////////////////////////////////////////////////
//
// EfficencyEval.cc
// Code extracted from ProbAssoc.cc, 29/04/08 R.Casagrande
//
//  
//
// --------------------------------------------------------------------------------
//
// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/EfficencyEval.h"

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

// #include "AnalysisExamples/AnalysisClasses/interface/MinimumDzVtx.h"

#include "AnalysisExamples/AnalysisClasses/interface/ProbEval.h"

#include "AnalysisExamples/AnalysisObjects/interface/BaseVertex.h"



//======ROBERTO========================
// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------


// Constructors and destructor
// ---------------------------
EfficencyEval::EfficencyEval(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  
  //--------------------
  // Jet
  //--------------------
  
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),

  //=====ROBERTO================================
   
  //   //----------------------------------
  //   // Simulated Tracks and Vertices
  //   //----------------------------------
  
  //   //Simulated vertices
  //   simVtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx" ) ),
  
  //   //Simulated tracks
  //   simTkLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimTk" ) ),
  
  //   simPUVtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimPUVtx" ) ),
 
  //----------------------------------
  // Reco Tracks and Vertices
  //----------------------------------  

  //  //Reco tracks impact parameter
  //   impactParameterTagInfos( iConfig.getUntrackedParameter<std::string>( "impactParameterTagInfos" ) ),

  //   //Reco tracks impact parameter
  //   simpleTrackLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "simpleTracks" ) ),


  //Reco vertexs
  vtxSample_( iConfig.getUntrackedParameter<edm::InputTag>( "vtxSample" ) ),



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
  

  //Histograms of # of vertexes with a specified jet multeplicity
  H_nVtxMapRJAlg_ = new TH1D("H_nVtxMapRJAlg","nVtxMapRJAlg", 20, 0, 20);
  H_nVtxMapProbEval_ = new TH1D("H_nVtxMapProbEval","nVtxMapProbEval", 20, 0, 20);

 
  //Histograms of # Tracks for jet with #DeltaZ best 
  H_nTrack_RJAlg_ =  new TH1D("H_nTrack_RJAlg","# Tracks for jet with #DeltaZ best < 5#sigma", 50, 0, 50);

  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after RJAlg
  H_nJMaxVsNjTot_RJAlg_ =  new TH2D("H_nJMaxVsNjTot_RJAlg","# jet for vertex with max multeplicity  vs # jet tot after RJAlg" ,20,0,20,20,0,20);

  //Histograms of # Tracks for jet with #DeltaZ best vs z jet (error)
  H_nTkRJAlgVSzJetError_ =  new TH2D("H_nTkRJAlgVSzJetError","# Tracks for jet after RJAlg vs z jet error" ,20,0,20,20,0.,2.);
  H_nTkRJAlgVSzJet_ = new TH2D("H_nTkRJAlgVSzJet","# Tracks for jet with after RJAlg best vs z jet" ,20,0.,20.,200,0.,20.);

  //Histograms of #jet with and without tracks (only cinematics cuts)
  H_nJet_ = new TH1D("H_nJet", "# of jet (Et>25)" ,40,0.,40.);
 
  //Histograms of #jet before, after RJAlg and SAAlg (with tracks)
  H_nJMax_RJAlg_  =  new TH1D("H_nJMax_RJAlg", "# of jet (Et>25) assigned after RJAlg" ,20,0.,20.);
  H_nJMax_BeforeAlg_  =  new TH1D("H_nJMax_BeforeAlg", "# of jet (Et>25) assigned BeforeAlg" ,20,0.,20.);
  H_nJMax_ProbEval_ =  new TH1D("H_nJMax_ProbEval", "# of jet (Et>25) assigned after ProbEval" ,20,0.,20.);


  H_ProfileVtxVsNjBeforeAlg_ = new TProfile("H_ProfileVtxVsNjBeforeAlg","#DeltaZ vtx - PV vs Nj Before Algorithm",80,-20.,20.);
  H_ProfileVtxVsNjRJAlg_ = new TProfile("H_ProfileVtxVsNjRJAlg","#DeltaZ vtx - PV vs Nj after RJAlgorithm",80,-20.,20.);
  //  H_ProfileVtxVsNjSAAlg_ = new TProfile("H_ProfileVtxVsNjSAAlg","#DeltaZ vtx - PV vs Nj after SAAlgorithm",80,-20.,20.);
  H_ProfileVtxVsNjProbEval_ = new TProfile("H_ProfileVtxVsNjProbEval","#DeltaZ vtx - PV vs Nj after ProbEval",80,-20.,20.);

  H_ZVtx_ = new TH1D("H_ZVtx", "Z Vertice" ,40,-20.,20.);
  H_ZVtxError_ = new TH1D("H_ZVtxError", "Error Z Vtx" ,200,0.,0.2);

  H_DzAfterProbEval_ = new TH1D("H_DzAfterProbEval", "#DeltaZ jet - vertex after ProbEval" ,800,-4.,4.);
  
  H_ZBeforevsAfterRJAlg_ = new TH2D("H_ZBeforevsAfterRJAlg","Z prima vs dopo rejezione tracce" ,400,-20,20,400,-20,20);

  H_Et4jetProbEval_ = new TH1D("H_Et4jetProbEval","Et4jetProbEval",2000,0.,2000.);
  H_Et5jetProbEval_ = new TH1D("H_Et5jetProbEval","Et5jetProbEval",2000,0.,2000.);
  H_Et6jetProbEval_ = new TH1D("H_Et6jetProbEval","Et6jetProbEval",2000,0.,2000.);
    		    
  H_Et4jetBeforeAlg_ = new TH1D("H_Et4jetBeforeAlg","Et4jetBeforeAlg",2000,0.,2000.);
  H_Et5jetBeforeAlg_ = new TH1D("H_Et5jetBeforeAlg","Et5jetBeforeAlg",2000,0.,2000.);
  H_Et6jetBeforeAlg_ = new TH1D("H_Et6jetBeforeAlg","Et6jetBeforeAlg",2000,0.,2000.);


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


EfficencyEval::~EfficencyEval()
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
void EfficencyEval::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


  //--------------------------
  // Reco vertex collection
  //--------------------------

  Handle<BaseVertexCollection> recVtxs;
  BaseVertexCollection recoVertexes;
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


  //======ROBERTO=======================
  
  
  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }
  
  //======ROBERTO================================
  
  ////////////////////////////////// MC analysis ////////////////////////////

  if(!SIM_){
    BaseVertexCollection::const_iterator recoVertexesIt = recoVertexes.begin();
    for(; recoVertexesIt != recoVertexes.end(); ++recoVertexesIt){
      H_ZVtx_->Fill(recoVertexesIt->z());
      H_ZVtxError_->Fill(recoVertexesIt->zError());
    }
  }  
   
  //--------------------------------
  //Reco tracks and Offline calojet (classe per simple calojet +z error e significanza traccie)
  //--------------------------------  

  std::map<int,int> map_Vmin; //map of id best Z vertex in association with dz minimum criterium 
  std::map<int,int> map_idRJAlg; //map of id best Z vertex in association with dz minimum criterium after tracks rejection

  //Et, eta, phi Offline caloJet
  double Et_Jet = 0.;
  double eta_Jet = 0.;
  double phi_Jet = 0.;

  SimpleCaloJetCollection vec_CaloJet; //Collection of simpleCalojet passing all cuts
  SimpleCaloJetCollection vec_CaloJetWithNoTracks; //Collection of simpleCalojet with no tracks passing all cuts
  SimpleCaloJetCollection vec_RowCaloJet; //Collection of simpleCalojet passing only cuts on Et and eta

  int Njet = 0 ; //# jet
  int NjetBeforeAlg = 0 ; //# jet before cuts
  int NJet_withTracks = 0; //# jet with assigned tracks
  int NJet_withNoTracks = 0; //# jet without assigned tracks
  int Ntk = 0; //# tracks

  map< int, pair< double , double > > probDzMap; //map of probability idJet is the key, value is a pair(Gaussian Prob best dz,Gaussian Prob second best dz)

  //loop on tracks (IptagInfo) and caloJets 
  anaobj::OfflineJetCollection::const_iterator caloJets_it = caloJets->begin();
  //   reco::TrackIPTagInfoCollection::const_iterator TkIpTagInfo_it = iPtagInfos->begin();
  //   for( ; TkIpTagInfo_it != iPtagInfos->end(); ++TkIpTagInfo_it ) {
  for( ;  caloJets_it !=  caloJets->end(); ++caloJets_it ) {   
    //vector of reco selected tracks' Pt, eta, phi, Z
    SimpleTrackCollection vec_RecoTk;

    //Et, eta, phi Offline caloJet
    Et_Jet = caloJets_it->et();
    eta_Jet = caloJets_it->eta();
    phi_Jet =  caloJets_it->phi();

    // Take the selectedTracks vector
    // Return the vector of tracks for which the IP information is available Quality cuts are applied to reject fake tracks.
    //      const reco::TrackRefVector & vec_TkColl = TkIpTagInfo_it->selectedTracks();
    
    // Take the IP vector (ordered as the selectedTracks vector)
    //   const vector<reco::TrackIPTagInfo::TrackIPData> & vec_TkIP = TkIpTagInfo_it->impactParameterData();
    
    //loop on selected tracks collection and tracks IP  
    //      vector<reco::TrackIPTagInfo::TrackIPData>::const_iterator TkIP_it = vec_TkIP.begin();

    SimpleTrackRefVector assocTk( caloJets_it->tkRefVec() ) ;

    SimpleTrackRefVector::const_iterator TkColl_it = assocTk.begin();
    //    RefVector<reco::TrackCollection>::const_iterator TkColl_it = vec_TkColl.begin();
    for (; TkColl_it != assocTk.end(); ++TkColl_it){ //, ++TkIP_it ) {
      
      //2D Significance
      double Tk_IpS2D = 0.;
      Tk_IpS2D =(*TkColl_it)->ip2Dsignificance();
      //      cout<<"significanza 2D IP traccia: "<<Tk_IpS2D<<endl;
 
      //Reco Track Z
      double Tk_dz = 0.;
      Tk_dz = (*TkColl_it)->z();
           
      //Reco Track Z error
      double Tk_dz_error = 0.;
      Tk_dz_error= (*TkColl_it)->zError();
         
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
	
      //----------------------------------------------------
      // Jet Z: Wieghted Average on tracks' Z, with tracks
      //        Z error
      //----------------------------------------------------
      
      //richiamo la funzione wAverager per calcolare la media pesata delle tracce ricostruite
      //richiede in input una simpletrackcollection e restituisce una coppia pair(contatore # 
      //tracce, pair (media, errore))
      //----------------------------------------------------------------------------------------
      pair< int, pair< double, double > > wAvgRecoTk( wAverager(vec_RecoTk));
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
	  
	//	  if(SIM_) dzFirstSecond =  dzAssociator(vec_Vertices, wAvgRecoTkValue) ;
	  
	if(!SIM_) dzFirstSecond = dzAssociator(recoVertexes, wAvgRecoTkValue) ;

	pair<int, double> * first = &(dzFirstSecond.first);
	pair<int, double> * second = &(dzFirstSecond.second);
	int Vmin = first->first;
	double Dz_Jet_Vtx = first->second;
	double Dz_Jet_Vtx_2Best = second->second;
	
	//Fill map with idVertex associated with minimum dz criterium with WAvg z jet 
	map_Vmin.insert(make_pair( NjetBeforeAlg, Vmin) );
	


	//Histograms minimum Dz Jet - Vertices	
	H_Vmin_->Fill(Vmin);
	H_DzJetVtxWAvg_->Fill(Dz_Jet_Vtx);
	H_DzJetVtxWAvg2Best_->Fill(Dz_Jet_Vtx_2Best);
	H_DzJetVtxWAvg_vs_2Best_->Fill(fabs(Dz_Jet_Vtx),fabs(Dz_Jet_Vtx_2Best));


// 	vec_RowCaloJet.push_back(SimpleCaloJet( Et_Jet,
// 						eta_Jet,
// 						phi_Jet,
// 						wAvgRecoTkValue,
// 						wAvgRecoTkError,
// 						NjetBeforeAlg,
// 						0 ,
// 						0 ,
// 						0 ,
// 						0.,
// 						0.,
// 						0. ) );

	NjetBeforeAlg++;
      
      }	//end if # of selected tracks with IP Significance < 3. != 0
	
	//----------------------------------------------
	//Algoritmo di reiezione tracce ricostruite nel jet
	//Z jet calcolata con la media pesata, scartando 
	//le tracce piu' distanti di 5 sigma dal valore 
	//della media pesata con algoritmo iterativo
	//----------------------------------------------
	//chiamo l'algoritmo di reiezione (e' una classe con il metodo eval(# sigma taglio)) sulla collezione di tracce ricostruite
	//agisce se il numero di trace nel jet e' > 2
      RejectionAlg tkRejection(vec_RecoTk);
      SimpleTrackCollection recoTkClosest(tkRejection.eval(sigmaCut));

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

	//	  if(SIM_) dzAfterRejection = dzAssociator(vec_Vertices,zJetRejectionAlg.second.first ) ;

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
      if(NTksRJAlg_S3 == 0){
	vec_CaloJetWithNoTracks.push_back(SimpleCaloJet( Et_Jet,
							 eta_Jet,
							 phi_Jet,
							 0.,
							 0.,
							 0 ,
							 0 ,
							 0 ,
							 0 ,
							 0.,
							 0.,
							 0. ) );
      }


      vec_RowCaloJet.push_back(SimpleCaloJet( Et_Jet,
					      eta_Jet,
					      phi_Jet,
					      wAvgRecoTkValue,
					      wAvgRecoTkError,
					      0 ,
					      0 ,
					      0 ,
					      0 ,
					      0.,
					      0.,
					      0. ) );

    }// end cut Et(jet>25) eta<3
    
    //      ++caloJets_it;
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

      double zVtx = 0. ;

      //	if(SIM_)  zVtx = vec_Vertices[idVtx].z();

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

    //-----------------------------
    //Dz vertices - primary vertex 
    //-----------------------------
    //cerco il vertice con max # jet assegnati	
    std::map<int , int> VtxMaxProbEval; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxBeforeAlg; //mappa con key il # di jet e value l'id del vertice
    std::map<int , int> VtxMaxRJAlg; //mappa con key il # di jet e value l'id del vertice
    //      std::map<int , int> VtxMaxSAAlg; //mappa con key il # di jet e value l'id del vertice

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

    //calcolo il dz
    int idVtx = 0;

    //id vertice con max num jet e sua Z (la mappa e' ordinata in valori crescenti 
    //della key che per queste mappe e' scelta esser il # di jet)
    int  idPVProbEval = VtxMaxProbEval.rbegin()->second;
    int nJProbEval = VtxMaxProbEval.rbegin()->first;
    double zPVProbEval = 0. ;
    //      if(SIM_)  zPVProbEval = vec_Vertices[idPVProbEval].z(); 
    if(!SIM_) zPVProbEval = recoVertexes[idPVProbEval].z(); 

    int  idPVBeforeAlg = VtxMaxBeforeAlg.rbegin()->second;
    int nJBeforeAlg = VtxMaxBeforeAlg.rbegin()->first;
    double zPVBeforeAlg = 0. ;
    //      if(SIM_)  zPVBeforeAlg = vec_Vertices[idPVBeforeAlg].z(); 
    if(!SIM_) zPVBeforeAlg = recoVertexes[idPVBeforeAlg].z(); 

    int  idPVRJAlg = VtxMaxRJAlg.rbegin()->second;
    int nJRJAlg = VtxMaxRJAlg.rbegin()->first;
    double zPVRJAlg = 0. ;
    //      if(SIM_)  zPVRJAlg = vec_Vertices[idPVRJAlg].z(); 
    if(!SIM_)  zPVRJAlg = recoVertexes[idPVRJAlg].z();

    H_ZBeforevsAfterRJAlg_->Fill(zPVBeforeAlg,zPVRJAlg);

    if(!SIM_){
      BaseVertexCollection::const_iterator recoVertexesIt = recoVertexes.begin();
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

      }//end loop on vertex collection
    }//end if(!SIM_)

    vector< int > vec_idJetRecoHMProbEval;//vettori di id jet assegnati al vertice con massima molteplicita'

    std::map<int,int>::const_iterator mapProbEvalIt1 = mapProbEval.begin();
    for(; mapProbEvalIt1 !=  mapProbEval.end(); ++mapProbEvalIt1 ){
      // se il vertice e' quello a massima molteplicita' mi salvo le id dei jet ad esso
      // associati per confrontarli con quelli assegnati al vertice con massima 
      // molteplicita' simulato
      if(mapProbEvalIt1->second == idPVProbEval){

	vec_idJetRecoHMProbEval.push_back(mapProbEvalIt1->first);
	
      }
    }
    //    cout<<"Reco HM #Jet"<<vec_idJetRecoHMProbEval.size()<<", "<<nJProbEval<<endl;
    
    SimpleCaloJetCollection jetHMVtx ; //vettore jet appartenenti al vertice primario
    //loop sulla collezione di jet da cui seleziono quelli appartenenti al vertice con massima molteplicita' 
    //e li metto in un vettore di simplecalojet 

    SimpleCaloJetCollection::const_iterator vec_CaloJetIt = vec_CaloJet.begin(); 
    for(; vec_CaloJetIt != vec_CaloJet.end(); vec_CaloJetIt++){
      std::vector< int >::const_iterator   vec_idJetRecoHMProbEvalIt = vec_idJetRecoHMProbEval.begin();
      for(; vec_idJetRecoHMProbEvalIt != vec_idJetRecoHMProbEval.end(); vec_idJetRecoHMProbEvalIt++){
	// se il jet appartiene al vertice HM 
	if(vec_CaloJetIt->idJet() == *vec_idJetRecoHMProbEvalIt) jetHMVtx.push_back( SimpleCaloJet(*vec_CaloJetIt ) ) ;
      }
    }

    //insert the jet with no tracks
    SimpleCaloJetCollection::const_iterator vec_CaloJetWithNoTracksIt = vec_CaloJetWithNoTracks.begin();
    for(; vec_CaloJetWithNoTracksIt != vec_CaloJetWithNoTracks.end(); vec_CaloJetWithNoTracksIt++){
      jetHMVtx.push_back( SimpleCaloJet(*vec_CaloJetWithNoTracksIt ) ) ;
    }
    
    
    //ordino il vettore di calojet per valori decrescenti di Et
    sort(jetHMVtx.rbegin(),jetHMVtx.rend(),EtSort);
    
    //se ho 4 o piu' jet assegnati riempio le distribuzioni di Et per 4, 5 e 6 jet
    if(jetHMVtx.size()>3){
      H_Et4jetProbEval_ ->Fill(jetHMVtx[3].et());
      if(jetHMVtx.size()>4){
	H_Et5jetProbEval_ ->Fill(jetHMVtx[4].et());
	if(jetHMVtx.size()>5){
	  H_Et6jetProbEval_ ->Fill(jetHMVtx[5].et());
	}else{
	  H_Et6jetProbEval_->Fill(0.);
	}
      }else{
	H_Et5jetProbEval_->Fill(0.);
	H_Et6jetProbEval_->Fill(0.);
      }
    }else{
      H_Et4jetProbEval_->Fill(0.);
      H_Et5jetProbEval_->Fill(0.);
      H_Et6jetProbEval_->Fill(0.);
    }
 
    //    cout<<jetHMVtx[0].et()<<", "<<jetHMVtx[1].et()<<", "<<jetHMVtx[2].et()<<endl;
    
    sort( vec_RowCaloJet.rbegin(), vec_RowCaloJet.rend(), EtSort );
    
    if(vec_RowCaloJet.size()>3){
      H_Et4jetBeforeAlg_ ->Fill(vec_RowCaloJet[3].et());
      if(vec_RowCaloJet.size()>4){
	H_Et5jetBeforeAlg_ ->Fill(vec_RowCaloJet[4].et());
	if(vec_RowCaloJet.size()>5){
	  H_Et6jetBeforeAlg_->Fill(vec_RowCaloJet[5].et());
	}else{
	  H_Et6jetBeforeAlg_->Fill(0.);
	}
      }else{
	H_Et5jetBeforeAlg_->Fill(0.);
	H_Et6jetBeforeAlg_->Fill(0.);
      }
    }else{
      H_Et4jetBeforeAlg_->Fill(0.);
      H_Et5jetBeforeAlg_->Fill(0.);
      H_Et6jetBeforeAlg_->Fill(0.);

    }



    //    cout<<vec_RowCaloJet[0].et()<<", "<<vec_RowCaloJet[1].et()<<", "<<vec_RowCaloJet[2].et()<<endl;
        
  } //end if !caloJet.empty()  

  //=====ROBERTO=================================
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void EfficencyEval::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void EfficencyEval::endJob() {
  
  
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
  //  H_idVtxSAAlg_->Write();
  
  //Histograms of difference of log of probability of best vertex association and second best vs 
  //difference of log of probability considering the vertexes jet multeplicity
  //  H_pVsp_->Write();
  //Histograms of # of vertexes with a specified jet multeplicity
  //  H_nVtxMap_->Write();
  H_nVtxMapRJAlg_->Write();
  //  H_nVtxMapSAAlg_->Write();
  H_nVtxMapProbEval_->Write();
  
  //Histograms of probability  of jet assigned to vertexes configuration
  //  H_ProbRJAlg_->Write();
  //  H_ProbAfterFAAlg_->Write();
  //  H_ProbAfterSAAlg_->Write();

  //Histograms of id vertex assigned to a jet after first association
  //  H_idVtxFAAlg_->Write();
  //Histograms of Dz jet - vertex after SAAlg
  //  H_DzAfterSAAlg_->Write();

  //Histograms of # Tracks for jet with #DeltaZ best
  H_nTrack_RJAlg_->Write();

  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after SAAlg
  //  H_nJMaxVsNjTot_SAAlg_->Write();
  //Histograms of # jet for vertex with max multeplicity  vs # jet tot after RJAlg
  H_nJMaxVsNjTot_RJAlg_->Write();

  //Histograms of # Tracks for jet with #DeltaZ best < 5#sigma vs z jet (error)
  H_nTkRJAlgVSzJetError_->Write();
  H_nTkRJAlgVSzJet_->Write();

  //Histograms of #jet with and without tracks (only cinematics cuts)
  H_nJet_->Write();

  //Histograms of #jet before, after RJAlg and SAAlg (with tracks)
  //  H_nJMax_SAAlg_->Write();
  H_nJMax_RJAlg_->Write();
  H_nJMax_BeforeAlg_->Write();
  H_nJMax_ProbEval_->Write();
  
  //Histograms of Minimum Dz between vertexes
  //  H_VtxMinDz_->Write();

  H_ProfileVtxVsNjBeforeAlg_->Write();
  H_ProfileVtxVsNjRJAlg_->Write();
  //  H_ProfileVtxVsNjSAAlg_->Write();
  H_ProfileVtxVsNjProbEval_->Write();

  H_ZVtx_->Write();
  H_ZVtxError_->Write();

  H_DzAfterProbEval_->Write();

  H_ZBeforevsAfterRJAlg_->Write();

  H_Et4jetProbEval_->Write();
  H_Et5jetProbEval_->Write();
  H_Et6jetProbEval_->Write();
  
  H_Et4jetBeforeAlg_->Write();
  H_Et5jetBeforeAlg_->Write();
  H_Et6jetBeforeAlg_->Write();
  


  //   H_ZPvVsZBeforeAlg_->Write();
  //   H_ZPvVsZRJAlg_->Write();
  //   H_ZPvVsZProbEval_->Write();

  //=====ROBERTO=============
  
}


// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(EfficencyEval);
