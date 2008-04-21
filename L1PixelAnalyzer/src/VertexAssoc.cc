//////////////////////////////////////////////////////////////////////////////
//
// VertexAssoc.cc
// Code extracted from Roberto.cc, 02/04/08 R.Casagrande
//
//  
//
// --------------------------------------------------------------------------------
//
// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/VertexAssoc.h"

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
#include "AnalysisExamples/L1PixelAnalyzer/interface/DzAssociator.h"
#include "AnalysisExamples/L1PixelAnalyzer/interface/D2Associator.h"
#include "AnalysisExamples/L1PixelAnalyzer/interface/WAverager.h"

//Implemented Classes
//------------------------
#include "AnalysisExamples/L1PixelAnalyzer/interface/RejectionAlg.h"

//======ROBERTO========================
// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------


// Constructors and destructor
// ---------------------------
VertexAssoc::VertexAssoc(const edm::ParameterSet& iConfig) :
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
  
  //----------------------------------
  // Reco Tracks and Vertices
  //----------------------------------  

  //Reco tracks impact parameter
  impactParameterTagInfos( iConfig.getUntrackedParameter<std::string>( "impactParameterTagInfos" ) ),

  
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

  //Histograms for primary vertex and PU vertices
  H_NPUVtx_ = new TH1D ("H_NPUVtx","# of pileup vertices of events",50,0.,50.);

  //Histograms for Z PV vertices
  H_primaryVtxZ_ =  new TH1D ("H_primaryVtxZ","Z of primary vertex of events",40,-20.,20.);
  //Histograms for Z PU vertices
  H_PUVtxZ_ = new TH1D ("H_PUVtxZ","Z of pileup vertices of event",40,-20.,20.);
  //Histograms for Z simTracks
  H_Z_simTks_ =  new TH1D("H_Z_simTks", "Z of sim tracks" ,310,-310.,310.);
  //Histograms for Z reco tracks
  H_TksZ_  =  new TH1D("H_TksZ", "Z of tracks" ,200,-20.,20.);
  H_TksZ_Error_  =  new TH1D("H_TksZ_Error", "Error Z of tracks" ,200,0.,0.2);

  //Histogram Z jet 
  H_TksZ_Weighted_Avg_ = new TH1D("H_TksZ_Weighted_Avg", "Z jet (W Average Z of tracks)" ,40,-20.,20.);
  H_TksZ_WAvg_Error_ = new TH1D("H_TksZ_WAvg_Error", "Error Z jet (W Average Z of tracks)" ,200,0.,0.2);

  //Histograms selected tracks  
  H_SelTk_ =  new TH1D("H_SelTk", "# of selected tracks assigned to a jet " ,40,0.,40.);

  //Histograms for  2D, 3D IP Significance
  H_S2DTks_ = new TH1D("H_S2DTks", "Significance of tracks" ,300,-15.,15.);
  H_S3DTks_ = new TH1D("H_S3DTks", "Significance of tracks" ,300,-15.,15.);
  H_Tks_Z_S3D_ = new TH2D("H_Tks_Z_S3D","z PV vs S3D of tracks" ,20,-20.,20.,20,-20.,20.);
  
  //Histograms for Njet   
  H_NJet_withTracks_ =  new TH1D("H_NJet_withTracks", "# of jet (Et>25) with assigned tracks " ,40,0.,40.);
  H_NJet_withNoTracks_ =  new TH1D("H_NJet_withNoTracks", "# of jet (Et>25) without assigned tracks " ,20,0.,20.);
  H_NJet_fromPV_ = new TH1D("H_NJet_fromPV", "# of jet (Et>25) from Primary Vertex" ,40,0.,40.);
  H_NJet_fromPU_ = new TH1D("H_NJet_fromPU", "# of jet (Et>25) from PU Vertices " ,40,0.,40.);
  H_NJet_NotAssigned_ = new TH1D("H_NJet_NotAssigned", "# of jet (Et>25) Not assigned to any vertices" ,40,0.,40.);

  //Histograms for Dz jet - vertices 
  H_Vmin_ =  new TH1D("H_Vmin", "vertex # assigned to a jet (W Avg) " ,50,0.,50.);
  H_Dz_WAvg_Jet_Vtx_tot_ =  new TH1D("H_Dz_Jet_Vtx_tot", "Delta Z jet - vertex (W Avg) " ,80,-4.,4.);
  H_Dz_WAvg_Jet_Vtx_2Best_tot_ =  new TH1D("H_Dz_Jet_Vtx_2Best_tot", "Delta Z second best jet - vertex(W Avg) " ,80,-4.,4.);
  H_Dz_WAvg_Jet_Vtx_vs_2Best_ =  new TH2D("H_Dz_Jet_Vtx_vs_2Best","DeltaZ vs DeltaZ second best jet - vertex(W Avg)" ,80,0.,4.,80,0.,4.);

  //Histograms for Dz vertices
  H_Dz_Vtx_ = new TH1D("H_Dz_Vtx", "Delta Z vertici " ,40,0.,4.);
  
  //Histogram for Dz jet - vertices 
  H_Dz_WAvg_Vtx_vs_Nj_ = new TProfile("H_Dz_WAvg_Vtx_vs_Nj ","Difference Z vtx vs Nj (WAvg)",80,-20.,20.);
  
  H_Z_WAvg_Vtx_vs_Nj_ = new TProfile("H_Z_WAvg_Vtx_vs_Nj","Z vtx vs Nj (WAvg)",80,-20.,20.);
    
  //Histogram (phases space) for Dz jet - PU vertices
  H_Profile_WAvg_DzVtx_ = new TProfile("H_Profile_WAvg_DzVtx","Dz vtx PU vs Nj (WAvg)",80,-20.,20.);
  
  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
  const int ntracks = 13;
  char title_1[30];
  for ( int i=0; i<ntracks; i++ ) {
    sprintf (title_1,"H_Z_WAvg_Jet_%dtracks", i+1 );
    H_Z_WAvg_Jet_Ntracks_[i]= new TH1D ( title_1, title_1, 100, -5., 5.);
  }
  
  //Histograms for eta, phi reco tracks
  H_TksEta_ =   new TH1D("H_TksEta", "Eta of recoTracks" ,40,-10.,10.);
  H_TksPhi_ =   new TH1D("H_TksPhi", "Phi of recoTracks" ,20,-5.,5.);
  //Histograms for Pt, Chi2 reco tracks
  H_TksPt_  =  new TH1D("H_TksPt", "Pt of tracks" ,80,0.,40.);
  H_TksChi2_ =   new TH1D("H_TksChi2", "Chi2 of tracks" ,20,0.,10.);
  
  //Histograms for Et, eta, phi of calojets
  H_Et_Jets_= new TH1D("H_Et_Jets ", "Et of jets" ,300,0.,600.);
  H_eta_Jets_= new TH1D("H_eta_Jets", "Eta of jets" ,40,-10.,10.);
  H_phi_Jets_= new TH1D("H_phi_Jets", "Phi of jets" ,20,-5.,5.);
  
  //Histograms Dphi, Deta, Delta2 jet - reco selected tracks
  H_Dphi_Jet_Tk_ = new TH1D("H_Dphi_Jet_Tk", "Dphi jet - tracks" ,20,-1.,1.);
  H_Deta_Jet_Tk_ = new TH1D("H_Deta_Jet_Tk", "Deta jet - tracks" ,20,-1.,1.);
  H_D2_Jet_Tk_ = new TH1D("H_D2_Jet_Tk", "D2 jet - tracks" ,80,0.,0.4);

  //Histograms for eta, phi sim tracks
  H_Phi_simTk_ = new TH1D("H_Phi_simTk", "Phi of simTracks" ,20,-5.,5.);
  H_Eta_simTk_ = new TH1D("H_Eta_simTk", "Eta of simTracks" ,40,-10.,10.);
  
  //Histograms for Delta2 minimum simTrack - reco selected track
  H_D2_Tk_sim_reco_= new TH1D("H_D2_Tk_sim_reco", "D2 simTracks - reco selected tracks" ,200,0.,0.0002);    
  H_D2_Tk_sim_reco_2Best_= new TH1D("H_D2_Tk_sim_reco_2Best", "D2 simTracks - reco selected tracks 2 Best" ,200,0.,0.045);
  H_D2_Tk_vs_2Best_= new TH2D("H_D2_Tk_vs_2Best", "D2 simTracks - reco selected tracks vs 2 Best " ,200,0.,0.045,200,0.,0.045);  
  
  //histograms of Dz simTrack - reco selected track
  H_Dz_Tk_sim_reco_=  new TH1D("H_Dz_Tk_sim_reco", "Delta Z sim track - reco selected track " ,400,-2.,2.);
  H_Z_sim_jet_ = new TH1D("H_Z_sim_jet", "Z jet sim  " ,80,-20.,20.);
  H_Dz_jet_sim_reco_ = new TH1D("H_Dz_jet_sim_reco", "Delta Z jet sim - reco " ,200,-10.,10.);
  H_Dz_jet_sim_reco_avg_ = new TH1D("H_Dz_jet_sim_reco_avg", "Delta Z jet sim - reco avg " ,200,-10.,10.);

  //Histogram of Dz selected sim tracks
  H_Dz_sim_selected_tk_ = new TH1D("H_Dz_sim_selected_tk", "Delta Z sim selected tracks " ,800,-20.,20.);
  //TProfile max rate of same Z in selected sim tracks in function of the # of tracks
  H_Profile_Ntk_vs_sameZ_ = new TProfile("H_Profile_Ntk_vs_sameZ","N tracks vs N tracks with same Z",20, 0.,20.);

  //Histograms of rate of same selected sim tracks Z and 
  //rate of same selected sim tracks Z vs # of tracks in the jet
  H_Frequenze_ =  new TH1D("H_Frequenze", "Frequenze sim selected tracks di = Z " ,40,-20,20);
  H_Ntk_vs_SameZ_ =  new TH2D("H_Ntk_vs_SameZ", "# tracks vs # same Z " ,20,0.,20,20,0.,20);
   
  //Histogram of selected simulated tracks' Z
  H_Z_selected_sim_tracks_ =  new TH1D("H_Z_selected_sim_tk", "Z sim selected tracks " ,200,-20.,20.);
  
  //Histograms of Dz jet - vertex: |Z(true) - Z(vtx)| 
  H_Dz_sim_Jet_Vtx_ = new TH1D("H_Dz_sim_Jet_Vtx", "Delta Z jet - vertex (sim) " ,80,-4.,4.);
  H_Dz_sim_Jet_Vtx_2Best_ =  new TH1D("H_Dz_sim_Jet_Vtx_2Best", "Delta Z jet - vertex (sim) " ,80,-4.,4.);
  H_Dz_sim_Jet_Vtx_vs_2Best_ =  new TH2D("H_Dz_sim_Jet_Vtx_vs_2Best","DeltaZ vs DeltaZ second best jet - vertex (sim)" ,40,0.,4.,40,0.,4.);
  
  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
  char title[30];
  for ( int i=0; i<ntracks; i++ ) {
    sprintf (title,"H_Dz_sim_Jet_Vtx_%dTk", i+1 );
    H_Dz_sim_Jet_Vtx_NTk_[i]= new TH1D ( title, title, 100, -5., 5.);
  }
  
  //Histogram of |Z(vera) - Z(reco)| jet
  H_Dz_jet_true_reco_ = new TH1D(" H_Dz_jet_true_reco", "Delta Z= |Z(vera) - Z(reco)| " ,80,-4.,4.);
  
  //Histograms of Dz vera jet - Z sim vertex associated
  H_Dz_sim_jet_simVtx_ = new TH1D("H_Dz_sim_jet_simVtx", "Delta Z sim jet - simvertex " ,80,-4.,4.);    
  H_Dz_sim_jet_simVtx_2Best_ =  new TH1D("H_Dz_sim_jet_simVtx_2Best", "Delta Z sim jet - simvertex " ,80,-0.4,0.4);
  H_Dz_sim_jet_simVtx_vs_2Best_ =  new TH2D("H_Dz_sim_jet_simVtx_vs_2Best","DeltaZ vs DeltaZ second best sim jet - simvertex" ,40,0.,4.,40,0.,4.);
  
  //Histograms of Z jet and error after the rejection algorithm (5 sigma)
  H_ZJet_RJAlg_ = new TH1D("H_ZJet_RJAlg", "Z jet (W Average Z of tracks after rejection alg)" ,40,-20.,20.);
  H_ZErrorJet_RJAlg_ = new TH1D("H_ZErrorJet_RJAlg", "Z jet (W Average Z of tracks after rejection alg)" ,200,0.,0.2);
  H_DzAfterRejection_ =  new TH1D("H_DzAfterRejection", "Delta Z jet after rejection - vertex " ,400,-2.,2.);
 
  //Histograms of dZ:  Z jet( weight average after rejection of i sigma ) - vertices 

  char title_2[30];
  for ( int i=0; i<20; i++ ) {
    sprintf (title_2,"H_dZJet_RJAlg_%dSigma", i );
    H_dZJet_RJAlg_SFunc_[i]= new TH1D ( title_2, title_2, 400, -0.2, 0.2);
  }
  
  //Histogram of Ntk in function # sigma rejection
  H_Ntk_vs_NS_= new TProfile("H_Ntk_vs_NS","# tracks vs # sigma rejection" ,20,0.,20.);
  
  
  H_NJetRJAlgFromPV_ =  new TH1D("H_NJetRJAlgFromPV", "# of jet (Et>25) from Primary Vertex" ,40,0.,40.);
  H_NJetRJAlgFromPU_ =  new TH1D("H_NJetRJAlgFromPU", "# of jet (Et>25) from PU Vertex" ,40,0.,40.);


  H_ProfileDzVtxVsNjRJAlg_ = new TProfile("H_ProfileDzVtxVsNjRJAlg","Difference Z vtx vs Nj after Rejection Algorithm",80,-20.,20.);
  
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


VertexAssoc::~VertexAssoc()
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
void VertexAssoc::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  // General definitions
  // -------------------
  //   double mhref=117.; 
  //   double mtref=186.;   // Massa top ricostruita, >172 per via dello shift in JEcorrs.
  //   double mwref=86.;

  
  //======ROBERTO=======================
  
  int nAccettanza = 20;

  //------------------------------
  // Pile-up vertices collection
  //------------------------------
  
  Handle<HepMCProduct> theFamosPileUp;
  bool famosPU = false;
  const HepMC::GenEvent* myGenPUEvent = 0;

  // Look for the GenEvent
  std::vector< Handle<HepMCProduct> > evts; 
  iEvent.getManyByType(evts);

  for ( unsigned i=0; i<evts.size(); ++i ) {
    //prendo solo le info sui PU, nella source ci sono i vertici primari non smeared, i vertici primari li prendo sotto dalla simvtx
    if (!famosPU &&  evts[i].provenance()->moduleLabel()=="famosPileUp") {
      famosPU = true;
      theFamosPileUp = evts[i];
    }
  }
   
  // cout<<"type= "<<theFamosPileUp.provenance()->moduleLabel()<<endl; 
 
  //-------------------------------
  // Simulated vertices collection
  //-------------------------------

  // Take the SimVertex collection
  Handle<SimVertexContainer> SimVtx;
  try{
    iEvent.getByLabel( simVtxLabel_, SimVtx );
    //  cout << "SimVtx size = " << SimVtx->size() << endl;
  }
  catch(...){
    std::cerr << "Could not find the SimVertex collection" << endl;
  }
  
  //Simvertex per associazione tracce simulate al vertice simulato
  std::vector<SimVertex> theSimVertexes(*(SimVtx.product()));
  
  //Equivalent methods
  //  std::vector<SimVertex> theSimVertexes(*SimVtx);
  //  theSimVertexes.insert(theSimVertexes.end(),SimVtx->begin(),SimVtx->end());

  //------------------------------  
  // Simulated tracks collection 
  //------------------------------

  //simtrack per l'associazione tracce simulate vertici simulati 
  // Take the SimTracks collection
  
  Handle<SimTrackContainer> SimTk;
  
  try {
    iEvent.getByLabel( simTkLabel_, SimTk );
    //  cout << "SimTk size = " << SimTk->size() << endl;
  }
  catch (...) {
    std::cerr << "Could not find the SimTracks collection" << std::endl;
    return;
  } 
  
  std::vector<SimTrack> theSimTracks(*(SimTk.product()));

  //inserisco prima della fine del vettore di SimTrack le tracce simulate
  // theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());

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
  
  //======ROBERTO=======================
  
  
  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }
  
  //======ROBERTO================================
  
  ////////////////////////////////// MC analysis ////////////////////////////

  //------------------
  // Pile-up vertices
  //------------------
 
  myGenPUEvent = theFamosPileUp->GetEvent();

  HepMC::GenEvent::vertex_const_iterator vPUiter;
  HepMC::GenEvent::vertex_const_iterator vPUbegin = myGenPUEvent->vertices_begin();
  HepMC::GenEvent::vertex_const_iterator vPUend = myGenPUEvent->vertices_end();
  
  int iPUevt = 0;
  
  //vettore di vertici di pile up dell'evento
  std::vector< math::XYZTLorentzVector > vec_PUVertex;
  //vettore di vertici, in cui il primo e' il vertice primario
  std::vector< math::XYZTLorentzVector > vec_Vertices;
  
  // Loop on all pile-up events
  for ( vPUiter=vPUbegin; vPUiter!=vPUend; ++vPUiter ) { 
    
    // std::cout << "Vertex n0 " << iPUevt << std::endl;
    
    // The origin vertex (turn it to cm's from GenEvent mm's)
    HepMC::GenVertex* vPU = *vPUiter;    
    math::XYZTLorentzVector PileUpVertex =  
      math::XYZTLorentzVector(vPU->position().x()/10.,vPU->position().y()/10.,
			      vPU->position().z()/10.,vPU->position().t()/10.);
    
    //riempio un vettore con i vertici di pileup dell'evento
    vec_PUVertex.push_back(PileUpVertex);

    //Histograms for Z PU vertices   
    H_PUVtxZ_->Fill( vPU->position().z()/10.);

    ++iPUevt;
  }
  
  H_NPUVtx_->Fill(iPUevt);
 
  // cout << "PU vertex number = " << iPUevt << endl;
  // cout << "vec_PUVertex size = " <<vec_PUVertex.size()<<endl;

  //------------------
  // Simulated vertices
  //------------------ 
  //verifica sempre true il primo dei simvertex
  //perche' assumiamo che il primo sia il vertice primario
  int c0 = -1;
  int c1 = -1;
  std::vector<SimVertex>::const_iterator isimvtx = theSimVertexes.begin();
  for(; isimvtx != theSimVertexes.end();isimvtx++){
    if((*isimvtx).noParent() == 1) c1++;
    if((*isimvtx).noParent() == 0) c0++;
  }

  //   cout << "SimVtx = " << (*SimVtx->begin()).noParent() << endl;
  //   cout<<"for primary vertex= "<<(*theSimVertexes.begin()).noParent()<<endl;
  //   cout<<"false= "<<c0<<", true= "<<c1<<endl;
  //   cout<<"theSimVertexes size= "<<theSimVertexes.size()<<endl;
  
  //----------------
  // Primary vertex
  //----------------

  //metto via un vettore col vertice primario dell'evento
  math::XYZTLorentzVector primarySimVertex =  
    math::XYZTLorentzVector(SimVtx->begin()->position().x() ,
			    SimVtx->begin()->position().y() ,
			    SimVtx->begin()->position().z() ,
			    SimVtx->begin()->position().t()
			    );

  //Histograms for Z PV vertices
  H_primaryVtxZ_->Fill( SimVtx->begin()->position().z());
  
  //riempio un vettore con il vertice primario e di seguito i vertici di pile up 
  //per poi fare il ciclo di associazione ai jet 
  vec_Vertices.push_back(primarySimVertex);
  vector<math::XYZTLorentzVector>::const_iterator PUVertex_it = vec_PUVertex.begin();
  for(; PUVertex_it != vec_PUVertex.end(); PUVertex_it++ ){
    // need *iterator to work
    vec_Vertices.push_back(*PUVertex_it);
  }
  
  //  cout << "Vertices size = " << vec_Vertices.size() << endl;

  
  //--------------------------
  //Simulated tracks vertices
  //-------------------------- 
  //faccio un loop sul vettore di tracce simulate per estrarre le informazioni sulle Z della traccia
  //Come da esempio in AnalysisExamples/SimTrackerAnalysis/src/SimTrackSimVertexAnalyzer.cc
  
  //Vector of simTracks' Z, eta, phi  
  SimpleTrackCollection vec_SimTk;

  double theSimVertexes_size = 0.;
  theSimVertexes_size =  theSimVertexes.size();

  for (int isimvtx = 0; isimvtx < theSimVertexes_size;isimvtx++){  //, simvtx_it++){  
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
	
	Z_simTk = theSimVertexes[isimvtx].position().z();

	//non considero le tracce senza vertice associato
	vec_SimTk.push_back(SimpleTrack(Pt_simTk, Eta_simTk, Phi_simTk, Z_simTk, 0.));

	//Histograms for Z simTracks
	H_Z_simTks_->Fill(Z_simTk);
	
	//Histograms for eta, phi sim tracks 
	H_Phi_simTk_->Fill(Phi_simTk);
	H_Eta_simTk_->Fill(Eta_simTk);


      }//end if track associated to this vertex
    }//end loop on sim tracks
  }//end loop on sim vertexes
  
  //--------------------------------
  //Reco tracks and Offline calojet (classe per simple calojet +z error e significanza traccie)
  //--------------------------------  
   
  //# of caloJets
  int caloJetsSize = 0;
  caloJetsSize = caloJets->size();
  //  cout<<"# calojet: "<<caloJetsSize<<endl;
 
  //vector of caloJets' Z
  std::vector<double> vec_Z_Jet;
  std::vector<double> vec_Z_error_Jet;

  //Vettori per storage informazioni sulla z vertici piu' vicina al jet
  //ordinati come i jet,con Vmin che dice quale e' il vertice piu' vicino
  std::vector<int> vec_Vmin;
  std::vector<int> vec_idRJAlg;

//   std::vector<double> vec_Dz_WAvg_Jet_Vtx_;
   
  //Et, eta, phi Offline caloJet
  double Et_Jet = 0.;
  double eta_Jet = 0.;
  double phi_Jet = 0.;

  //vectors of Et, eta, phi Offline caloJet
  std::vector<double> vec_Et_Jet; 
  std::vector<double> vec_eta_Jet;
  std::vector<double> vec_phi_Jet;

  //Come si fa a fare un vettore di vettori di double?
  //std::vector<std::vector(double)> vec_eta_Tk_vectors;
  
  int Njet = 0 ; //# jet
  int NJet_withTracks = 0; //# jet with assigned tracks
  int NJet_withNoTracks = 0; //# jet without assigned tracks
  int NJet_fromPV = 0; //#jet assigned to Primary Vertex
  int NJet_fromPU = 0; //#jet assigned to PU vertices
  int NJet_NotAssigned = 0; //#jet (Et >25, # of tracks with S<3  >0) not assigned to any vertex
  int Ntk = 0; //# tracks

  int NJetRJAlgFromPV = 0;
  int NJetRJAlgFromPU = 0;


  //loop on tracks (IptagInfo) and caloJets 
  anaobj::OfflineJetCollection::const_iterator caloJets_it = caloJets->begin();
  reco::TrackIPTagInfoCollection::const_iterator TkIpTagInfo_it = iPtagInfos->begin();
  for( ; TkIpTagInfo_it != iPtagInfos->end(); ++TkIpTagInfo_it ) {
   
    
    //vector of reco selected tracks' significance 2D and 3D of IP
    std::vector<double>  vec_Tk_IpS2D; //vector of reco selected tracks' significance 2D of IP
    std::vector<double>  vec_Tk_IpS3D; //vector of reco selected tracks' significance  3D of IP
    //    std::vector<double> vec_Tk_dz_error;  //vector of reco selected tracks' dz error
    std::vector<double> vec_Tk_Chi2; //vector of reco selected tracks' chi2
    //    std::vector<double> sorted_vec_Z_selected_sim_tracks;//vector of sorted selected simulated tracks' Z

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
    
    //# of selected tracks
    int SelTk_size = 0; 
    SelTk_size = vec_TkColl.size();
    // cout<<"# Selected Tracks= "<<SelTk_size<<endl;
      
    //Histogram for # of tracks assigned to a single jet without S<3 cut
    //Histograms selected tracks  
    H_SelTk_->Fill(SelTk_size);

    //# of selected tracks counter
    int STk_count = 0;
    
    //loop on selected tracks collection and tracks IP  
    vector<reco::TrackIPTagInfo::TrackIPData>::const_iterator TkIP_it = vec_TkIP.begin();
    RefVector<reco::TrackCollection>::const_iterator TkColl_it = vec_TkColl.begin();
    for ( ; TkColl_it != vec_TkColl.end(); ++TkColl_it, ++TkIP_it ) {
      
      //3D Significance
      double Tk_IpS3D = 0.; 
      Tk_IpS3D = TkIP_it->ip3d.significance();
      //      cout<<"significanza 3D IP traccia: "<<Tk_IpS3D<<endl;
      vec_Tk_IpS3D.push_back(Tk_IpS3D);
      
      //2D Significance
      double Tk_IpS2D = 0.;
      Tk_IpS2D =TkIP_it->ip2d.significance();
      //      cout<<"significanza 2D IP traccia: "<<Tk_IpS2D<<endl;
      vec_Tk_IpS2D.push_back(Tk_IpS2D);

      //uso dz invece di vz, ha l'errore sul parametro e e' un'approssimazione lineare buona
      //-------DataFormats/TrackReco/interface/TrackBase.h-----------------------------------
      // dz parameter with respect to a user-given beamSpot (WARNING: this quantity can only 
      // be interpreted as the track z0, if the beamSpot is reasonably close to the refPoint, since linear 
      // approximations are involved). This is a good approximation for Tracker tracks.
      //--------------------------------------------------------------------------------------
      
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
    
      //Reco Track Chi2
      double Tk_Chi2 = 0.;
      Tk_Chi2 = (*TkColl_it)->normalizedChi2();
      // cout<<"Chi2 traccia: "<<Tk_Chi2<<endl;
    
      //cut
      //-------

      //Fill vector with info if  Et(jet)>25, IP significance < 3. , eta (jet) < 3.   
      if(Et_Jet>25 && fabs(Tk_IpS2D) < 3.  && eta_Jet < 3){
	//	vec_Tk_dz_error.push_back(Tk_dz_error);
	vec_Tk_Chi2.push_back(Tk_Chi2);
	vec_RecoTk.push_back(SimpleTrack(Tk_Pt, Tk_eta, Tk_phi, Tk_dz, Tk_dz_error));
      }//end cut Et(jet)>25, IP significance < 3.

      //---------------------------------
      // Delta2 jet - track
      // (Verifing jet - selected tracks  (ok)
      // association)
      //---------------------------------
      
      double D2_Jet_Tk = 0.;
      double Dphi_Jet_Tk = 0.;
      double Deta_Jet_Tk = 0.;
      
      //Dphi jet - track
      Dphi_Jet_Tk  = M_PI - fabs(fabs(phi_Jet - Tk_phi) - M_PI);
      // cout<< Dphi_Jet_Tk <<endl;
     
      //Deta jet - track
      Deta_Jet_Tk = eta_Jet -  Tk_eta;
      // cout<< Deta_Jet_Tk <<endl;
    
      //Delta2 jet - track
      D2_Jet_Tk = pow(Dphi_Jet_Tk,2)+pow(Deta_Jet_Tk,2);
      // cout<<D2_Jet_Tk<<endl;
      
      //Fill histograms with info if  Et(jet)>25, IP significance < 3.,eta (jet) < 3.
      if(Et_Jet>25 && fabs(Tk_IpS2D) < 3. && eta_Jet < 3 ) {
	//Histograms for Z reco tracks
	H_TksZ_->Fill(Tk_dz);
	H_TksZ_Error_->Fill(Tk_dz_error);
	//Histograms for eta, phi reco tracks
	H_TksEta_->Fill(Tk_eta);
	H_TksPhi_->Fill(Tk_phi);
	//Histograms for Pt, Chi2 reco tracks
	H_TksPt_->Fill(Tk_Pt);
	H_TksChi2_->Fill(Tk_Chi2);
	//Histograms 2D, 3D IP Significance
	H_S2DTks_->Fill(Tk_IpS2D);
	H_S3DTks_->Fill(Tk_IpS3D);
	H_Tks_Z_S3D_->Fill(SimVtx->begin()->position().z() ,Tk_IpS3D);
	//Histograms Dphi, Deta, Delta2 jet - reco selected tracks
	H_Dphi_Jet_Tk_->Fill(Dphi_Jet_Tk);
	H_Deta_Jet_Tk_->Fill(Deta_Jet_Tk);
	H_D2_Jet_Tk_->Fill(D2_Jet_Tk);
      }//end cut Et(jet)>25, IP significance < 3.
      
      ++STk_count;

    }// end loop on selected tracks
    
 
    //    cout<<vec_Tk_dz_error.size()<<endl;  
    //    cout<<vec_Tk_Chi2.size()<<endl;
    //    cout<<"# track processed= "<<STk_count<<endl;

    //------coordinata Z del jet---------------------------------------------------
    //calcolo il vertice medio delle tracce per questo gruppo di selected tracks
    //e poi metto i vertici cosi' calcolati in un vettore per un confronto con i 
    //vertici di PU e primario definisco un iteratore che corre sul vettore di 
    //vertici delle tracce
    //-----------------------------------------------------------------------------

    //----------------------------------
    //Taglio Et(jet)>25 , eta(jet) < 3
    //----------------------------------
    
    if(Et_Jet>25 && eta_Jet < 3){
 
      //----------------------------------------------------
      // Jet Z: Wieghted Average on tracks' Z, with tracks (ok)
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


      if (NTks_S3 != 0){ 	
	vec_Z_Jet.push_back(wAvgRecoTkValue);	
	vec_Z_error_Jet.push_back(wAvgRecoTkError);
	
	//vectors of Et, eta, phi Offline caloJets
	vec_Et_Jet.push_back(Et_Jet); 
	vec_eta_Jet.push_back(eta_Jet);
	vec_phi_Jet.push_back(phi_Jet);
	
	//Histogram Z Jet
	H_TksZ_Weighted_Avg_->Fill(wAvgRecoTkValue);
	H_TksZ_WAvg_Error_->Fill(wAvgRecoTkError);
	
	//Histograms for Et, eta, phi of calojets
	H_Et_Jets_->Fill(Et_Jet);
	H_eta_Jets_->Fill(eta_Jet);
	H_phi_Jets_->Fill(phi_Jet);
	
	NJet_withTracks++;
      } 
      else {
	//	cout<<" 0 selected tracks "<<endl;
	NJet_withNoTracks++;
      }

      //---------------------------------
      // Delta Z minimo tra jet e vertici (ok)
      // con Z media pesata tracce
      //---------------------------------
      //Se il numero di tracce non e' 0
      if (NTks_S3 != 0){
	
	//Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con l amedia pesata e i vertici (PV e PU)	
	//la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
	//la funzione vuole in input un XYZTLorentz vector e un double
	//----------------------------------------------------------------------------------------------------------------------

	pair<pair<int, double>, pair<int, double> > dzFirstSecond( dzAssociator(vec_Vertices, wAvgRecoTkValue) );
	pair<int, double> * first = &(dzFirstSecond.first);
	pair<int, double> * second = &(dzFirstSecond.second);
        int Vmin = first->first;
	double Dz_Jet_Vtx = first->second;
	double Dz_Jet_Vtx_2Best = second->second;

	vec_Vmin.push_back(Vmin);

// 	vec_Dz_WAvg_Jet_Vtx_.push_back(Dz_Jet_Vtx);

	//incremento contatori per il conteggio di quanti jet sono assegnati 
 	if(Vmin == 0) NJet_fromPV++; //il vertice primario occupa il primo posto nel vettore di vertici
       	
 	if(Vmin > 0) NJet_fromPU++; //i vertici di PU occupano i posti nel vettore di vertici successivi al primo

 	if(Vmin < 0) NJet_NotAssigned++; //se non e' assegnato nel loop il Vmin non e' modificato (e' un check per l'algoritmo)

	//Histograms minimum Dz Jet - Vertices	
	H_Vmin_->Fill(Vmin);
        H_Dz_WAvg_Jet_Vtx_tot_->Fill(Dz_Jet_Vtx);
	H_Dz_WAvg_Jet_Vtx_2Best_tot_->Fill(Dz_Jet_Vtx_2Best);
	H_Dz_WAvg_Jet_Vtx_vs_2Best_->Fill(fabs(Dz_Jet_Vtx),fabs(Dz_Jet_Vtx_2Best));
	
	//---------------------------------------------------------------
	//Istogrammi della z media pesata delle tracce (cioe' il vertice del jet) 
	//per numero di tracce nel jet da 1 a 12 e maggiore di 12
	//---------------------------------------------------------------
	//riempio l'istogramma corrispondente al # di tracce nel jet con la delta Z 
	
	//Histograms Dz jet - vertices for exclusive # of tracks in the jet 
	if ( NTks_S3>0 && NTks_S3<13) {
	  H_Z_WAvg_Jet_Ntracks_[NTks_S3 - 1]->Fill(Dz_Jet_Vtx);
	} 
	//se il numero di tracce e' maggiore o uguale a ntracks-1 (=12) riempio l'ultimo istogramma
	if ( NTks_S3>=13 ) H_Z_WAvg_Jet_Ntracks_[12]->Fill(Dz_Jet_Vtx);

	//----------------------------------------------
	//Simtracks - reco selected tracks association 
	//Delta2 minimum simTrack - reco selected track (ok)
	//----------------------------------------------

	double Z_sim_jet = 0.; //true Z jet (avg)
	double sum_Z_simTk = 0.; //sum Z sim tracks
	double Dz_jet_sim_reco = 0.; //Delta Z sim - reco jet

	double sum_Z_recoTk = 0.; //sum Z reco tracks
	double Z_reco_jet = 0.; //reco Z jet (avg)
	double Dz_jet_sim_reco_avg =0.; //Delta Z sim reco jet calculated with avg


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
	  double d2SimTk = d2SimReco.first.second ;
	  double d2SimTk2Best = d2SimReco.second.second;

	  //fill vector collection of selected simulated tracks' Z 
	  vec_SimTk_selected.push_back( vec_SimTk[idSimTk ] );
	  
	  //Compute sum of sim tracks' Z to use for average
	  sum_Z_recoTk += vec_RecoTk_iter->z();
	  sum_Z_simTk += vec_SimTk[ idSimTk ].z();	  

 	  //Histogram of selected simulated tracks' Z
 	  H_Z_selected_sim_tracks_->Fill( vec_SimTk[ idSimTk ].z() );
	  //Histograms for Delta2 minimum simTrack - reco selected track
	  H_D2_Tk_sim_reco_->Fill( d2SimTk );
	  H_D2_Tk_sim_reco_2Best_->Fill( d2SimTk2Best );
	  H_D2_Tk_vs_2Best_->Fill( d2SimTk, d2SimTk2Best );
	  //Histograms of Dz simTrack - reco selected track
	  H_Dz_Tk_sim_reco_->Fill( vec_RecoTk_iter->z() - vec_SimTk[ idSimTk ].z() );

	}//end loop reco selected tracks 

	Z_sim_jet = sum_Z_simTk/ NTks_S3; //Avg Z sim tracks = Z jet true
	Dz_jet_sim_reco = Z_sim_jet - wAvgRecoTkValue;
  
	Z_reco_jet = sum_Z_recoTk/ NTks_S3; //Avg Z reco tracks = Z jet reco
	Dz_jet_sim_reco_avg =  Z_sim_jet - Z_reco_jet;

	//--------------------------------------------
	//Delta Z tra le tracce simulate selezionate (ok)
	//--------------------------------------------
	
	// 	cout<<"Sim  selected tracks= "<<vec_SimTk_selected.size()<<endl;
	// 	cout<<"Reco selected tracks= "<<vec_RecoTk.size()<<endl;


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
	
	//--------------------------------------------
	//Conteggio frequenze tracce simulate uguali (ok)
	//--------------------------------------------

	// Utilizzo una mappa per fare l'ordinamento decrescente e il conteggio dei vari rate, uso come 
	// key la z della traccia e come value il rate, che incremento ogni volta che trovo una key uguale
	// sfrutto il fatto che se la key che sto inserendo e' gia' inserita la mappa restituisce 
	// un bool false quando si tenta l'inserimento
	//---------------------------------------------------------------------------------------------------

	map<double, int> mapFreq;

       	anaobj::SimpleTrackCollection::const_iterator vec_SimTk_it = vec_SimTk_selected.begin();
	for(; vec_SimTk_it !=  vec_SimTk_selected.end(); vec_SimTk_it++ ){      
	  pair<map<double, int>::iterator, bool> mapInsert(mapFreq.insert( make_pair( (*vec_SimTk_it).z(), 1 ) ));
	  if ( !(mapInsert.second) ) {
	    //	    cout << "before mapInsert = " << (mapInsert.first)->second << endl;
	    ++((mapInsert.first)->second);
	    //	    cout << "after mapInsert = " << (mapInsert.first)->second << endl;
	  }
	}
	
	// Determino la frequenza massima e il valore corrispondente di Z inserendo come valore di partenza la prima 
	// coppia della mappa e iterando sulla mappa dal secondo valore in poi. Il primo valore e' la z, il secondo 
	// la frequenza.
	//------------------------------------------------------------------------------------------------------------ 


	double zMax = 0.;
	int freq = 0;
	if ( !(mapFreq.empty()) ) {
	  zMax = mapFreq.begin()->first;
	  freq = mapFreq.begin()->second;
	  map<double, int>::const_iterator mapFreqIt = mapFreq.begin();
	  ++mapFreqIt;
	  for ( ; mapFreqIt != mapFreq.end(); ++mapFreqIt ) {
	    if ( freq < mapFreqIt->second ) {
	      zMax = mapFreqIt->first;
	      freq = mapFreqIt->second;
	    }
	  }
	}

// 	  //Histograms of rate of same selected sim tracks Z and 
// 	  //rate of same selected sim tracks Z vs # of tracks in the jet
// 	  H_Frequenze_->Fill( mapFreqIt->second);
// 	  H_Ntk_vs_SameZ_->Fill(NTks_S3, mapFreqIt->second);

	//if ho almeno 2 tracce simulate nel jet che puntano nello stesso punto

	if( freq > 1){
	  //------------------------------------
	  //Associazione Z jet vera - Z vertici (ok)
	  //------------------------------------
	  // 	  cout << "freq = " << freq << endl;
	  // 	  cout << "zMax = " << zMax << endl;

	  //richiamo la funzione dzAssociator per eseguire l'associazione tra zMax (Z vera edl jet) e Z dei vertici (PV e PU)

	  pair<pair<int, double>, pair<int, double> > dzJetTrueVertex( dzAssociator( vec_Vertices, zMax ) );
	  pair<int, double> * firstJV = &(dzJetTrueVertex.first);
	  pair<int, double> * secondJV = &(dzJetTrueVertex.second);
	  //	  int Vmin_sim = firstJV->first;
	  double Dz_sim_Jet_Vtx = firstJV->second;
	  double Dz_sim_Jet_Vtx_2Best = secondJV->second;
	
	  
	  //Histograms of Dz jet - vertex: |Z(true) - Z(vtx)| 
	  H_Dz_sim_Jet_Vtx_->Fill(Dz_sim_Jet_Vtx);
	  H_Dz_sim_Jet_Vtx_2Best_->Fill(Dz_sim_Jet_Vtx_2Best);
	  H_Dz_sim_Jet_Vtx_vs_2Best_->Fill(fabs(Dz_sim_Jet_Vtx),fabs(Dz_sim_Jet_Vtx_2Best));
	  
	  //---------------------------------------------------------------
	  //Istogrammi della Z vera del jet - Z vertice associato 
	  //per numero di tracce nel jet da 1 a 12 e maggiore di 12
	  //---------------------------------------------------------------
	  //riempio l'istogramma corrispondente al # di tracce nel jet con la delta Z 
	  
	  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
	  if (freq > 0 && freq <13){
	    H_Dz_sim_Jet_Vtx_NTk_[freq - 1]->Fill(Dz_sim_Jet_Vtx);
	  } 
	  //se il numero di tracce e' maggiore o uguale a ntracks-1 (=12) riempio l'ultimo istogramma
	  if ( freq >=13 ) H_Dz_sim_Jet_Vtx_NTk_[12]->Fill(Dz_sim_Jet_Vtx);

	  //------------------------
	  //|Z(vera) - Z(reco)| jet (ok)
	  //------------------------

	  //Histogram of |Z(vera) - Z(reco)| jet
	  H_Dz_jet_true_reco_->Fill( zMax - wAvgRecoTkValue );

	  //--------------------------------------------------------------------
	  //assegnazione jet - simVertex - Delta Z(vera) jet - vertice simulato (overloded function dzAssoc)
	  //--------------------------------------------------------------------
	 
	  //La parte commentata da valori per Dz2Best, la funzione dzAssociator da sempre 0 anche per il 2Best,
	  //cosa che non mi aspetto ( BUG? ) --> no, nei simVertex ci sono ripetizioni => uso una mappa per il conto 
	  //la funzione e' overloded per ricevere in input anche vector<SimVertex>
	  //---------------------------------------------------------------------------------------------------
	  pair<pair<int, double>, pair<int, double> > dzSimJetSimVtx( dzAssociator( theSimVertexes, zMax ) );
	  pair<int, double> * firstSJSV = &(dzSimJetSimVtx.first);
	  pair<int, double> * secondSJSV = &(dzSimJetSimVtx.second);
	  
	  double  Dz_sim_jet_simVtx  = firstSJSV->second;
	  double  Dz_sim_jet_simVtx_2Best = secondSJSV->second;

	  
	  //Histograms of Dz zMaxFreq jet - Z sim vertex associated
	  //Mi aspetto i best siano esattamente 0!!
	  //----------------------------------------------------
	  H_Dz_sim_jet_simVtx_->Fill(Dz_sim_jet_simVtx);
	  H_Dz_sim_jet_simVtx_2Best_->Fill(Dz_sim_jet_simVtx_2Best);
	  H_Dz_sim_jet_simVtx_vs_2Best_->Fill(Dz_sim_jet_simVtx, Dz_sim_jet_simVtx_2Best);


	}//end if ho almeno 2 tracce simulate nel jet che puntano nello stesso punto
	  
	//TProfile max rate of same Z in selected sim tracks in function of the # of tracks 
	H_Profile_Ntk_vs_sameZ_->Fill(NTks_S3,freq);
	
	//Histograms of sim tracks' Z, sim - reco tracks' Dz 
	H_Z_sim_jet_->Fill(Z_sim_jet);
	H_Dz_jet_sim_reco_->Fill(Dz_jet_sim_reco);
	H_Dz_jet_sim_reco_avg_->Fill(Dz_jet_sim_reco_avg);

	//----------------------------------------------
	//Algoritmo di reiezione tracce ricostruite nel jet
	//Z jet calcolata con la media pesata, scartando 
	//le tracce piu' distanti di 5 sigma dal valore 
	//della media pesata con algoritmo iterativo
	//----------------------------------------------
	//chiamo l'algoritmo di reiezione (e' una classe con il metodo eval(# sigma taglio)) sulla collezione di tracce ricostruite
	RejectionAlg tkRejection(vec_RecoTk);
	SimpleTrackCollection recoTkClosest(tkRejection.eval(5));
	cout<<"before rejection= "<<vec_RecoTk.size()<<endl;
	cout<<"after rejection= "<<recoTkClosest.size()<<endl<<endl;
	//chiamo la funzione che fa la media pesata wAvenger, che restituisce una pair(# tracce, pair(value, error))
	pair< int, pair< double, double > > zJetRejectionAlg( wAverager(recoTkClosest));
	
	//riempio gli istogrammi di value e error (Z del jet calcolata e suo errore) 
	//Histograms of Z jet and error after the rejection algorithm (5 sigma)
	H_ZJet_RJAlg_->Fill(zJetRejectionAlg.second.first);
	H_ZErrorJet_RJAlg_->Fill(zJetRejectionAlg.second.second);

	//---------------------------------
	// Delta Z minimo tra jet e vertici
	// con Z media pesata tracce
	//---------------------------------
	//Se il numero di tracce non e' 0
	
	//Richiamo la funzione dzAssociator per minimizzare il dz tra Z jet calcolata con la media pesata dopo la reiezione e i vertici (PV e PU)	
	//la prima coppia contiene il l'indice e il valore minimo, la seconda coppia l'indice e il valore del second best
	//la funzione vuole in input un XYZTLorentz vector e un double
	//----------------------------------------------------------------------------------------------------------------------
	
	pair<pair<int, double>, pair<int, double> > dzAfterRejection( dzAssociator(vec_Vertices,zJetRejectionAlg.second.first ) );
        
	int  idRJAlg = dzAfterRejection.first.first;
	double dzRJAlg = dzAfterRejection.first.second;
	//	double dzRJAlg2Best = dzAfterRejection.second.second;

	vec_idRJAlg.push_back(idRJAlg);

	//incremento contatori per il conteggio di quanti jet sono assegnati 
 	if(idRJAlg == 0) NJetRJAlgFromPV++; //il vertice primario occupa il primo posto nel vettore di vertici
       	
 	if(idRJAlg > 0) NJetRJAlgFromPU++; //i vertici di PU occupano i posti nel vettore di vertici successivi al primo
 
	//Histogram of minimum dz jet - vertices after rejection (5 sigma)
	H_DzAfterRejection_->Fill(dzRJAlg);


	//-----------------------------------
	//In funzione del taglio di reiezione
	//-----------------------------------
	
	for(int i = 0; i < nAccettanza ; i++){
	  //Rejection algorithm
	  RejectionAlg tkRejectionSFunc(vec_RecoTk);
	  SimpleTrackCollection recoTkSFunc(tkRejectionSFunc.eval(i));
	  //Weight averager function
	  pair< int, pair< double, double > > zJetRJAlg_SFunc( wAverager(recoTkSFunc));
	  //dZ associator: Z( weight average after rejection ) - vertices 
	  pair<pair<int, double>, pair<int, double> > dzRejectionNS( dzAssociator( vec_Vertices, zJetRJAlg_SFunc.second.first ) );
	  //Histograms of dZ:  Z jet( weight average after rejection of i sigma ) - vertices 
	  H_dZJet_RJAlg_SFunc_[i]->Fill(dzRejectionNS.first.second);
	  H_Ntk_vs_NS_->Fill( i, zJetRJAlg_SFunc.first );
	}
	
	Ntk += NTks_S3;
	Njet++;
      }//end cut NTks_S3 != 0
    }// end cut Et(jet>25)
    
    ++caloJets_it;
  }//end loop calojet

  cout<<"# jet= "<<Njet<<endl;
  cout<<"# tracks= "<<Ntk<<endl;

  
  //Histograms # of jet with (without) tracks in the sample  
  H_NJet_withTracks_->Fill(NJet_withTracks); //# jet nell'evento con tracce assegnate
  H_NJet_withNoTracks_->Fill(NJet_withNoTracks); //# jet nell'evento che non ha tracce assegnate
  
  //---------------------------------
  // # jet assigned to a vertex (ok)
  //---------------------------------
  
  //Primary vertex Z
  double PV_Z = 0.;
  PV_Z =  vec_Vertices.begin()->z(); 
  //vec_Vertices size  
  int Vertices_size = 0;
  Vertices_size = vec_Vertices.size();
  //vettore del # di jet assegnati al corrispondente vertice
  std::vector<int> vec_Nj_Vtx_WAvg;
  
  // Conto quanti jet ho per ogni vertice identificato dal numero Vmin
  // riempio un vettore con questi valori
  for(int N_Vtx_WAvg_id = 0 ;  N_Vtx_WAvg_id < Vertices_size;  N_Vtx_WAvg_id++){
    int Nj_Vtx_temp = 0;
    std::vector<int>::const_iterator Vmin_it = vec_Vmin.begin();
    for(; Vmin_it != vec_Vmin.end() ; Vmin_it++){
      if((*Vmin_it) ==  N_Vtx_WAvg_id) Nj_Vtx_temp++;
    }
    vec_Nj_Vtx_WAvg.push_back(Nj_Vtx_temp);
  }
  //  cout<<"Vector NjVtx size=  "<<vec_Nj_Vtx_WAvg.size()<<endl;

  //___________________________________
  //
  //  REJECTION ALGORITHM
  //___________________________________
  

  //---------------------------------
  // # jet assigned to a vertex 
  // After Rejection Algorithm
  //---------------------------------
  
  //vettore del # di jet assegnati al corrispondente vertice
  std::vector<int> vecNjRJAlgVtx;
  
  // Conto quanti jet ho per ogni vertice identificato dal numero Vmin
  // riempio un vettore con questi valori
  for(int idVtx = 0 ;  idVtx < Vertices_size;  idVtx++){
    int tempNjRJAlg = 0;
    std::vector<int>::const_iterator vec_idRJAlgIt = vec_idRJAlg.begin();
    for(; vec_idRJAlgIt != vec_idRJAlg.end() ; vec_idRJAlgIt++){
      if((*vec_idRJAlgIt) ==  idVtx) tempNjRJAlg++;
    }
    vecNjRJAlgVtx.push_back(tempNjRJAlg);
  }
  cout<<"Size= "<<vecNjRJAlgVtx.size()<<endl;
  
  //-----------------------------
  //Dz vertices - primary vertex (ok)
  //-----------------------------

  //In un loop sui vertici considero le differenze tra la z di questi e quella del vertice primario
  //associando il # di jet assegnati a quel vertice, informazione ricavata sopra e contenuta nel vettore su cui loop
  //----------------------------------------------------------------------------------------------------------------
//   //loop on vector # of jet assigned to a vertex and vertices 
//   std::vector<int>::const_iterator vec_Nj_Vtx_WAvg_it = vec_Nj_Vtx_WAvg.begin();
//   vector<math::XYZTLorentzVector>::const_iterator Vertices_iter_3 = vec_Vertices.begin();
//   for(; vec_Nj_Vtx_WAvg_it != vec_Nj_Vtx_WAvg.end(); vec_Nj_Vtx_WAvg_it++, Vertices_iter_3++){
//     if( Vertices_iter_3 != vec_Vertices.begin()){ // if not primary vertex
      
//       double Vtx_Z_WAvg_temp = 0.;
//       double Diff_Vtx_WAvg_temp = 0.;
//       int Nj_Vtx_WAvg_temp = 0;    
      
//       Vtx_Z_WAvg_temp = Vertices_iter_3->z();
//       Diff_Vtx_WAvg_temp =  Vtx_Z_WAvg_temp - PV_Z;
//       Nj_Vtx_WAvg_temp = (*vec_Nj_Vtx_WAvg_it);
      
//       //TProfile Dz vertex - primary vertex VS # jet assigned to the vertex     
//       H_Dz_WAvg_Vtx_vs_Nj_->Fill(Diff_Vtx_WAvg_temp,(double)Nj_Vtx_WAvg_temp);
//       H_Z_WAvg_Vtx_vs_Nj_->Fill(Vtx_Z_WAvg_temp,(double)Nj_Vtx_WAvg_temp);

      
//     } //end if not primary vertex
//   } //end loop on vector # of jet assigned to a vertex and vertices  

  //loop on vector # of jet assigned to a vertex and vertices 
  std::vector<int>::const_iterator vec_Nj_Vtx_WAvg_it = vec_Nj_Vtx_WAvg.begin();
  std::vector<int>::const_iterator vecNjRJAlgVtxIt = vecNjRJAlgVtx.begin();
  vector<math::XYZTLorentzVector>::const_iterator Vertices_iter_3 = vec_Vertices.begin();
  for(; Vertices_iter_3 != vec_Vertices.end(); vec_Nj_Vtx_WAvg_it++, Vertices_iter_3++, vecNjRJAlgVtxIt++ ){
    if( Vertices_iter_3 != vec_Vertices.begin()){ // if not primary vertex
      
      double Vtx_Z_WAvg_temp = 0.;
      int Nj_Vtx_WAvg_temp = 0;    
      
      double tempDzVtx =  Vertices_iter_3->z() - PV_Z;
      
      Vtx_Z_WAvg_temp = Vertices_iter_3->z();
      Nj_Vtx_WAvg_temp = (*vec_Nj_Vtx_WAvg_it);
      
      //TProfile Dz vertex - primary vertex VS # jet assigned to the vertex     
      H_Dz_WAvg_Vtx_vs_Nj_->Fill(tempDzVtx,(double)Nj_Vtx_WAvg_temp);
      H_Z_WAvg_Vtx_vs_Nj_->Fill(Vtx_Z_WAvg_temp,(double)Nj_Vtx_WAvg_temp);
            
      int tempNjRJAlgVtx = (*vecNjRJAlgVtxIt);
      
      //TProfile Dz vertex - primary vertex VS # jet assigned to the vertex     
      H_ProfileDzVtxVsNjRJAlg_->Fill( tempDzVtx,(double)tempNjRJAlgVtx);
      
    } //end if not primary vertex   
  } //end loop on vector # of jet assigned to a vertex and vertices  
  


//   std::vector<int>::const_iterator vecNjRJAlgVtxIt = vecNjRJAlgVtx.begin();
//   vector<math::XYZTLorentzVector>::const_iterator VerticesIter4 = vec_Vertices.begin();
//   for(; vecNjRJAlgVtxIt != vecNjRJAlgVtx.end();vecNjRJAlgVtxIt++, VerticesIter4++){
//     if( VerticesIter4 != vec_Vertices.begin()){ // if not primary vertex
      
//       double tempDzVtx =  VerticesIter4->z() - PV_Z;
//       int tempNjRJAlgVtx = (*vecNjRJAlgVtxIt);
      
//       //TProfile Dz vertex - primary vertex VS # jet assigned to the vertex     
//       H_ProfileDzVtxVsNjRJAlg_->Fill( tempDzVtx,(double)tempNjRJAlgVtx);
      
//     }
//   }


  //-------------------------
  // Minimum Dz vertices (ok)
  //-------------------------
  vector<double> zVertices;
  //Loop on vertices
  vector<math::XYZTLorentzVector>::const_iterator Vertices_iter_1 = vec_Vertices.begin();
  for(; Vertices_iter_1 != vec_Vertices.end() ; Vertices_iter_1++){
    //riempio un vettore con le z del vertice    
    zVertices.push_back(Vertices_iter_1->z());
  }//end loop on vertices

  //lo ordino con ordine crescente
  sort(zVertices.begin(),zVertices.end());
  
  //inizializzo alla prima differenza tra z
  double Dz_Vtx =  *(zVertices.begin()+1) -  *(zVertices.begin()); 
  
  //  cout<<"dz_vtx= "<< Dz_Vtx<<endl;
  
  //loop sulle altre differenze fino all'elemento prima della fine
  //parto dal secondo
  std::vector<double>::const_iterator zVerticesIt =  zVertices.begin()+1;
  for(; zVerticesIt != zVertices.end() ; zVerticesIt++){
    //se il seguente non e' la fine, quindi se l'attuale non e' l'ultimo calcolo la differenza
    if(zVerticesIt+1 != zVertices.end()){
      double Dz_Vtx_temp = (*(zVerticesIt+1)  - *(zVerticesIt));
      //se la differenza tra attuale e seguente e' minore di quella tra attuale e precedente riempio con questa
      if(fabs ( Dz_Vtx_temp ) < fabs( Dz_Vtx )){
	//Histograms for Dz vertices 
  	H_Dz_Vtx_->Fill( Dz_Vtx_temp );
      }else { //se invece e' maggiore riempio con la precedente
	H_Dz_Vtx_->Fill( Dz_Vtx  );
      } //end else
      //aggiorno alla distanza tra attuale e successivo prima di passare al successivo 
      Dz_Vtx =  Dz_Vtx_temp ;
    }else if(zVerticesIt+1 == zVertices.end()){ // se il successivo e' l'ultimo il precedente era la distanza minore e riempio con quella
      H_Dz_Vtx_->Fill( Dz_Vtx  );
    } //end else if
  }//end loop on zVertices
 
  //---------------------------------------
  // Phases Space: Dz PU vertices Z (ok)
  //---------------------------------------

  //Considero tutte le differenze in z, una volta sola imponendo la condizione che l'iteratore del loop interno sia 
  //sempre maggiore di quello del loop esterno
  //----------------------------------------------------------------------------------------------------------------
  //first loop on vertices
  vector<math::XYZTLorentzVector>::const_iterator  i_vec_PUVertex_it = vec_Vertices.begin(); 
  for(; i_vec_PUVertex_it != vec_Vertices.end(); i_vec_PUVertex_it++){
    if(i_vec_PUVertex_it != vec_Vertices.begin()){

      //second loop on vertices
      std::vector<int>::const_iterator vec_Nj_Vtx_WAvg_it = vec_Nj_Vtx_WAvg.begin();
      vector<math::XYZTLorentzVector>::const_iterator j_vec_PUVertex_it = vec_Vertices.begin();
      for(; j_vec_PUVertex_it != vec_Vertices.end(); j_vec_PUVertex_it++,vec_Nj_Vtx_WAvg_it++){
	//if vertex 1 != and < vertex 2 and # jet don't  primary vertex 
	if ( i_vec_PUVertex_it < j_vec_PUVertex_it ){ //&& vec_Nj_Vtx_WAvg_it != vec_Nj_Vtx_WAvg.begin()){
	  double Zi = 0.;
	  double Zj = 0.;
	
	  Zi= i_vec_PUVertex_it->z();
	  Zj= j_vec_PUVertex_it->z();
	  
	  //Histogram (phases space) for Dz jet - PU vertices	
	  H_Profile_WAvg_DzVtx_->Fill(Zj - Zi, (*vec_Nj_Vtx_WAvg_it));
	  
	}//end if
      }//end secondo loop sui vertici di pileup


    }//end if 
  }//end primo loop sui vertici di pileup
  
  //Histograms # jet assigned to PV end PU vertices in the event
  H_NJet_fromPV_->Fill(NJet_fromPV); //# jet dell'evento assegnati al vertice primario
  H_NJet_fromPU_->Fill(NJet_fromPU); //# jet dell'evento assegnati ai vertici di PU
  H_NJet_NotAssigned_->Fill(NJet_NotAssigned); //# jet dell'evento non assegnati 
  
  
  H_NJetRJAlgFromPV_->Fill(NJetRJAlgFromPV); //# jet dell'evento assegnati al vertice primario
  H_NJetRJAlgFromPU_->Fill(NJetRJAlgFromPU); //# jet dell'evento assegnati ai vertici di PU


  
  //=====ROBERTO=================================
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void VertexAssoc::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void VertexAssoc::endJob() {
  
  
  //=====ROBERTO=============

  //Histograms for primary vertices and PU vertices
  H_NPUVtx_->Write();

  //Histograms for Z PV vertices
  H_primaryVtxZ_->Write(); 
  //Histograms for Z PU vertices  
  H_PUVtxZ_->Write(); 
  //Histograms for Z simTracks
  H_Z_simTks_->Write();
  //Histograms for Z reco tracks 
  H_TksZ_->Write(); 
  H_TksZ_Error_->Write();

  //Histogram Z jet  
  H_TksZ_Weighted_Avg_->Write();
  H_TksZ_WAvg_Error_->Write();

  //Histograms selected tracks  
  H_SelTk_->Write();
  
  //Histograms for 2D, 3D IP Significance
  H_S2DTks_->Write();
  H_S3DTks_->Write();
  H_Tks_Z_S3D_->Write();

  //Histograms for Njet  
  H_NJet_withTracks_->Write();
  H_NJet_withNoTracks_->Write();
  H_NJet_fromPV_->Write();
  H_NJet_fromPU_->Write();
  H_NJet_NotAssigned_->Write();
  
  //Histograms for Dz jet - vertices 
  H_Vmin_->Write();
  H_Dz_WAvg_Jet_Vtx_tot_->Write();
  H_Dz_WAvg_Jet_Vtx_2Best_tot_ ->Write();
  H_Dz_WAvg_Jet_Vtx_vs_2Best_->Write();

  //Histograms for Dz vertices
  H_Dz_Vtx_->Write();
  
  //Histogram for Dz jet - vertices 
  H_Dz_WAvg_Vtx_vs_Nj_->Write();
  
  H_Z_WAvg_Vtx_vs_Nj_ ->Write();

  //Histogram (phases space) for Dz jet - PU vertices
  H_Profile_WAvg_DzVtx_->Write();
    
  //Histograms Dz jet - vertices for exclusive # of tracks in the jet 
  for (int i = 0; i < 13 ; i++){
    H_Z_WAvg_Jet_Ntracks_[i]->Write();
  }
  
  //Histograms for eta, phi reco tracks
  H_TksEta_->Write();
  H_TksPhi_->Write();
  //Histograms for Pt, Chi2 reco tracks
  H_TksPt_->Write();
  H_TksChi2_->Write();
 
  //Histograms for Et, eta, phi of calojets
  H_Et_Jets_->Write();
  H_eta_Jets_->Write();
  H_phi_Jets_->Write();

  //Histograms Dphi, Deta, Delta2 jet - reco selected tracks 
  H_Dphi_Jet_Tk_->Write();
  H_Deta_Jet_Tk_->Write();
  H_D2_Jet_Tk_->Write();

  //Histograms for eta, phi sim tracks
  H_Phi_simTk_->Write();
  H_Eta_simTk_->Write();

  //Histograms for Delta2 minimum simTrack - reco selected track
  H_D2_Tk_sim_reco_->Write();     
  H_D2_Tk_sim_reco_2Best_->Write();
  H_D2_Tk_vs_2Best_->Write(); 
  //Histograms of Dz simTrack - reco selected track   
  H_Dz_Tk_sim_reco_->Write();

  H_Z_sim_jet_->Write();
  H_Dz_jet_sim_reco_->Write();
  H_Dz_jet_sim_reco_avg_->Write();

  //Histogram of Dz selected sim tracks
  H_Dz_sim_selected_tk_->Write();
  //TProfile max rate of same Z in selected sim tracks in function of the # of tracks
  H_Profile_Ntk_vs_sameZ_->Write();
  //Histograms of rate of same selected sim tracks Z and 
  //rate of same selected sim tracks Z vs # of tracks in the jet
  H_Frequenze_->Write();
  H_Ntk_vs_SameZ_->Write();
  //Histogram of selected simulated tracks' Z
  H_Z_selected_sim_tracks_->Write();

  //Histograms of Dz jet - vertex: |Z(true) - Z(vtx)| 
  H_Dz_sim_Jet_Vtx_->Write();
  H_Dz_sim_Jet_Vtx_2Best_->Write();
  H_Dz_sim_Jet_Vtx_vs_2Best_->Write();

  //Histograms Dz jet - vertices for exclusive # of tracks in the jet    
  for (int i = 0; i < 13 ; i++){
    H_Dz_sim_Jet_Vtx_NTk_[i]->Write();
  }

  //Histogram of |Z(vera) - Z(reco)| jet
  H_Dz_jet_true_reco_->Write();
 
  //Histograms of Dz vera jet - Z sim vertex associated
  H_Dz_sim_jet_simVtx_ ->Write();        
  H_Dz_sim_jet_simVtx_2Best_ ->Write();  
  H_Dz_sim_jet_simVtx_vs_2Best_->Write();

  //Histograms of Z jet and error after the rejection algorithm (5 sigma)
  H_ZJet_RJAlg_->Write();
  H_ZErrorJet_RJAlg_->Write();
  H_DzAfterRejection_->Write();

  for(int i = 0; i < 20; i++){
    H_dZJet_RJAlg_SFunc_[i]->Write();  
    //    H_ZErrorJet_RJAlg_SFunc_[i]->Write();  
  }
  //Histogram of Ntk in function # sigma rejection
  H_Ntk_vs_NS_->Write();

  
  H_NJetRJAlgFromPV_->Write(); 
  H_NJetRJAlgFromPU_->Write(); 

  H_ProfileDzVtxVsNjRJAlg_->Write();

  //=====ROBERTO=============
  
}


// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(VertexAssoc);
