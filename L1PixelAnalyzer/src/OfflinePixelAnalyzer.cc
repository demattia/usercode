//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/OfflinePixelAnalyzer.h"

// Classes to be accessed
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimplePixelJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/Summary.h"

#include "AnalysisExamples/PixelJet/interface/PixelJet.h"

// For file output
#include <fstream>
#include <sstream>

#include <cmath>

#include <memory>

//
// constants, enums and typedefs
//

//
// static data member definitions
//

L1Trig OfflinePixelAnalyzer::L1Trigger;

//
// constructors and destructor
//
OfflinePixelAnalyzer::OfflinePixelAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  cenJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "CenJets" ) ),
  forJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "ForJets" ) ),
  tauJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "TauJets" ) ),
  l1MEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "L1MEt" ) ),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  simplePixelJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimplePixelJets" ) ),
  globalMuonLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons" ) ),
  simpleElectronLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons" ) ),
  simpleTauLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus" ) ),
  summaryLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "Summary" ) ),
  simVtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx" ) ),
  numTkCut_( iConfig.getUntrackedParameter<unsigned int>( "TracksMinimumNum_in_PixelJet" ) ),
  minDz_( iConfig.getUntrackedParameter<double>( "MinDz" ) ),
  maxDz_( iConfig.getUntrackedParameter<double>( "MaxDz" ) ),
  doTrigger_( iConfig.getUntrackedParameter<bool>( "doTrigger" ) ),
  extendedInfo_( iConfig.getUntrackedParameter<bool>( "extendedInfo" ) ),
  OutputEffFileName( iConfig.getUntrackedParameter<string>( "OutputEffFileName" ) )
{
  //now do what ever initialization is needed

  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,"RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  OutputFile->cd();

  // White background for the canvases
  gROOT->SetStyle("Plain");

  PixelJet_dz_ = new TH1F( "PixelJet_dz", "PixelJet minimum dz", 100, 0, 10 );
  PixelJet_Num_ = new TH1F( "PixelJet_Num", "Number of pixeljets", 30, 0, 30 );
  PixelJet_Track_Num_ = new TH1F( "PixelJet_Track_Num", "Number of pixeltracks in pixeljets", 15, 0, 15 );

  // Generate histograms for the verteces
  // ------------------------------------
  dz_ = 0.1;
  bins_ = 30;

  dzmax_ = dz_ + dz_*(bins_);

  // Multiple histograms, including mean and stacked histograms
  Multi_Vertex_Dz_ = new MultiTH1F( "Vertex_Dz", "Minimum distance between verteces", 100, 0., 10., bins_, dz_, dzmax_, OutputFile );
  Multi_Prim_Second_Vertex_Dz_ = new MultiTH1F( "Prim_Second_Vertex_Dz", "Minimum distance between primary and secondary verteces", 100, 0., 10., bins_, dz_, dzmax_, OutputFile );
  Multi_Vertex_Num_ = new MultiTH1F( "Vertex_Num", "Number of verteces", 30, 0., 30., bins_, dz_, dzmax_, OutputFile );

  Multi_PrimVNum_ = new MultiTH1F( "PrimVNum", "Number of pixeljets in primary vertex", 30, 0, 30, bins_, dz_, dzmax_, OutputFile );
  Multi_PrimVPt_ = new MultiTH1F( "PrimVPt" , "Pt of primary vertex", 100, 0, 1000, bins_, dz_, dzmax_, OutputFile );
  Multi_PrimVEta_ = new MultiTH1F( "PrimVEta", "Eta of primary vertex", 100, -3, 3, bins_, dz_, dzmax_, OutputFile );
  Multi_PrimVPhi_ = new MultiTH1F( "PrimVPhi", "Phi of primary vertex", 100, -3.15, 3.15, bins_, dz_, dzmax_, OutputFile );

  Multi_SecVNum_ = new MultiTH1F( "SecVNum", "Number of pixeljets in secondary vertex" , 30, 0, 30, bins_, dz_, dzmax_, OutputFile );
  Multi_SecVPt_ = new MultiTH1F( "SecVPt", "Pt of secondary vertex", 100, 0, 1000, bins_, dz_, dzmax_, OutputFile );
  Multi_SecVEta_ = new MultiTH1F( "SecVEta", "Eta of secondary vertex", 100, -3, 3, bins_, dz_, dzmax_, OutputFile );
  Multi_SecVPhi_ = new MultiTH1F( "SecVPhi", "Phi of secondary vertex", 100, -3.15, 3.15, bins_, dz_, dzmax_, OutputFile );

  Multi_AllSecVNum_ = new MultiTH1F( "AllSecVNum", "Number of pixeljets in all secondary verteces", 30, 0, 30, bins_, dz_, dzmax_, OutputFile );
  Multi_AllSecVPt_ = new MultiTH1F( "AllSecVPt_Name", "Pt of all secondary verteces", 100, 0, 1000, bins_, dz_, dzmax_, OutputFile );
  Multi_AllSecVEta_ = new MultiTH1F( "AllSecVEta_Name", "Eta of all secondary verteces", 100, -3., 3, bins_, dz_, dzmax_, OutputFile );
  Multi_AllSecVPhi_ = new MultiTH1F( "AllSecVPhi_Name", "Phi of all secondary verteces", 100, -3.15, 3.15, bins_, dz_, dzmax_, OutputFile );

  MultiVertexDeltaZ_ = new MultiTH1F( "vertexDeltaZ", "DeltaZ between SimVertex and PJVertex", 100, -3, 3, bins_, dz_, dzmax_, OutputFile );
  MultiVertexDeltaZres_ = new MultiTH1F( "vertexDeltaZres", "Resolution of Z of PJVertex", 100, -3, 3, bins_, dz_, dzmax_, OutputFile );

  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 100, 0, 3.15 );

  eventcounter_ = 0;

  numgoodpjeff_ = new int[20];
  numgoodpjeff_3_ = new int[20];
  numgoodpjeff_4_ = new int[20];
  numgoodpjeff_5_ = new int[20];
  numgoodpjeff_6_ = new int[20];

  EffNumGoodPj_ = new TH1F ( "EffNumGoodPj", "Efficiency of cut on number of good pixeljets", 20, 0, 20 );
  EffNumGoodPj_3_ = new TH1F ( "EffNumGoodPj_3", "Efficiency of cut on number of good pixeljets && 3 pj from pv", 20, 0, 20 );
  EffNumGoodPj_4_ = new TH1F ( "EffNumGoodPj_4", "Efficiency of cut on number of good pixeljets && 4 pj from pv", 20, 0, 20 );
  EffNumGoodPj_5_ = new TH1F ( "EffNumGoodPj_5", "Efficiency of cut on number of good pixeljets && 5 pj from pv", 20, 0, 20 );
  EffNumGoodPj_6_ = new TH1F ( "EffNumGoodPj_6", "Efficiency of cut on number of good pixeljets && 6 pj from pv", 20, 0, 20 );

  for (int i=0; i<20; ++i) {
    numgoodpjeff_[i] = 0;
    numgoodpjeff_3_[i] = 0;
    numgoodpjeff_4_[i] = 0;
    numgoodpjeff_5_[i] = 0;
    numgoodpjeff_6_[i] = 0;
  }
}


OfflinePixelAnalyzer::~OfflinePixelAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // Draw the histograms
  //  HiVar->Plot();

  OutputFile->Write();

  delete[] numgoodpjeff_;
  delete[] numgoodpjeff_3_;
  delete[] numgoodpjeff_4_;
  delete[] numgoodpjeff_5_;
  delete[] numgoodpjeff_6_;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
OfflinePixelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace anaobj;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // ------- //
  // Offline //
  // ------- //

//   // Take genParticleCandidates
//   edm::Handle < MCParticleCollection > MCpartons;
//   iEvent.getByLabel( MCParticleLabel_, MCpartons );
//   MCParticleCollection::const_iterator mc_it = MCpartons->begin();
//   int mcCounter = 0;
//   for( ; mc_it != MCpartons->end(); ++mc_it, ++mcCounter ) {
//     cout << "parton["<<mcCounter<<"] vertex z = " << mc_it->vz() << endl;
//   }

  // Take genParticleCandidates
    Handle< CandidateCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  int mcCounter = 0;
  for( CandidateCollection::const_iterator mc_it = MCpartons->begin(); mc_it != MCpartons->end(); ++mc_it, ++mcCounter ) {
    const Candidate* temp_mother = mc_it->mother();
    if ( temp_mother != 0 ) {

      int pid = abs( mc_it->pdgId() );
      int Mpid = abs( mc_it->mother()->pdgId() );

      // If is a top(6), Z(23), W(24) or Higgs(25) or a daughter of those
      if ( pid  == 6 || pid  == 23 || pid  == 24 || pid  == 25 ||
           Mpid == 6 || Mpid == 23 || Mpid == 24 || Mpid == 25 ) {

//        vec_MC_ptr->push_back( MCParticle( mc_it->pt(), mc_it->eta(), mc_it->phi(), mc_it->mass(), mc_it->pdgId(), mc_it->mother()->pdgId() ) );



        cout << "parton id = " << pid << " and vertex z = " << mc_it->vz() << endl;




      }
      // Store also the case in which it is a b with mother != b, or a c with mother != b and c
      else if ( ( pid == 5 && Mpid != 5 ) || ( pid == 4 && Mpid != 5 && Mpid != 4) ) {
//        vec_MC_ptr->push_back( MCParticle( mc_it->pt(), mc_it->eta(), mc_it->phi(), mc_it->mass(), mc_it->pdgId(), mc_it->mother()->pdgId() ) );
      }
    }
  } // end loop on MC particles

  // Take the SimVertex collection
  std::vector<SimVertex> theSimVerteces;
  Handle<SimVertexContainer> SimVtx;
  iEvent.getByLabel( simVtxLabel_, SimVtx );

//  theSimVerteces.insert(theSimVerteces.end(),SimVtx->begin(),SimVtx->end());
//  cout << "SimVertex size = " << theSimVerteces.size() << endl;

  cout << "SimVertex size = " << SimVtx->size() << endl;


  // Store the z of the simVertex
  float simVertexZ = SimVtx->begin()->position().z();


  //  for( vector<SimVertex>::const_iterator vec_it = theSimVerteces.begin(); vec_it != theSimVerteces.end(); ++vec_it ) {
//   for( vector<SimVertex>::const_iterator vec_it = SimVtx->begin(); vec_it != SimVtx->end(); ++vec_it ) {
//     cout << "Sim Vertex z = " << vec_it->position().z() << endl;
//   }

  // Pixel jets
  // ----------

  Handle<PixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );

  const PixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // PixelJets with at least numTkCut_ tracks
  PixelJetCollection goodPixelJets;

  // Count the number of pixeljets with at least numTkCut_ tracks
  PixelJetCollection::const_iterator pj_it = pixeljets.begin();
  int goodPjNum = 0;
  for ( ; pj_it != pixeljets.end(); ++pj_it ) {
    if ( pj_it->tkNum() >= int(numTkCut_) ) {
      ++goodPjNum;
      goodPixelJets.push_back( *pj_it );
    }
  }

  // Pixel trigger requiring at least 2 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)

  for( int numdz = 0; numdz < bins_; ++numdz ) {

    double dz = dz_ + dz_*(numdz);

    L1PixelTrig<PixelJet> PixelTrig_2(2, dz, numTkCut_);

    // Performs the cut on the number of track inside
    PixelTrig_2.Fill( pixeljets );

    // The vector of verteces is already sorted in ascending Pt
    auto_ptr<vector<Vertex<PixelJet> > > vec_vertex_aptr( PixelTrig_2.VevVertex() );

    // Number of verteces
    int numVertex = 0;
    if ( vec_vertex_aptr.get() != 0 ) {
      numVertex = vec_vertex_aptr->size();
    }
    Multi_Vertex_Num_->Fill( numVertex, numdz );

    if ( numVertex > 0 ) {
      //Primary vertex characteristics
      Vertex<PixelJet> prim_vertex( vec_vertex_aptr->back() );

      Multi_PrimVNum_->Fill( prim_vertex.number(), numdz );
      Multi_PrimVPt_->Fill( prim_vertex.pt(), numdz );
      Multi_PrimVEta_->Fill( prim_vertex.eta(), numdz );
      Multi_PrimVPhi_->Fill( prim_vertex.phi(), numdz );

      MultiVertexDeltaZ_->Fill( simVertexZ - prim_vertex.z(), numdz );
      MultiVertexDeltaZres_->Fill( (simVertexZ - prim_vertex.z())/simVertexZ, numdz );

      // Loop on the remaining verteces and fill their characteristics
      if ( numVertex > 1 ) {
        // End returns just just after the last element, to take the one before the last go two times back
        vector<Vertex<PixelJet> >::const_iterator svec_it = (vec_vertex_aptr->end()-2);
        Multi_SecVNum_->Fill( svec_it->number(), numdz );
        Multi_SecVPt_->Fill( svec_it->pt(), numdz );
        Multi_SecVEta_->Fill( svec_it->eta(), numdz );
        Multi_SecVPhi_->Fill( svec_it->phi(), numdz );

        vector<Vertex<PixelJet> >::const_iterator all_svec_it = (vec_vertex_aptr->begin());
        for ( ; all_svec_it != vec_vertex_aptr->end()-1; ++all_svec_it ) {
          Multi_AllSecVNum_->Fill( all_svec_it->number(), numdz );
          Multi_AllSecVPt_->Fill( all_svec_it->pt(), numdz );
          Multi_AllSecVEta_->Fill( all_svec_it->eta(), numdz );
          Multi_AllSecVPhi_->Fill( all_svec_it->phi(), numdz );
        }
      }

      // Now sort the vertex in z and evaluate closest distance
      sort( vec_vertex_aptr->begin(), vec_vertex_aptr->end(), Sort_Greater_Z<Vertex<PixelJet> >() );

      if ( numVertex > 1 ) {
        vector<Vertex<PixelJet> >::const_iterator svec_z_it = (vec_vertex_aptr->begin());
        for ( ; svec_z_it != vec_vertex_aptr->end()-1; ++svec_z_it ) {
          Multi_Vertex_Dz_->Fill( fabs(svec_z_it->z() - (svec_z_it-1)->z()), numdz );
        }
      }

      // Evaluate the distance between primary and secondary vertex
      if ( numVertex > 1 ) {
        // Take the last and the one before
        vector<Vertex<PixelJet> >::const_iterator primV_secondV_z_it = (vec_vertex_aptr->end()-1);
        Multi_Prim_Second_Vertex_Dz_->Fill( fabs(primV_secondV_z_it->z() - (primV_secondV_z_it-1)->z()), numdz );      
      }
    }
  } // end dz loop

  // Fill histograms for all goodPixelJets in the event (those which have at least numTkCut_ pixelTracks)
  int numPixelJet = goodPixelJets.size();
  PixelJet_Num_->Fill( numPixelJet );
  if ( numPixelJet > 1 ) {
    vector<PixelJet>::const_iterator pj_z_it (goodPixelJets.begin()+1);
    PixelJet_Track_Num_->Fill(goodPixelJets.begin()->tkNum());
    for ( ; pj_z_it != goodPixelJets.end(); ++pj_z_it ) {
      PixelJet_Track_Num_->Fill(pj_z_it->tkNum());
      PixelJet_dz_->Fill( fabs(pj_z_it->z() - (pj_z_it-1)->z()) );
    }
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void OfflinePixelAnalyzer::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void OfflinePixelAnalyzer::endJob() {

  Multi_Vertex_Dz_->Write();
  Multi_Vertex_Num_->Write();
  Multi_Prim_Second_Vertex_Dz_->Write();

  Multi_PrimVPt_->Write();
  Multi_PrimVNum_->Write();
  Multi_PrimVEta_->Write();
  Multi_PrimVPhi_->Write();

  Multi_SecVPt_->Write();
  Multi_SecVNum_->Write();
  Multi_SecVEta_->Write();
  Multi_SecVPhi_->Write();

  Multi_AllSecVPt_->Write();
  Multi_AllSecVNum_->Write();
  Multi_AllSecVEta_->Write();
  Multi_AllSecVPhi_->Write();

  MultiVertexDeltaZ_->Write();
  MultiVertexDeltaZres_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflinePixelAnalyzer);
