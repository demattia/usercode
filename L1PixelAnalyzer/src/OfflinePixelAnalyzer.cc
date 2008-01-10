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

  HiVar = new HiVariables( ( conf_.getUntrackedParameter<string>("HiVarName") ).c_str() );

  uncorr_JetPt_IC5_ = new TH1F( "uncorr_JetPt_IC5", "uncorrected JetPt IC5", 100, 0, 200 );
  corr_JetPt_IC5_ = new TH1F( "corr_JetPt_IC5", "corrected JetPt IC5", 100, 0, 200 );
  JetNumber_IC5_ = new TH1F( "JetNumber_IC5", "Number of IC5 jets", 100, 0, 100 );

  MEt_CorrIC5_Pt_ = new TH1F( "MEt_CorrIC5_Pt", "MEt corrected with IC5 jets", 100, 0, 150 );
  MEt_CorrIC5_Phi_ = new TH1F( "MEt_CorrIC5_Phi", "MEt Phi corrected with IC5 jets", 100, -3.15, 3.15 );
  MEt_CorrIC5_SumEt_ = new TH1F( "MEt_CorrIC5_SumEt", "SumEt corrected with IC5 jets", 100, 0, 2500 );
  MEt_CorrIC5_mEtSig_ = new TH1F( "MEt_CorrIC5_mEtSig", "MEt significance corrected with IC5 jets", 100, 0, 10 );

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


  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 100, 0, 3.15 );

  eventcounter_ = 0;
  //  PI_ = 3.141593;

  // PixelTrigger alone efficiency
//   pixelTrig_3_ = 0;
//   pixelTrig_4_ = 0;
//   pixelTrig_5_ = 0;
//   pixelTrig_6_ = 0;

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
  using namespace std;
  using namespace anaobj;

  // L1 Calo
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  // should get rid of this try/catch?
  try {
    iEvent.getByLabel(cenJetLabel_, l1eCenJets);
    iEvent.getByLabel(forJetLabel_, l1eForJets);
    iEvent.getByLabel(tauJetLabel_, l1eTauJets);
    iEvent.getByLabel(l1MEtLabel_, l1eEtMiss);
  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  // Global muons
  edm::Handle < GlobalMuonCollection > globalMuons;
  // SimpleElectrons
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  // SimpleTaus
  edm::Handle < SimpleTauCollection > simpleTaus;
  // Summary
  edm::Handle < Summary > summary;
  try {
    iEvent.getByLabel( globalMuonLabel_, globalMuons );
    iEvent.getByLabel( simpleElectronLabel_, simpleElectrons );
    iEvent.getByLabel( simpleTauLabel_, simpleTaus );
    iEvent.getByLabel( summaryLabel_, summary );
  }
  catch (...) {
    std::cerr << "One of the remaining collections cannot be found" << endl;
  }

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // ------- //
  // Offline //
  // ------- //

  // MEt
  // ---

  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  MEt_CorrIC5_Pt_->Fill( caloMET->et() );
  double MET_phi = caloMET->phi();
  MEt_CorrIC5_Phi_->Fill( MET_phi );
  MEt_CorrIC5_SumEt_->Fill( caloMET->sumEt() );
  MEt_CorrIC5_mEtSig_->Fill( caloMET->mEtSig() );

  //  Associator<CaloMET, CaloJet> associator( 0.5 );
  //  std::auto_ptr<std::map<const reco::CaloMET*, const reco::CaloJet*> > AssocMap( associator.Associate( *caloMET, *caloJets ) );

  // Take the first (and it should be the only one) element in the map and get the eta of the closest jet (in DeltaR)
  //   double JetPhi = (*AssocMap).begin()->second->phi();

  // Do not use DR association. The MEt has no z component. Use DeltaPhi association.

  // HiVariables
  // -----------

  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Count IC5 jets with Et>=30GeV and |eta| < 3.0
  int goodIc5Jets = 0;

  // Offline cuts
  bool offline = false;

  if ( caloJets->size() != 0 ) {

    std::vector<double> vec_DPhi;
    vec_DPhi.reserve(caloJets->size());

    vector<SimpleJet> vec_calojet; 
    vec_calojet.reserve(caloJets->size());
    // Correct offline jets on the fly
    for( OfflineJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal ) {
      vec_calojet.push_back( SimpleJet( cal->et(), cal->eta(), cal->phi() ) );
      uncorr_JetPt_IC5_->Fill( cal->uncorrEt() );
      corr_JetPt_IC5_->Fill( cal->et() );

      // Evaluate DeltaPhi between MET and calo-jets
      // Consider only high-luminosity--good-jets
      //       if ( cal->et() >= 40. && fabs( cal->eta() ) < 3.0 ) {
      //         vec_DPhi.push_back( DeltaPhi( MET_phi, cal->phi() ) );
      //       }
      if ( cal->et() >= 30. && fabs( cal->eta() ) < 3.0 ) {
        ++goodIc5Jets;
      }
    }

    // Offline cut
    // -----------
    if ( goodIc5Jets >= 5 ) offline = true;


    JetNumber_IC5_->Fill( vec_calojet.size() );

    HiVar->Fill( vec_calojet );

    // Minimum DeltaPhi
    std::sort( vec_DPhi.begin(), vec_DPhi.end() );
    DPhimin_->Fill( *(vec_DPhi.begin()) );

  }
  else {
    std::cout << "ATTENTION: Jet collection empty" << std::endl;
  }


  // Take genParticleCandidates
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  // Take the mothers
  // MCParticleCollection::const_iterator MCp = MCpartons->begin();

  // Pixel jets
  // ----------

  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );

  const SimplePixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // PixelJets with at least numTkCut_ tracks
  SimplePixelJetCollection goodPixelJets;

  // Count the number of pixeljets with at least numTkCut_ tracks
  SimplePixelJetCollection::const_iterator pj_it = pixeljets.begin();
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

    L1PixelTrig<SimplePixelJet> PixelTrig_2(2, dz, numTkCut_);

    // Performs the cut on the number of track inside
    PixelTrig_2.Fill( pixeljets );

    // Store the response for the different requirements on the number of pixelJets from the PV
    //   bool pixelTrigResponse[10];

    //   pixelTrigResponse[0] = PixelTrig_2.Response();


    // The number of pixelJets from the PV is stored in the following Multi_PrimVNum_


    // The vector of verteces is already sorted in ascending Pt
    auto_ptr<vector<Vertex<SimplePixelJet> > > vec_vertex_aptr( PixelTrig_2.VevVertex() );

    // Number of verteces
    int numVertex = 0;
    if ( vec_vertex_aptr.get() != 0 ) {
      numVertex = vec_vertex_aptr->size();
    }
    Multi_Vertex_Num_->Fill( numVertex, numdz );

    if ( numVertex > 0 ) {
      //Primary vertex characteristics
      Vertex<SimplePixelJet> prim_vertex( vec_vertex_aptr->back() );

      Multi_PrimVNum_->Fill( prim_vertex.number(), numdz );
      Multi_PrimVPt_->Fill( prim_vertex.pt(), numdz );
      Multi_PrimVEta_->Fill( prim_vertex.eta(), numdz );
      Multi_PrimVPhi_->Fill( prim_vertex.phi(), numdz );

      // Loop on the remaining verteces and fill their characteristics
      if ( numVertex > 1 ) {
        // End returns just just after the last element, to take the one before the last go two times back
        vector<Vertex<SimplePixelJet> >::const_iterator svec_it = (vec_vertex_aptr->end()-2);
        Multi_SecVNum_->Fill( svec_it->number(), numdz );
        Multi_SecVPt_->Fill( svec_it->pt(), numdz );
        Multi_SecVEta_->Fill( svec_it->eta(), numdz );
        Multi_SecVPhi_->Fill( svec_it->phi(), numdz );

        vector<Vertex<SimplePixelJet> >::const_iterator all_svec_it = (vec_vertex_aptr->begin());
        for ( ; all_svec_it != vec_vertex_aptr->end()-1; ++all_svec_it ) {
          Multi_AllSecVNum_->Fill( all_svec_it->number(), numdz );
          Multi_AllSecVPt_->Fill( all_svec_it->pt(), numdz );
          Multi_AllSecVEta_->Fill( all_svec_it->eta(), numdz );
          Multi_AllSecVPhi_->Fill( all_svec_it->phi(), numdz );
        }
      }

      // Now sort the vertex in z and evaluate closest distance
      sort( vec_vertex_aptr->begin(), vec_vertex_aptr->end(), Sort_Greater_Z<Vertex<SimplePixelJet> >() );

      if ( numVertex > 1 ) {
        vector<Vertex<SimplePixelJet> >::const_iterator svec_z_it = (vec_vertex_aptr->begin());
        for ( ; svec_z_it != vec_vertex_aptr->end()-1; ++svec_z_it ) {
          Multi_Vertex_Dz_->Fill( fabs(svec_z_it->z() - (svec_z_it-1)->z()), numdz );
        }
      }

      // Evaluate the distance between primary and secondary vertex
      if ( numVertex > 1 ) {
        // Take the last and the one before
        vector<Vertex<SimplePixelJet> >::const_iterator primV_secondV_z_it = (vec_vertex_aptr->end()-1);
        Multi_Prim_Second_Vertex_Dz_->Fill( fabs(primV_secondV_z_it->z() - (primV_secondV_z_it-1)->z()), numdz );      
      }
    }

//     // PixelTrigger alone efficiency
//     // -----------------------------

//     L1PixelTrig<SimplePixelJet> PixelTrig_3(3, dz, numTkCut_);
//     L1PixelTrig<SimplePixelJet> PixelTrig_4(4, dz, numTkCut_);
//     L1PixelTrig<SimplePixelJet> PixelTrig_5(5, dz, numTkCut_);
//     L1PixelTrig<SimplePixelJet> PixelTrig_6(6, dz, numTkCut_);

//     pixelTrigResponse[0] = PixelTrig_3.Response();
//     pixelTrigResponse[1] = PixelTrig_4.Response();
//     pixelTrigResponse[2] = PixelTrig_5.Response();
//     pixelTrigResponse[3] = PixelTrig_6.Response();

//     if (pixelTrigResponse[0]) ++pixelTrig_3_;
//     if (pixelTrigResponse[1]) ++pixelTrig_4_;
//     if (pixelTrigResponse[2]) ++pixelTrig_5_;
//     if (pixelTrigResponse[3]) ++pixelTrig_6_;

//     // Cut on the number of good pixeljets
//     for ( int numgoodpjcut=0; numgoodpjcut<20; ++numgoodpjcut ) {
//       if ( goodPjNum >= numgoodpjcut ) {
//         ++numgoodpjeff_[numgoodpjcut];
//         if (pixelTrigResponse[0]) ++numgoodpjeff_3_[numgoodpjcut];
//         if (pixelTrigResponse[1]) ++numgoodpjeff_4_[numgoodpjcut];
//         if (pixelTrigResponse[2]) ++numgoodpjeff_5_[numgoodpjcut];
//         if (pixelTrigResponse[3]) ++numgoodpjeff_6_[numgoodpjcut];
//       }
//     }

  } // end dz loop

  // Fill histograms for all goodPixelJets in the event (those which have at least numTkCut_ pixelTracks)
  int numPixelJet = goodPixelJets.size();
  PixelJet_Num_->Fill( numPixelJet );
  if ( numPixelJet > 1 ) {
    vector<SimplePixelJet>::const_iterator pj_z_it (goodPixelJets.begin()+1);
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

// void OfflinePixelAnalyzer::meanHistoSetup( TH1F * meanhisto_ptr, vector<TH1F *> vec_histo ) {
//   meanhisto_ptr->SetBinContent( numdz+1, mean[numdz]->GetMean() );
//   meanhisto_ptr->SetBinError( numdz+1, Vertex_Dz_[numdz]->GetMeanError() );
//   meanhisto_ptr->GetXaxis()->CenterLabels();
//   meanhisto_ptr->GetXaxis()->SetNdivisions(Vertex_Dz_Mean_->GetSize()-2, false);
// }

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

//  // PixelTrigger alone efficiency
//   for ( int i=0; i<20; ++i ) {
//     EffNumGoodPj_->SetBinContent(i+1, float(numgoodpjeff_[i])/float(eventcounter_));
//     EffNumGoodPj_->SetBinError(i+1, sqrt(float(numgoodpjeff_[i]))/float(eventcounter_));
//     EffNumGoodPj_3_->SetBinContent(i+1, float(numgoodpjeff_3_[i])/float(eventcounter_));
//     EffNumGoodPj_3_->SetBinError(i+1, sqrt(float(numgoodpjeff_3_[i]))/float(eventcounter_));
//     EffNumGoodPj_4_->SetBinContent(i+1, float(numgoodpjeff_4_[i])/float(eventcounter_));
//     EffNumGoodPj_4_->SetBinError(i+1, sqrt(float(numgoodpjeff_4_[i]))/float(eventcounter_));
//     EffNumGoodPj_5_->SetBinContent(i+1, float(numgoodpjeff_5_[i])/float(eventcounter_));
//     EffNumGoodPj_5_->SetBinError(i+1, sqrt(float(numgoodpjeff_5_[i]))/float(eventcounter_));
//     EffNumGoodPj_6_->SetBinContent(i+1, float(numgoodpjeff_6_[i])/float(eventcounter_));
//     EffNumGoodPj_6_->SetBinError(i+1, sqrt(float(numgoodpjeff_6_[i]))/float(eventcounter_));
//   }

}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflinePixelAnalyzer);
