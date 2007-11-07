//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/L1TrigPixelAnalyzer.h"

#include "../../PJVERTEX_CMSSW/Classes/SimpleJet/SimpleJet.h"
#include "../../PJVERTEX_CMSSW/Classes/Associator/Associator.h"
#include "../../PJVERTEX_CMSSW/Classes/DeltaPhi/DeltaPhi.h"
#include "../../PJVERTEX_CMSSW/Classes/L1PixelTrig/L1PixelTrig.h"

// For the offline jets and corrections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

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

L1Trig L1TrigPixelAnalyzer::L1Trigger;

//
// constructors and destructor
//
L1TrigPixelAnalyzer::L1TrigPixelAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
//  HiVar( ( iConfig.getUntrackedParameter<std::string> ("HiVarName") ).c_str() ),
  CaloJetAlgorithm( iConfig.getUntrackedParameter<string>( "CaloJetAlgorithm" ) ),
  JetCorrectionService( iConfig.getUntrackedParameter<string>( "JetCorrectionService" ) ),
  METCollection( iConfig.getUntrackedParameter<string>( "METCollection" ) ),
  genParticleCandidates( iConfig.getUntrackedParameter<string>( "genParticleCandidates" ) ),
  numTkCut( iConfig.getUntrackedParameter<unsigned int>( "TracksMinimumNum_in_PixelJet" ) ),
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

  // Pixel trigger efficiency
  EffMultijetPixelSize_ = 10000;
  EffMEtJetPixelSize_ = 100;
  TString EffMultijetName("EffMultijetPixel");
  TString EffMultijetTitle("Efficiency of multijet + pixel trigger for PVnum = ");
  TString EffMEtJetName("EffMEtJetPixel");
  TString EffMEtJetTitle("Efficiency of MEt + Jet + pixel trigger for PVnum = ");
  TString PVnumString[4];
  PVnumString[0] = "3";
  PVnumString[1] = "4";
  PVnumString[2] = "5";
  PVnumString[3] = "6";
  EffMultijetPixelArray_ = new int*[4];
  EffMEtJetPixelArray_ = new int*[4];
  EffMultijetPixel_ = new TH1F*[4];
  EffMEtJetPixel_ = new TH1F*[4];
  for ( int PVnum=0; PVnum<4; ++PVnum ) {
    EffMultijetPixelArray_[PVnum] = new int[EffMultijetPixelSize_];
    EffMEtJetPixelArray_[PVnum]= new int[EffMEtJetPixelSize_];
    // Initialize the efficiency counters
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSize_; ++iEffMultijetNum ) {
      (EffMultijetPixelArray_[PVnum])[iEffMultijetNum] = 0;
    }
    for ( int iEffMEtJetNum=0; iEffMEtJetNum<EffMEtJetPixelSize_; ++iEffMEtJetNum ) {
      (EffMEtJetPixelArray_[PVnum])[iEffMEtJetNum] = 0;
    }

    EffMultijetPixel_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum], EffMultijetTitle + PVnumString[PVnum], EffMultijetPixelSize_, 0, EffMultijetPixelSize_ );
    EffMEtJetPixel_[PVnum] = new TH1F( EffMEtJetName + PVnumString[PVnum], EffMEtJetTitle + PVnumString[PVnum], EffMEtJetPixelSize_, 0, EffMEtJetPixelSize_ );
  }

  // Generate histograms for the verteces
  // ------------------------------------
  dz_ = 0.04;
  bins_ = 20;

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

  // Trigger efficiency counters
  // Multijet
  Eff_ = 0;
  Eff_cen_ = 0;
  Eff_tau_ = 0;
  Eff_for_ = 0;
  Eff_nofor_ = 0;
  // MEt+Jet
  Eff_MEtJet_ = 0;
  Eff_MEtJet_cen_ = 0;
  Eff_MEtJet_tau_ = 0;
  Eff_MEtJet_for_ = 0;
  Eff_MEtJet_nofor_ = 0;
  // Tau
  Eff_tautrig_ = 0;

  eventcounter_ = 0;
//  PI_ = 3.141593;
}


L1TrigPixelAnalyzer::~L1TrigPixelAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // Draw the histograms
  //  HiVar->Plot();

  OutputFile->Write();
  for (int i=0; i<4; ++i) {
    delete[] EffMultijetPixelArray_[i];
    delete[] EffMEtJetPixelArray_[i];
  }
  delete[] EffMultijetPixelArray_;
  delete[] EffMEtJetPixelArray_;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
L1TrigPixelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace l1extra;
  using namespace std;

  // L1 Calo
  edm::Handle < L1JetParticleCollection > l1eCenJets;
  edm::Handle < L1JetParticleCollection > l1eForJets;
  edm::Handle < L1JetParticleCollection > l1eTauJets;
  edm::Handle < L1EtMissParticle > l1eEtMiss;


  // should get rid of this try/catch?
  try {
    edm::InputTag L1CJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eCentralJetsSource");
    edm::InputTag L1FJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eForwardJetsSource");
    edm::InputTag L1TauJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eTauJetsSource");
    edm::InputTag L1EtMissLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eEtMissSource");

    iEvent.getByLabel(L1CJetLabel, l1eCenJets);
    iEvent.getByLabel(L1FJetLabel, l1eForJets);
    iEvent.getByLabel(L1TauJetLabel, l1eTauJets);
    iEvent.getByLabel(L1EtMissLabel, l1eEtMiss);

  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }



  // Pixel jets
  // ----------

  edm::Handle<PixelJetCollection> pixeljetshandle;
  edm::InputTag PixelJetsLabel = conf_.getUntrackedParameter<edm::InputTag>("PixelJetSource");
  iEvent.getByLabel(PixelJetsLabel, pixeljetshandle);

  const PixelJetCollection pixeljets = *(pixeljetshandle.product());

  // Pixel trigger requiring at least 2 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)

  for( int numdz = 0; numdz < bins_; ++numdz ) {

    double dz = dz_ + dz_*(numdz);

#ifdef DEBUG
    std::cout << "numTkCut = " << numTkCut << std::endl;
#endif

    L1PixelTrig<PixelJet> PJtrig(2, dz, numTkCut);

    PJtrig.Fill( pixeljets );

#ifdef DEBUG
    std::cout << "Pixel trigger response = " << PJtrig.Response() << std::endl;
#endif

    // The vector of verteces is already sorted in ascending Pt
    auto_ptr<vector<Vertex<PixelJet> > > vec_vertex_aptr( PJtrig.VevVertex() );

    // Number of verteces
    int numVertex = 0;
    if ( vec_vertex_aptr.get() != 0 ) {
      numVertex = vec_vertex_aptr->size();
    }
    Multi_Vertex_Num_->Fill( numVertex, numdz );

#ifdef DEBUG
    std::cout << "Loop on primary vertex" << std::endl;
#endif
    if ( numVertex > 0 ) {
      //Primary vertex characteristics
      Vertex<PixelJet> prim_vertex( vec_vertex_aptr->back() );

      Multi_PrimVNum_->Fill( prim_vertex.number(), numdz );
      Multi_PrimVPt_->Fill( prim_vertex.pt(), numdz );
      Multi_PrimVEta_->Fill( prim_vertex.eta(), numdz );
      Multi_PrimVPhi_->Fill( prim_vertex.phi(), numdz );

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

#ifdef DEBUG
      std::cout << "Verteces minimum distance in z" << std::endl;
#endif

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
  }

#ifdef DEBUG
  std::cout << "PixelJet minimum distance in z" << std::endl;
#endif

  // Take another pixeljet vector and evaluate the minimum dz distance
  const PixelJetCollection temp_pixeljets = *(pixeljetshandle.product());
  int numPixelJet = temp_pixeljets.size();
  PixelJet_Num_->Fill( numPixelJet );
  if ( numPixelJet > 1 ) {
    vector<PixelJet>::const_iterator pj_z_it (temp_pixeljets.begin()+1);
    PixelJet_Track_Num_->Fill(temp_pixeljets.begin()->NumTk());
    for ( ; pj_z_it != temp_pixeljets.end(); ++pj_z_it ) {
      PixelJet_Track_Num_->Fill(pj_z_it->NumTk());
      PixelJet_dz_->Fill( fabs(pj_z_it->z() - (pj_z_it-1)->z()) );
    }
  }

  // Level 1 trigger
  // ---------------

  // All the jets together fot the L1Trigger
  vector<SimpleJet> vec_TriggerCenJet;
  vector<SimpleJet> vec_TriggerForJet;
  vector<SimpleJet> vec_TriggerTauJet;
  for ( L1JetParticleCollection::const_iterator tcj = l1eCenJets->begin(); tcj != l1eCenJets->end(); ++tcj ) {
    vec_TriggerCenJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
  }
  int fjcount = 0;
  for ( L1JetParticleCollection::const_iterator tfj = l1eForJets->begin(); tfj != l1eForJets->end(); ++tfj ) {
    vec_TriggerForJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
    std::cout << "ForwardJet Et["<<fjcount<<"] = " << tfj->et() << std::endl;
    std::cout << "ForwardJet Eta["<<fjcount<<"] = " << tfj->eta() << std::endl;
    std::cout << "ForwardJet Phi["<<fjcount<<"] = " << tfj->phi() << std::endl;
    ++fjcount;
  }
  // Tau jets
  for ( L1JetParticleCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
    vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
  }

  // Multijet
  // --------
  // Central
  bool response_cen = false;
  L1Trigger.Fill( vec_TriggerCenJet );
  response_cen = L1Trigger.Response();
  // Forward
  bool response_for = false;
  L1Trigger.Fill( vec_TriggerForJet );
  response_for = L1Trigger.Response();
  // Tau
  bool response_tau = false;
  L1Trigger.Fill( vec_TriggerTauJet );
  response_tau = L1Trigger.Response();
  // Full and no-forward
  bool response = ( response_cen || response_tau || response_for );
  bool response_nofor = ( response_cen || response_tau );

  // MEt + Jet
  // ---------
  // Central
  bool response_MEtJet_cen = false;
  sort( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
  reverse( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
  if ( vec_TriggerCenJet.size() != 0 ) {
    if ( (vec_TriggerCenJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_cen = true;
    }
  }
  // Tau
  bool response_MEtJet_tau = false;
  sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  if ( vec_TriggerTauJet.size() != 0 ) {
    if ( (vec_TriggerTauJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_tau = true;
    }
  }
  // Forward
  bool response_MEtJet_for = false;
  sort( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
  reverse( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
  if ( vec_TriggerForJet.size() != 0 ) {
    if ( (vec_TriggerForJet[0].pt() >= 80.) && (l1eEtMiss->et() >= 100.) ) {
      response_MEtJet_for = true;
    }
  }
  // Full and no-forward
  bool response_MEtJet = ( response_MEtJet_cen || response_MEtJet_tau || response_MEtJet_for );
  bool response_MEtJet_nofor = ( response_MEtJet_cen || response_MEtJet_tau );

  // Tau trigger
  // -----------
  bool response_tautrig = false;
  // Already sorted for the previous case
//   sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
//   reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );

  // Single tau trigger
  if ( vec_TriggerTauJet.size() != 0 ) {
    if ( vec_TriggerTauJet[0].pt() >= 150. ) {
      response_tautrig = true;
    }
  }
  // Di-tau trigger
  if ( vec_TriggerTauJet.size() >= 1 ) {
    if ( vec_TriggerTauJet[1].pt() >= 80. ) {
      response_tautrig = true;
    }
  }

#ifdef DEBUG
  std::cout << "L1Trigger response cen = " << response_cen << std::endl;
  std::cout << "L1Trigger response for = " << response_for << std::endl;
  std::cout << "L1Trigger response tau = " << response_tau << std::endl;
  std::cout << "L1Trigger response = " << response << std::endl;

  std::cout << "Number of L1Jets central = " << vec_TriggerJetCen.size() << std::endl;
  std::cout << "Number of L1Jets forward = " << vec_TriggerJetFor.size() << std::endl;
  std::cout << "Number of L1Jets tau = " << vec_TriggerJetTau.size() << std::endl;
#endif

  // Count the trigger efficiency
  // ----------------------------

  // Multijet
  if ( response_cen ) {
    ++Eff_cen_;
  }
  if ( response_tau ) {
    ++Eff_tau_;
  }
  if ( response_for ) {
    ++Eff_for_;
  }
  if ( response ) {
    ++Eff_;
  }
  if ( response_nofor ) {
    ++Eff_nofor_;
  }
  // MEt+jet
  if ( response_MEtJet_cen ) {
    ++Eff_MEtJet_cen_;
  }
  if ( response_MEtJet_tau ) {
    ++Eff_MEtJet_tau_;
  }
  if ( response_MEtJet_for ) {
    ++Eff_MEtJet_for_;
  }
  if ( response_MEtJet ) {
    ++Eff_MEtJet_;
  }
  if ( response_MEtJet_nofor ) {
    ++Eff_MEtJet_nofor_;
  }

  // Tau trigger
  if ( response_tautrig ) {
    ++Eff_tautrig_;
  }

  // PixelTrigger
  // ------------
  // dz = 0.4 and numtk = 3
  L1PixelTrig<PixelJet> Pixeltrig_3(3, 0.4, 3);
  L1PixelTrig<PixelJet> Pixeltrig_4(4, 0.4, 3);
  L1PixelTrig<PixelJet> Pixeltrig_5(5, 0.4, 3);
  L1PixelTrig<PixelJet> Pixeltrig_6(6, 0.4, 3);

  Pixeltrig_3.Fill( pixeljets );
  Pixeltrig_4.Fill( pixeljets );
  Pixeltrig_5.Fill( pixeljets );
  Pixeltrig_6.Fill( pixeljets );

  bool pixelTrigResponse[4];

  pixelTrigResponse[0] = Pixeltrig_3.Response();
  pixelTrigResponse[1] = Pixeltrig_4.Response();
  pixelTrigResponse[2] = Pixeltrig_5.Response();
  pixelTrigResponse[3] = Pixeltrig_6.Response();

  // Multijet initializations
  // ------------------------
  bool Response_cen = false;
  bool Response_tau = false;

  sort( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );
  reverse( vec_TriggerCenJet.begin(), vec_TriggerCenJet.end() );

  float cenJet1 = 0.;
  float cenJet2 = 0.;
  float cenJet3 = 0.;
  float cenJet4 = 0.;
  int cenSize = vec_TriggerCenJet.size();
  if ( cenSize != 0 ) cenJet1 = vec_TriggerCenJet[0].pt();
  if ( cenSize > 1 ) cenJet2 = vec_TriggerCenJet[1].pt();
  if ( cenSize > 2 ) cenJet3 = vec_TriggerCenJet[2].pt();
  if ( cenSize > 3 ) cenJet4 = vec_TriggerCenJet[3].pt();

  sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
  reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );

  float tauJet1 = 0.;
  float tauJet2 = 0.;
  float tauJet3 = 0.;
  float tauJet4 = 0.;
  int tauSize = vec_TriggerTauJet.size();
  if ( tauSize != 0 ) tauJet1 = vec_TriggerTauJet[0].pt();
  if ( tauSize > 1 ) tauJet2 = vec_TriggerTauJet[1].pt();
  if ( tauSize > 2 ) tauJet3 = vec_TriggerTauJet[2].pt();
  if ( tauSize > 3 ) tauJet4 = vec_TriggerTauJet[3].pt();

  // MEt + Jet initializations
  // -------------------------
  bool Response_MEtJet_cen = false;
  bool Response_MEtJet_tau = false;
  float L1MEt = l1eEtMiss->et();

  // Loop on all the values
  // ----------------------

  // Do it for the different PV cut values for the pixel trigger
  for ( int PVnum=0; PVnum<4; ++PVnum ) {

    // Loop on all the cut values
    int multijetCount = 0;
    for ( int Et1=160; Et1<260; Et1+=10 ) {
      for ( int Et2=110; Et2<210; Et2+=10 ) {
        for ( int Et3=55; Et3<105; Et3+=5 ) {
          for ( int Et4=35; Et4<85; Et4+=5 ) {
            Response_cen = false;
            Response_tau = false;
            // Central
            if ( (cenJet1 >= Et1) || (cenJet2 >= Et2) || (cenJet3 >= Et3) || (cenJet4 >= Et4) ) Response_cen = true;
            // Tau
            if ( (tauJet1 >= Et1) || (tauJet2 >= Et2) || (tauJet3 >= Et3) || (tauJet4 >= Et4) ) Response_tau = true;
            // Full and no-forward
            if ( (Response_cen || Response_tau) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArray_[PVnum])[multijetCount]);
            // increase the index, in the inner loop
            ++multijetCount;
          }
        }
      }
    }

    // MEt + Jet
    // ---------
    int mEtJetCount = 0;
    for ( int MEt=35; MEt<85; MEt+=5 ) {
      for ( int Jet=55; Jet<105; Jet+=5 ) {
        Response_MEtJet_cen = false;
        Response_MEtJet_tau = false;
        // Central
        if ( ( L1MEt >= float(MEt) ) && ( cenJet1 >= float(Jet) ) ) Response_MEtJet_cen = true;
        // Tau
        if ( ( L1MEt >= float(MEt) ) && ( tauJet1 >= float(Jet) ) ) Response_MEtJet_tau = true;
        // Full and no-forward
        if ( (Response_MEtJet_cen || Response_MEtJet_tau) && pixelTrigResponse[PVnum] ) ++((EffMEtJetPixelArray_[PVnum])[mEtJetCount]);
        ++mEtJetCount;
      }
    }
  } // end loop on PV cuts

  // ------- //
  // Offline //
  // ------- //

  // MEt
  // ---

  edm::Handle<reco::CaloMETCollection> caloMET;
  iEvent.getByLabel( METCollection, caloMET );

  const reco::CaloMET * MET = &( *(caloMET->begin()) );

  MEt_CorrIC5_Pt_->Fill( MET->pt() );
  double MET_phi = MET->phi();
  MEt_CorrIC5_Phi_->Fill( MET_phi );
  MEt_CorrIC5_SumEt_->Fill( MET->sumEt() );
  MEt_CorrIC5_mEtSig_->Fill( MET->mEtSig() );

//  Associator<reco::CaloMET, reco::CaloJet> associator( 0.5 );
//  std::auto_ptr<std::map<const reco::CaloMET*, const reco::CaloJet*> > AssocMap( associator.Associate( *caloMET, *caloJets ) );

  // Take the first (and it should be the only one) element in the map and get the eta of the closest jet (in DeltaR)
//   double JetPhi = (*AssocMap).begin()->second->phi();

  // Do not use DR association. The MEt has no z component. Use DeltaPhi association.

  // HiVariables
  // -----------

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel( CaloJetAlgorithm, caloJets );

  if ( caloJets->size() != 0 ) {

    std::vector<double> vec_DPhi;
    vec_DPhi.reserve(caloJets->size());

    vector<SimpleJet> vec_calojet; 
    vec_calojet.reserve(caloJets->size());
    // Correct offline jets on the fly
    const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService, iSetup);
    for( reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal ) {
      double scale = corrector->correction( *cal );
      double corPt = scale*cal->pt();
      vec_calojet.push_back( SimpleJet( corPt, cal->eta(), cal->phi() ) );
      uncorr_JetPt_IC5_->Fill( cal->pt() );
      corr_JetPt_IC5_->Fill( corPt );

      // Evaluate DeltaPhi between MET and calo-jets
      // Consider only high-luminosity--good-jets
      if ( corPt > 40 && fabs( cal->eta() ) < 3.0 ) {
        vec_DPhi.push_back( DeltaPhi( MET_phi, cal->phi() ) );
      }
    }

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
  edm::Handle < CandidateCollection > MCpartons;
  iEvent.getByLabel( genParticleCandidates, MCpartons );

  // Take the mothers
  CandidateCollection::const_iterator Mothers = MCpartons->begin() + 6;
//  std::cout << "First of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
  const Candidate* mother_1 = &*Mothers;
  ++Mothers;
//  std::cout << "Second of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
  const Candidate* mother_2 = &*Mothers;
  ++Mothers;
  // Take the other mother only if it is a Higgs
  const Candidate* mother_3 = 0;
  if ( Mothers->pdgId() == 25 ) {
//    std::cout << "Third of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
    mother_3 = &*Mothers;
  }

  std::vector<const Candidate*> vec_Partons;
  int counter = 0;
  ofstream Partons( "Partons.txt", ios::app );
  Partons << std::endl;
  int nu_count = 0;
  for( CandidateCollection::const_iterator MCp = MCpartons->begin()+8; MCp != MCpartons->end(); ++MCp, ++counter ) {
    const Candidate* temp_mother = MCp->mother();
    if ( temp_mother != 0 ) {

      // Store the status = 3 partons
      if ( MCp->status() == 3 ) {
        vec_Partons.push_back( &*MCp );
      }

      int pid = abs( MCp->pdgId() );
      int Mpid = abs( MCp->mother()->pdgId() );
      if ( ( pid == 12 || pid == 14 || pid == 16 ) && ( Mpid == 24 ) ) {
        ++nu_count;
        std::cout << "Neutrino from W(from top) number " << nu_count << std::endl;;
      }

      if ( MCp->status() == 3 && ( temp_mother == mother_1 || temp_mother == mother_2 || temp_mother == mother_3 ) ){
#ifdef DEBUG
        std::cout << "For parton number = " << counter << std::endl;
        std::cout << "status = " << MCp->status() << std::endl;
        std::cout << "pdgId = " << MCp->pdgId() << std::endl;
        std::cout << "Et = " << MCp->et() << std::endl;
        std::cout << "Eta = " << MCp->eta() << std::endl;
        std::cout << "Phi = " << MCp->phi() << std::endl;
        std::cout << "Number of mothers = " << MCp->numberOfMothers() << std::endl;
        std::cout << "first mother = " << MCp->mother() << std::endl;
        std::cout << "Mother pdgId = " << MCp->mother()->pdgId() << std::endl;

#endif // DEBUG
        Partons << "For parton number = " << counter << std::endl;
        Partons << "status = " << MCp->status() << std::endl;
        Partons << "pdgId = " << MCp->pdgId() << std::endl;
        Partons << "Et = " << MCp->et() << std::endl;
        Partons << "Eta = " << MCp->eta() << std::endl;
        Partons << "Phi = " << MCp->phi() << std::endl;
        Partons << "Number of mothers = " << MCp->numberOfMothers() << std::endl;
        Partons << "first mother = " << MCp->mother() << std::endl;
        Partons << "Mother pdgId = " << MCp->mother()->pdgId() << std::endl;
      }
    }
  }
  Partons.close();

  // Reverse the Partons, this we start from the last and we can exclude the mothers
  reverse( vec_Partons.begin(), vec_Partons.end() );

  vector<const Candidate*> vec_Mothers;
  for( vector<const Candidate*>::const_iterator Par_it = vec_Partons.begin(); Par_it != vec_Partons.end(); ++Par_it ) {

    vec_Mothers.push_back( (*Par_it)->mother() );

  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void L1TrigPixelAnalyzer::beginJob(const edm::EventSetup&) {
}

// void L1TrigPixelAnalyzer::meanHistoSetup( TH1F * meanhisto_ptr, vector<TH1F *> vec_histo ) {
//   meanhisto_ptr->SetBinContent( numdz+1, mean[numdz]->GetMean() );
//   meanhisto_ptr->SetBinError( numdz+1, Vertex_Dz_[numdz]->GetMeanError() );
//   meanhisto_ptr->GetXaxis()->CenterLabels();
//   meanhisto_ptr->GetXaxis()->SetNdivisions(Vertex_Dz_Mean_->GetSize()-2, false);
// }

// ------------ method called once each job just after ending the event loop  ------------
void L1TrigPixelAnalyzer::endJob() {

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

  // Pixel trigger efficiency
  ofstream multijetPixelEffFile( "multijetPixelEff.txt" );
  ofstream mEtJetPixelEffFile( "mEtJetPixelEff.txt" );
  for ( int PVnum=0; PVnum<4; ++PVnum ) {
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSize_; ++multijetCount ) {
      float multijetPixelEff = float((EffMultijetPixelArray_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErr = sqrt(float((EffMultijetPixelArray_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixel_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEff );
      EffMultijetPixel_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErr );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff = " << multijetPixelEff << " +- " << multijetPixelEffErr << endl;
    }
    for ( int mEtJetCount=0; mEtJetCount<EffMEtJetPixelSize_; ++mEtJetCount ) {
      float mEtJetPixelEff = float((EffMEtJetPixelArray_[PVnum])[mEtJetCount])/float(eventcounter_);
      float mEtJetPixelEffErr = sqrt(float((EffMEtJetPixelArray_[PVnum])[mEtJetCount]))/float(eventcounter_);
      EffMEtJetPixel_[PVnum]->SetBinContent( mEtJetCount+1, mEtJetPixelEff );
      EffMEtJetPixel_[PVnum]->SetBinError( mEtJetCount+1, mEtJetPixelEffErr );
      mEtJetPixelEffFile << "For PVnum = "<<PVnum<<" and i = "<< mEtJetCount <<" Eff = " << mEtJetPixelEff << " +- " << mEtJetPixelEffErr << endl;
    }
  }

  ofstream Effoutputfile( OutputEffFileName.c_str() );

  Effoutputfile << "Multijet trigger" << endl;
  Effoutputfile << "----------------" << endl;
  Effoutputfile << "Eff multijet = " << float(Eff_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff multijet no-forward = " <<  float(Eff_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet = " << float(Eff_cen_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet = " << float(Eff_tau_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff for jet = " << float(Eff_for_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "MEt + Jet trigger" << endl;
  Effoutputfile << "-----------------" << endl;
  Effoutputfile << "Eff MEt + Jet = " << float(Eff_MEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + Jet no-forward = " << float(Eff_MEtJet_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + central jet = " << float(Eff_MEtJet_cen_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + tau jet = " << float(Eff_MEtJet_tau_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff MEt + forward jet = " << float(Eff_MEtJet_for_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Tau trigger" << endl;
  Effoutputfile << "-----------" << endl;
  Effoutputfile << "Eff tau trigger = " << float(Eff_tautrig_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Total events = " << eventcounter_ << endl;
  Effoutputfile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrigPixelAnalyzer);
