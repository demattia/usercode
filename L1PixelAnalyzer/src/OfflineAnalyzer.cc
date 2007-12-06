//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/OfflineAnalyzer.h"

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

L1Trig OfflineAnalyzer::L1Trigger;

//
// constructors and destructor
//
OfflineAnalyzer::OfflineAnalyzer(const edm::ParameterSet& iConfig) :
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
  EffMultijetPixelSizeEt1_ = 10;
  EffMultijetPixelSizeEt2_ = 10;
  EffMultijetPixelSizeEt3_ = 10;
  EffMultijetPixelSizeEt4_ = 10;
  offlineEffMultijetPixelSize_ = 10000;
  offlineEffMEtJetPixelSize_ = 100;
  TString EffMultijetName("EffMultijetPixel");
  TString EffMultijetTitle("Efficiency of multijet + pixel trigger for PVnum = ");
  TString EffMEtJetName("EffMEtJetPixel");
  TString EffMEtJetTitle("Efficiency of MEt + Jet + pixel trigger for PVnum = ");
  TString offlineEffMultijetName("offlineEffMultijetPixel");
  TString offlineEffMultijetTitle("Efficiency of multijet + pixel trigger + offline for PVnum = ");
  TString offlineEffMEtJetName("offlineEffMEtJetPixel");
  TString offlineEffMEtJetTitle("Efficiency of MEt + Jet + pixel + offline trigger for PVnum = ");
  TString PVnumString[4];
  PVnumString[0] = "3";
  PVnumString[1] = "4";
  PVnumString[2] = "5";
  PVnumString[3] = "6";
  TString Et1String("Et1");
  TString Et2String("Et2");
  TString Et3String("Et3");
  TString Et4String("Et4");
  EffMultijetPixelArray_ = new int*[4];
  EffMultijetPixelArrayEt1_ = new int*[4];
  EffMultijetPixelArrayEt2_ = new int*[4];
  EffMultijetPixelArrayEt3_ = new int*[4];
  EffMultijetPixelArrayEt4_ = new int*[4];
  EffMEtJetPixelArray_ = new int*[4];
  EffMultijetPixel_ = new TH1F*[4];
  EffMultijetPixelEt1_ = new TH1F*[4];
  EffMultijetPixelEt2_ = new TH1F*[4];
  EffMultijetPixelEt3_ = new TH1F*[4];
  EffMultijetPixelEt4_ = new TH1F*[4];
  EffMEtJetPixel_ = new TH1F*[4];
  // Offline
  offlineEffMultijetPixel_ = new TH1F*[4];
  offlineEffMultijetPixelArray_ = new int*[4];
  offlineEffMEtJetPixelArray_ = new int*[4];
  offlineEffMEtJetPixel_ = new TH1F*[4];
  for ( int PVnum=0; PVnum<4; ++PVnum ) {
    EffMultijetPixelArray_[PVnum] = new int[EffMultijetPixelSize_];
    EffMultijetPixelArrayEt1_[PVnum] = new int[EffMultijetPixelSizeEt1_];
    EffMultijetPixelArrayEt2_[PVnum] = new int[EffMultijetPixelSizeEt2_];
    EffMultijetPixelArrayEt3_[PVnum] = new int[EffMultijetPixelSizeEt3_];
    EffMultijetPixelArrayEt4_[PVnum] = new int[EffMultijetPixelSizeEt4_];

    EffMEtJetPixelArray_[PVnum]= new int[EffMEtJetPixelSize_];

    offlineEffMultijetPixelArray_[PVnum] = new int[EffMultijetPixelSize_];
    offlineEffMEtJetPixelArray_[PVnum]= new int[EffMEtJetPixelSize_];

    // Initialize the efficiency counters
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSize_; ++iEffMultijetNum ) {
      (EffMultijetPixelArray_[PVnum])[iEffMultijetNum] = 0;
      (offlineEffMultijetPixelArray_[PVnum])[iEffMultijetNum] = 0;
    }
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSizeEt1_; ++iEffMultijetNum ) {
      (EffMultijetPixelArrayEt1_[PVnum])[iEffMultijetNum] = 0;
    }
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSizeEt2_; ++iEffMultijetNum ) {
      (EffMultijetPixelArrayEt2_[PVnum])[iEffMultijetNum] = 0;
    }
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSizeEt3_; ++iEffMultijetNum ) {
      (EffMultijetPixelArrayEt3_[PVnum])[iEffMultijetNum] = 0;
    }
    for ( int iEffMultijetNum=0; iEffMultijetNum<EffMultijetPixelSizeEt4_; ++iEffMultijetNum ) {
      (EffMultijetPixelArrayEt4_[PVnum])[iEffMultijetNum] = 0;
    }

    for ( int iEffMEtJetNum=0; iEffMEtJetNum<EffMEtJetPixelSize_; ++iEffMEtJetNum ) {
      (EffMEtJetPixelArray_[PVnum])[iEffMEtJetNum] = 0;
      (offlineEffMEtJetPixelArray_[PVnum])[iEffMEtJetNum] = 0;
    }

    EffMultijetPixel_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum], EffMultijetTitle + PVnumString[PVnum], EffMultijetPixelSize_, 0, EffMultijetPixelSize_ );
    EffMultijetPixelEt1_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum] + Et1String, EffMultijetTitle + PVnumString[PVnum] + Et1String, EffMultijetPixelSizeEt1_, 0, EffMultijetPixelSizeEt1_ );
    EffMultijetPixelEt2_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum] + Et2String, EffMultijetTitle + PVnumString[PVnum] + Et2String, EffMultijetPixelSizeEt2_, 0, EffMultijetPixelSizeEt2_ );
    EffMultijetPixelEt3_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum] + Et3String, EffMultijetTitle + PVnumString[PVnum] + Et3String, EffMultijetPixelSizeEt3_, 0, EffMultijetPixelSizeEt3_ );
    EffMultijetPixelEt4_[PVnum] = new TH1F( EffMultijetName + PVnumString[PVnum] + Et4String, EffMultijetTitle + PVnumString[PVnum] + Et4String, EffMultijetPixelSizeEt4_, 0, EffMultijetPixelSizeEt4_ );

    EffMEtJetPixel_[PVnum] = new TH1F( EffMEtJetName + PVnumString[PVnum], EffMEtJetTitle + PVnumString[PVnum], EffMEtJetPixelSize_, 0, EffMEtJetPixelSize_ );

    offlineEffMultijetPixel_[PVnum] = new TH1F( offlineEffMultijetName + PVnumString[PVnum], offlineEffMultijetTitle + PVnumString[PVnum], EffMultijetPixelSize_, 0, EffMultijetPixelSize_ );
    offlineEffMEtJetPixel_[PVnum] = new TH1F( offlineEffMEtJetName + PVnumString[PVnum], offlineEffMEtJetTitle + PVnumString[PVnum], EffMEtJetPixelSize_, 0, EffMEtJetPixelSize_ );
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
  Eff_et1_ = 0;
  Eff_et2_ = 0;
  Eff_et3_ = 0;
  Eff_et4_ = 0;
  Eff_cen_ = 0;
  Eff_cen_et1_ = 0;
  Eff_cen_et2_ = 0;
  Eff_cen_et3_ = 0;
  Eff_cen_et4_ = 0;
  Eff_tau_ = 0;
  Eff_tau_et1_ = 0;
  Eff_tau_et2_ = 0;
  Eff_tau_et3_ = 0;
  Eff_tau_et4_ = 0;
  Eff_for_ = 0;
  Eff_for_et1_ = 0;
  Eff_for_et2_ = 0;
  Eff_for_et3_ = 0;
  Eff_for_et4_ = 0;
  Eff_nofor_ = 0;
  Eff_nofor_et1_ = 0;
  Eff_nofor_et2_ = 0;
  Eff_nofor_et3_ = 0;
  Eff_nofor_et4_ = 0;
  // MEt+Jet
  Eff_MEtJet_ = 0;
  Eff_MEtJet_cen_ = 0;
  Eff_MEtJet_tau_ = 0;
  Eff_MEtJet_for_ = 0;
  Eff_MEtJet_nofor_ = 0;
  // Tau
  Eff_tautrig_ = 0;
  Eff_tautrig_single_ = 0;
  Eff_tautrig_ditau_ = 0;

  // Offline
  offlineEffMultijet_ = 0;
  offlineEffMEtJet_ = 0;
  offlineEffTauTrig_ = 0;

  eventcounter_ = 0;
//  PI_ = 3.141593;
}


OfflineAnalyzer::~OfflineAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // Draw the histograms
  //  HiVar->Plot();

  OutputFile->Write();
  for (int i=0; i<4; ++i) {
    delete[] EffMultijetPixelArray_[i];
    delete[] EffMultijetPixelArrayEt1_[i];
    delete[] EffMultijetPixelArrayEt2_[i];
    delete[] EffMultijetPixelArrayEt3_[i];
    delete[] EffMultijetPixelArrayEt4_[i];
    delete[] EffMEtJetPixelArray_[i];
    delete[] offlineEffMultijetPixelArray_[i];
    delete[] offlineEffMEtJetPixelArray_[i];
  }
  delete[] EffMultijetPixelArray_;
  delete[] EffMultijetPixelArrayEt1_;
  delete[] EffMultijetPixelArrayEt2_;
  delete[] EffMultijetPixelArrayEt3_;
  delete[] EffMultijetPixelArrayEt4_;
  delete[] EffMEtJetPixelArray_;
  delete[] offlineEffMultijetPixelArray_;
  delete[] offlineEffMEtJetPixelArray_;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
OfflineAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  MCParticleCollection::const_iterator MCp = MCpartons->begin();
  for ( ; MCp != MCpartons->end(); ++MCp ) {

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
  }

  // Pixel jets
  // ----------

  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );

  const SimplePixelJetCollection pixeljets = *(pixelJetsHandle.product());

  // Pixel trigger requiring at least 2 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)

  for( int numdz = 0; numdz < bins_; ++numdz ) {

    double dz = dz_ + dz_*(numdz);

#ifdef DEBUG
    std::cout << "numTkCut = " << numTkCut << std::endl;
#endif

    L1PixelTrig<SimplePixelJet> PJtrig(2, dz, numTkCut);

    PJtrig.Fill( pixeljets );

#ifdef DEBUG
    std::cout << "Pixel trigger response = " << PJtrig.Response() << std::endl;
#endif

    // The vector of verteces is already sorted in ascending Pt
    auto_ptr<vector<Vertex<SimplePixelJet> > > vec_vertex_aptr( PJtrig.VevVertex() );

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

#ifdef DEBUG
      std::cout << "Verteces minimum distance in z" << std::endl;
#endif

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
  }

#ifdef DEBUG
  std::cout << "PixelJet minimum distance in z" << std::endl;
#endif

  // Take another pixeljet vector and evaluate the minimum dz distance
  const SimplePixelJetCollection temp_pixeljets = *(pixelJetsHandle.product());
  int numPixelJet = temp_pixeljets.size();
  PixelJet_Num_->Fill( numPixelJet );
  if ( numPixelJet > 1 ) {
    vector<SimplePixelJet>::const_iterator pj_z_it (temp_pixeljets.begin()+1);
    PixelJet_Track_Num_->Fill(temp_pixeljets.begin()->tkNum());
    for ( ; pj_z_it != temp_pixeljets.end(); ++pj_z_it ) {
      PixelJet_Track_Num_->Fill(pj_z_it->tkNum());
      PixelJet_dz_->Fill( fabs(pj_z_it->z() - (pj_z_it-1)->z()) );
    }
  }

  // Level 1 trigger
  // ---------------

  // All the jets together fot the L1Trigger
  vector<SimpleJet> vec_TriggerCenJet;
  vector<SimpleJet> vec_TriggerForJet;
  vector<SimpleJet> vec_TriggerTauJet;
  for ( BaseJetCollection::const_iterator tcj = l1eCenJets->begin(); tcj != l1eCenJets->end(); ++tcj ) {
    vec_TriggerCenJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
  }
  int fjcount = 0;
  for ( BaseJetCollection::const_iterator tfj = l1eForJets->begin(); tfj != l1eForJets->end(); ++tfj ) {
    vec_TriggerForJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
#ifdef DEBUG
    std::cout << "ForwardJet Et["<<fjcount<<"] = " << tfj->et() << std::endl;
    std::cout << "ForwardJet Eta["<<fjcount<<"] = " << tfj->eta() << std::endl;
    std::cout << "ForwardJet Phi["<<fjcount<<"] = " << tfj->phi() << std::endl;
#endif
    ++fjcount;
  }
  // Tau jets
  for ( BaseJetCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
    vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
  }

  // Multijet
  // --------
  // Central
  bool response_cen = false;
  bool response_cen_et1 = false;
  bool response_cen_et2 = false;
  bool response_cen_et3 = false;
  bool response_cen_et4 = false;
  L1Trigger.Fill( vec_TriggerCenJet );
  response_cen = L1Trigger.Response();
  response_cen_et1 = L1Trigger.Et1Cut();
  response_cen_et2 = L1Trigger.Et2Cut();
  response_cen_et3 = L1Trigger.Et3Cut();
  response_cen_et4 = L1Trigger.Et4Cut();
  // Forward
  bool response_for = false;
  bool response_for_et1 = false;
  bool response_for_et2 = false;
  bool response_for_et3 = false;
  bool response_for_et4 = false;
  L1Trigger.Fill( vec_TriggerForJet );
  response_for = L1Trigger.Response();
  response_for_et1 = L1Trigger.Et1Cut();
  response_for_et2 = L1Trigger.Et2Cut();
  response_for_et3 = L1Trigger.Et3Cut();
  response_for_et4 = L1Trigger.Et4Cut();
  // Tau
  bool response_tau = false;
  bool response_tau_et1 = false;
  bool response_tau_et2 = false;
  bool response_tau_et3 = false;
  bool response_tau_et4 = false;
  L1Trigger.Fill( vec_TriggerTauJet );
  response_tau = L1Trigger.Response();
  response_tau_et1 = L1Trigger.Et1Cut();
  response_tau_et2 = L1Trigger.Et2Cut();
  response_tau_et3 = L1Trigger.Et3Cut();
  response_tau_et4 = L1Trigger.Et4Cut();
  // Full and no-forward
  bool response = ( response_cen || response_tau || response_for );
  bool response_Et1 = ( response_cen_et1 || response_tau_et1 || response_for_et1 );
  bool response_Et2 = ( response_cen_et2 || response_tau_et2 || response_for_et2 );
  bool response_Et3 = ( response_cen_et3 || response_tau_et3 || response_for_et3 );
  bool response_Et4 = ( response_cen_et4 || response_tau_et4 || response_for_et4 );
  bool response_nofor = ( response_cen || response_tau );
  bool response_Et1_nofor = ( response_cen_et1 || response_tau_et1 );
  bool response_Et2_nofor = ( response_cen_et2 || response_tau_et2 );
  bool response_Et3_nofor = ( response_cen_et3 || response_tau_et3 );
  bool response_Et4_nofor = ( response_cen_et4 || response_tau_et4 );

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
  bool response_tautrig_single = false;
  bool response_tautrig_ditau = false;
  // Already sorted for the previous case
//   sort( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );
//   reverse( vec_TriggerTauJet.begin(), vec_TriggerTauJet.end() );

  // Single tau trigger
  if ( vec_TriggerTauJet.size() != 0 ) {
    if ( vec_TriggerTauJet[0].pt() >= 150. ) {
      response_tautrig_single = true;
    }
  }
  // Di-tau trigger
  if ( vec_TriggerTauJet.size() >= 1 ) {
    if ( vec_TriggerTauJet[1].pt() >= 80. ) {
      response_tautrig_ditau = true;
    }
  }
  if ( response_tautrig_single || response_tautrig_ditau ) response_tautrig = true;

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
  if ( response_cen ) ++Eff_cen_;
  if ( response_cen_et1 ) ++Eff_cen_et1_;
  if ( response_cen_et2 ) ++Eff_cen_et2_;
  if ( response_cen_et3 ) ++Eff_cen_et3_;
  if ( response_cen_et4 ) ++Eff_cen_et4_;

  if ( response_tau ) ++Eff_tau_;
  if ( response_tau_et1 ) ++Eff_tau_et1_;
  if ( response_tau_et2 ) ++Eff_tau_et2_;
  if ( response_tau_et3 ) ++Eff_tau_et3_;
  if ( response_tau_et4 ) ++Eff_tau_et4_;

  if ( response_for ) ++Eff_for_;
  if ( response_for_et1 ) ++Eff_for_et1_;
  if ( response_for_et2 ) ++Eff_for_et2_;
  if ( response_for_et3 ) ++Eff_for_et3_;
  if ( response_for_et4 ) ++Eff_for_et4_;

  if ( response ) ++Eff_;
  if ( response_Et1 ) ++Eff_et1_;
  if ( response_Et2 ) ++Eff_et2_;
  if ( response_Et3 ) ++Eff_et3_;
  if ( response_Et4 ) ++Eff_et4_;

  if ( response_nofor ) ++Eff_nofor_;
  if ( response_Et1_nofor ) ++Eff_nofor_et1_;
  if ( response_Et2_nofor ) ++Eff_nofor_et2_;
  if ( response_Et3_nofor ) ++Eff_nofor_et3_;
  if ( response_Et4_nofor ) ++Eff_nofor_et4_;

  // MEt+jet
  if ( response_MEtJet_cen ) ++Eff_MEtJet_cen_;
  if ( response_MEtJet_tau ) ++Eff_MEtJet_tau_;
  if ( response_MEtJet_for ) ++Eff_MEtJet_for_;
  if ( response_MEtJet ) ++Eff_MEtJet_;
  if ( response_MEtJet_nofor ) ++Eff_MEtJet_nofor_;

  // Tau trigger
  if ( response_tautrig ) ++Eff_tautrig_;
  if ( response_tautrig_single ) ++Eff_tautrig_single_;
  if ( response_tautrig_ditau ) ++Eff_tautrig_ditau_;

  // Offline
  if ( response && offline ) ++offlineEffMultijet_;
  if ( response_MEtJet && offline ) ++offlineEffMEtJet_;
  if ( response_tautrig && offline ) ++offlineEffTauTrig_;

  // PixelTrigger
  // ------------
  // dz = 0.4 and numtk = 3
  L1PixelTrig<SimplePixelJet> Pixeltrig_3(3, 0.4, 3);
  L1PixelTrig<SimplePixelJet> Pixeltrig_4(4, 0.4, 3);
  L1PixelTrig<SimplePixelJet> Pixeltrig_5(5, 0.4, 3);
  L1PixelTrig<SimplePixelJet> Pixeltrig_6(6, 0.4, 3);

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
  bool Response_cen_et1 = false;
  bool Response_cen_et2 = false;
  bool Response_cen_et3 = false;
  bool Response_cen_et4 = false;
  bool Response_for = false;
  bool Response_for_et1 = false;
  bool Response_for_et2 = false;
  bool Response_for_et3 = false;
  bool Response_for_et4 = false;
  bool Response_tau = false;
  bool Response_tau_et1 = false;
  bool Response_tau_et2 = false;
  bool Response_tau_et3 = false;
  bool Response_tau_et4 = false;

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

  sort( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );
  reverse( vec_TriggerForJet.begin(), vec_TriggerForJet.end() );

  float forJet1 = 0.;
  float forJet2 = 0.;
  float forJet3 = 0.;
  float forJet4 = 0.;
  int forSize = vec_TriggerForJet.size();
  if ( forSize != 0 ) forJet1 = vec_TriggerForJet[0].pt();
  if ( forSize > 1 ) forJet2 = vec_TriggerForJet[1].pt();
  if ( forSize > 2 ) forJet3 = vec_TriggerForJet[2].pt();
  if ( forSize > 3 ) forJet4 = vec_TriggerForJet[3].pt();

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
  bool Response_MEtJet_for = false;
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
            Response_for = false;
            Response_tau = false;
            // Central
            if ( (cenJet1 >= Et1) || (cenJet2 >= Et2) || (cenJet3 >= Et3) || (cenJet4 >= Et4) ) Response_cen = true;
            // Forward
            if ( (forJet1 >= Et1) || (forJet2 >= Et2) || (forJet3 >= Et3) || (forJet4 >= Et4) ) Response_for = true;
            // Tau
            if ( (tauJet1 >= Et1) || (tauJet2 >= Et2) || (tauJet3 >= Et3) || (tauJet4 >= Et4) ) Response_tau = true;
            // Full
            if ( (Response_cen || Response_tau || Response_for) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArray_[PVnum])[multijetCount]);
            if ( (Response_cen || Response_tau || Response_for) && pixelTrigResponse[PVnum] && offline ) ++((offlineEffMultijetPixelArray_[PVnum])[multijetCount]);
            // increase the index, in the inner loop
            ++multijetCount;
          }
        }
      }
    }

    // Multijet single cuts efficiencies
    int et1Count = 0;
    for ( int Et1=160; Et1<260; Et1+=10 ) {
      Response_cen_et1 = false;
      Response_for_et1 = false;
      Response_tau_et1 = false;
      if ( cenJet1 >= Et1 ) Response_cen_et1 = true;
      if ( forJet1 >= Et1 ) Response_for_et1 = true;
      if ( tauJet1 >= Et1 ) Response_tau_et1 = true;
      if ( (Response_cen_et1 || Response_for_et1 || Response_for_et1) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArrayEt1_[PVnum])[et1Count]);
      ++et1Count;
    }
    int et2Count = 0;
    for ( int Et2=110; Et2<210; Et2+=10 ) {
      Response_cen_et2 = false;
      Response_for_et2 = false;
      Response_tau_et2 = false;
      if ( cenJet2 >= Et2 ) Response_cen_et2 = true;
      if ( forJet2 >= Et2 ) Response_for_et2 = true;
      if ( tauJet2 >= Et2 ) Response_tau_et2 = true;
      if ( (Response_cen_et2 || Response_for_et2 || Response_for_et2) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArrayEt2_[PVnum])[et2Count]);
      ++et2Count;
    }
    int et3Count = 0;
    for ( int Et3=55; Et3<105; Et3+=5 ) {
      Response_cen_et3 = false;
      Response_for_et3 = false;
      Response_tau_et3 = false;
      if ( cenJet3 >= Et3 ) Response_cen_et3 = true;
      if ( forJet3 >= Et3 ) Response_for_et3 = true;
      if ( tauJet3 >= Et3 ) Response_tau_et3 = true;
      if ( (Response_cen_et3 || Response_for_et3 || Response_for_et3) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArrayEt3_[PVnum])[et3Count]);
      ++et3Count;
    }
    int et4Count = 0;
    for ( int Et4=35; Et4<85; Et4+=5 ) {
      Response_cen_et4 = false;
      Response_for_et4 = false;
      Response_tau_et4 = false;
      if ( cenJet4 >= Et4 ) Response_cen_et4 = true;
      if ( forJet4 >= Et4 ) Response_for_et4 = true;
      if ( tauJet4 >= Et4 ) Response_tau_et4 = true;
      if ( (Response_cen_et4 || Response_for_et4 || Response_for_et4) && pixelTrigResponse[PVnum] ) ++((EffMultijetPixelArrayEt4_[PVnum])[et4Count]);
      ++et4Count;
    }

    // MEt + Jet
    // ---------
    int mEtJetCount = 0;
    for ( int MEt=35; MEt<85; MEt+=5 ) {
      for ( int Jet=55; Jet<105; Jet+=5 ) {
        Response_MEtJet_cen = false;
        Response_MEtJet_for = false;
        Response_MEtJet_tau = false;
        // Central
        if ( ( L1MEt >= float(MEt) ) && ( cenJet1 >= float(Jet) ) ) Response_MEtJet_cen = true;
        // Forward
        if ( ( L1MEt >= float(MEt) ) && ( forJet1 >= float(Jet) ) ) Response_MEtJet_for = true;
        // Tau
        if ( ( L1MEt >= float(MEt) ) && ( tauJet1 >= float(Jet) ) ) Response_MEtJet_tau = true;
        // Full
        if ( (Response_MEtJet_cen || Response_MEtJet_tau || Response_MEtJet_for) && pixelTrigResponse[PVnum] ) ++((EffMEtJetPixelArray_[PVnum])[mEtJetCount]);
        if ( (Response_MEtJet_cen || Response_MEtJet_tau || Response_MEtJet_for) && pixelTrigResponse[PVnum] && offline ) ++((offlineEffMEtJetPixelArray_[PVnum])[mEtJetCount]);
        ++mEtJetCount;
      }
    }
  } // end loop on PV cuts


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void OfflineAnalyzer::beginJob(const edm::EventSetup&) {
}

// void OfflineAnalyzer::meanHistoSetup( TH1F * meanhisto_ptr, vector<TH1F *> vec_histo ) {
//   meanhisto_ptr->SetBinContent( numdz+1, mean[numdz]->GetMean() );
//   meanhisto_ptr->SetBinError( numdz+1, Vertex_Dz_[numdz]->GetMeanError() );
//   meanhisto_ptr->GetXaxis()->CenterLabels();
//   meanhisto_ptr->GetXaxis()->SetNdivisions(Vertex_Dz_Mean_->GetSize()-2, false);
// }

// ------------ method called once each job just after ending the event loop  ------------
void OfflineAnalyzer::endJob() {

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
  ofstream offlineMultijetPixelEffFile( "offlineMultijetPixelEff.txt" );
  ofstream offlineMEtJetPixelEffFile( "offlineMEtJetPixelEff.txt" );
  for ( int PVnum=0; PVnum<4; ++PVnum ) {
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSize_; ++multijetCount ) {
      float multijetPixelEff = float((EffMultijetPixelArray_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErr = sqrt(float((EffMultijetPixelArray_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixel_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEff );
      EffMultijetPixel_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErr );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff = " << multijetPixelEff << " +- " << multijetPixelEffErr << endl;

      float offlineMultijetPixelEff = float((offlineEffMultijetPixelArray_[PVnum])[multijetCount])/float(eventcounter_);
      float offlineMultijetPixelEffErr = sqrt(float((offlineEffMultijetPixelArray_[PVnum])[multijetCount]))/float(eventcounter_);
      offlineEffMultijetPixel_[PVnum]->SetBinContent( multijetCount+1, offlineMultijetPixelEff );
      offlineEffMultijetPixel_[PVnum]->SetBinError( multijetCount+1, offlineMultijetPixelEffErr );
      offlineMultijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" offline eff = " << offlineMultijetPixelEff << " +- " << offlineMultijetPixelEffErr << endl;

    }
    // Et1
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSizeEt1_; ++multijetCount ) {
      float multijetPixelEffEt1 = float((EffMultijetPixelArrayEt1_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErrEt1 = sqrt(float((EffMultijetPixelArrayEt1_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixelEt1_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEffEt1 );
      EffMultijetPixelEt1_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErrEt1 );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff for Et1 = " << multijetPixelEffEt1 << " +- " << multijetPixelEffErrEt1 << endl;
    }
    // Et2
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSizeEt2_; ++multijetCount ) {
      float multijetPixelEffEt2 = float((EffMultijetPixelArrayEt2_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErrEt2 = sqrt(float((EffMultijetPixelArrayEt2_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixelEt2_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEffEt2 );
      EffMultijetPixelEt2_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErrEt2 );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff for Et2 = " << multijetPixelEffEt2 << " +- " << multijetPixelEffErrEt2 << endl;
    }
    // Et3
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSizeEt3_; ++multijetCount ) {
      float multijetPixelEffEt3 = float((EffMultijetPixelArrayEt3_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErrEt3 = sqrt(float((EffMultijetPixelArrayEt3_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixelEt3_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEffEt3 );
      EffMultijetPixelEt3_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErrEt3 );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff for Et3 = " << multijetPixelEffEt3 << " +- " << multijetPixelEffErrEt3 << endl;
    }
    // Et4
    for ( int multijetCount=0; multijetCount<EffMultijetPixelSizeEt4_; ++multijetCount ) {
      float multijetPixelEffEt4 = float((EffMultijetPixelArrayEt4_[PVnum])[multijetCount])/float(eventcounter_);
      float multijetPixelEffErrEt4 = sqrt(float((EffMultijetPixelArrayEt4_[PVnum])[multijetCount]))/float(eventcounter_);
      EffMultijetPixelEt4_[PVnum]->SetBinContent( multijetCount+1, multijetPixelEffEt4 );
      EffMultijetPixelEt4_[PVnum]->SetBinError( multijetCount+1, multijetPixelEffErrEt4 );
      multijetPixelEffFile << "For PVnum ="<<PVnum<<" and i = "<< multijetCount <<" Eff for Et4 = " << multijetPixelEffEt4 << " +- " << multijetPixelEffErrEt4 << endl;
    }

    for ( int mEtJetCount=0; mEtJetCount<EffMEtJetPixelSize_; ++mEtJetCount ) {
      float mEtJetPixelEff = float((EffMEtJetPixelArray_[PVnum])[mEtJetCount])/float(eventcounter_);
      float mEtJetPixelEffErr = sqrt(float((EffMEtJetPixelArray_[PVnum])[mEtJetCount]))/float(eventcounter_);
      EffMEtJetPixel_[PVnum]->SetBinContent( mEtJetCount+1, mEtJetPixelEff );
      EffMEtJetPixel_[PVnum]->SetBinError( mEtJetCount+1, mEtJetPixelEffErr );
      mEtJetPixelEffFile << "For PVnum = "<<PVnum<<" and i = "<< mEtJetCount <<" Eff = " << mEtJetPixelEff << " +- " << mEtJetPixelEffErr << endl;

      float offlineMEtJetPixelEff = float((offlineEffMEtJetPixelArray_[PVnum])[mEtJetCount])/float(eventcounter_);
      float offlineMEtJetPixelEffErr = sqrt(float((offlineEffMEtJetPixelArray_[PVnum])[mEtJetCount]))/float(eventcounter_);
      offlineEffMEtJetPixel_[PVnum]->SetBinContent( mEtJetCount+1, offlineMEtJetPixelEff );
      offlineEffMEtJetPixel_[PVnum]->SetBinError( mEtJetCount+1, offlineMEtJetPixelEffErr );
      offlineMEtJetPixelEffFile << "For PVnum = "<<PVnum<<" and i = "<< mEtJetCount <<" offline eff = " << offlineMEtJetPixelEff << " +- " << offlineMEtJetPixelEffErr << endl;
    }
  }

  ofstream Effoutputfile( OutputEffFileName.c_str() );

  Effoutputfile << "Multijet trigger" << endl;
  Effoutputfile << "----------------" << endl;
  Effoutputfile << "Eff multijet = " << float(Eff_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et1 = " << float(Eff_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et2 = " << float(Eff_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et3 = " << float(Eff_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et4 = " << float(Eff_et4_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff multijet no-forward = " <<  float(Eff_nofor_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et1 no-forward = " << float(Eff_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et2 no-forward = " << float(Eff_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et3 no-forward = " << float(Eff_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff Et4 no-forward = " << float(Eff_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff central jet = " << float(Eff_cen_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et1 = " << float(Eff_cen_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et2 = " << float(Eff_cen_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et3 = " << float(Eff_cen_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff central jet Et4 = " << float(Eff_cen_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff forward jet = " << float(Eff_for_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et1 = " << float(Eff_for_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et2 = " << float(Eff_for_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et3 = " << float(Eff_for_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff forward jet Et4 = " << float(Eff_for_et4_)/float(eventcounter_) << endl;
  Effoutputfile << endl;
  Effoutputfile << "Eff tau jet = " << float(Eff_tau_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et1 = " << float(Eff_tau_et1_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et2 = " << float(Eff_tau_et2_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et3 = " << float(Eff_tau_et3_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff tau jet Et4 = " << float(Eff_tau_et4_)/float(eventcounter_) << endl;
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
  Effoutputfile << "Eff single tau trigger = " << float(Eff_tautrig_single_)/float(eventcounter_) << endl;
  Effoutputfile << "Eff ditau trigger = " << float(Eff_tautrig_ditau_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "Offline efficiency" << endl;
  Effoutputfile << "------------------" << endl;
  Effoutputfile << "offline efficiency after multijet trigger = " << float(offlineEffMultijet_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after MEt + Jet trigger = " << float(offlineEffMEtJet_)/float(eventcounter_) << endl;
  Effoutputfile << "offline efficiency after Tau trigger = " << float(offlineEffTauTrig_)/float(eventcounter_) << endl;
  Effoutputfile << endl;

  Effoutputfile << "Total events = " << eventcounter_ << endl;

  Effoutputfile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(OfflineAnalyzer);
