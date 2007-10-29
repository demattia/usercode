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
  // Create directory to hold the multiple histograms
  DirVertexDz_ = OutputFile->mkdir("Vertex_dz");

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
  dz_ = 0.04;
  bins_ = 20;

  TString Vertex_Dz_Name("Vertex_Dz_");
  TString Vertex_Dz_Title("Minimum distance between verteces for dz = ");
  TString Vertex_Num_Name("Vertex_Num_");
  TString Vertex_Num_Title("Number of verteces for dz = ");
  TString Prim_Second_Vertex_Dz_Name("Prim_Second_Vertex_Dz_");
  TString Prim_Second_Vertex_Dz_Title("Minimum distance between primary and secondary vertex for dz = ");

  TString PrimVNum_Name( "PrimVNum" );
  TString PrimVNum_Title( "Number of pixeljets in primary vertex for dz = " );
  TString PrimVPt_Name( "PrimVPt" );
  TString PrimVPt_Title( "Pt of primary vertex" );
  TString PrimVEta_Name( "PrimVEta" );
  TString PrimVEta_Title( "Eta of primary vertex" );
  TString PrimVPhi_Name( "PrimVPhi" );
  TString PrimVPhi_Title( "Phi of primary vertex" );

  TString SecVNum_Name( "SecVNum" );
  TString SecVNum_Title( "Number of pixeljets in secondary vertex for dz = " );
  TString SecVPt_Name( "SecVPt" );
  TString SecVPt_Title( "Pt of secondary vertex" );
  TString SecVEta_Name( "SecVEta" );
  TString SecVEta_Title( "Eta of secondary vertex" );
  TString SecVPhi_Name( "SecVPhi" );
  TString SecVPhi_Title( "Phi of secondary vertex" );

  TString AllSecVNum_Name( "AllSecVNum" );
  TString AllSecVNum_Title( "Number of pixeljets in all secondary verteces for dz = " );
  TString AllSecVPt_Name( "AllSecVPt" );
  TString AllSecVPt_Title( "Pt of all secondary verteces" );
  TString AllSecVEta_Name( "AllSecVEta" );
  TString AllSecVEta_Title( "Eta of all secondary verteces" );
  TString AllSecVPhi_Name( "AllSecVPhi" );
  TString AllSecVPhi_Title( "Phi of all secondary verteces" );

  TString Mean("Mean");

  dzmax_ = dz_ + dz_*(bins_);
  Multi_Vertex_Dz_ = new MultiTH1F( "Vertex_Dz", "Minimum distance between verteces", 100, 0., 10., bins_, dz_, dzmax_, OutputFile );

  Vertex_Dz_Mean_ = new TH1F( Vertex_Dz_Name + Mean, Vertex_Dz_Title + Mean, bins_, dz_, dzmax_ );
  Vertex_Num_Mean_ = new TH1F( Vertex_Num_Name + Mean, Vertex_Num_Title + Mean, bins_, dz_, dzmax_ );
  Prim_Second_Vertex_Dz_Mean_ = new TH1F( Prim_Second_Vertex_Dz_Name + Mean, Prim_Second_Vertex_Dz_Title + Mean, bins_, dz_, dzmax_ );

  PrimVNum_Mean_ = new TH1F( PrimVNum_Name + Mean, PrimVNum_Title + Mean, bins_, dz_, dzmax_ );
  PrimVPt_Mean_ = new TH1F( PrimVPt_Name + Mean, PrimVPt_Title + Mean, bins_, dz_, dzmax_ );
  PrimVEta_Mean_ = new TH1F( PrimVEta_Name + Mean, PrimVEta_Title + Mean, bins_, dz_, dzmax_ );
  PrimVPhi_Mean_ = new TH1F( PrimVPhi_Name + Mean, PrimVPhi_Title + Mean, bins_, dz_, dzmax_ );

  SecVNum_Mean_ = new TH1F( SecVNum_Name + Mean, SecVNum_Title + Mean, bins_, dz_, dzmax_ );
  SecVPt_Mean_ = new TH1F( SecVPt_Name + Mean, SecVPt_Title + Mean, bins_, dz_, dzmax_ );
  SecVEta_Mean_ = new TH1F( SecVEta_Name + Mean, SecVEta_Title + Mean, bins_, dz_, dzmax_ );
  SecVPhi_Mean_ = new TH1F( SecVPhi_Name + Mean, SecVPhi_Title + Mean, bins_, dz_, dzmax_ );

  AllSecVNum_Mean_ = new TH1F( AllSecVNum_Name + Mean, AllSecVNum_Title + Mean, bins_, dz_, dzmax_ );
  AllSecVPt_Mean_ = new TH1F( AllSecVPt_Name + Mean, AllSecVPt_Title + Mean, bins_, dz_, dzmax_ );
  AllSecVEta_Mean_ = new TH1F( AllSecVEta_Name + Mean, AllSecVEta_Title + Mean, bins_, dz_, dzmax_ );
  AllSecVPhi_Mean_ = new TH1F( AllSecVPhi_Name + Mean, AllSecVPhi_Title + Mean, bins_, dz_, dzmax_ );

  // Generate the multiple histograms
  ostringstream snum;
  // Put these histograms in the subdir
  DirVertexDz_->cd();
  for ( int num = 0; num < bins_; ++num ) {
    double dz_for_name = dz_ + dz_*(num);
    snum << dz_for_name;

//     string teststring("test");
//     string temp = (teststring + snum.str()).c_str();

    Vertex_Dz_.push_back( new TH1F( Vertex_Dz_Name + snum.str(), Vertex_Dz_Title + snum.str(), 100, 0, 10 ) );
    Vertex_Num_.push_back( new TH1F( Vertex_Num_Name + snum.str(), Vertex_Num_Title + snum.str(), 30, 0, 30 ) );
    Prim_Second_Vertex_Dz_.push_back( new TH1F( Prim_Second_Vertex_Dz_Name + snum.str(), Prim_Second_Vertex_Dz_Title + snum.str(), 100, 0, 10 ) );

    // Use big values for pt or the means will be biased (underflow and overflow are not used)

    PrimVNum_.push_back( new TH1F( PrimVNum_Name + snum.str(), PrimVNum_Title + snum.str(), 30, 0, 30 ) );
    PrimVPt_.push_back( new TH1F( PrimVPt_Name + snum.str(), PrimVPt_Title + snum.str(), 100, 0, 1000 ) );
    PrimVEta_.push_back( new TH1F( PrimVEta_Name + snum.str(), PrimVEta_Title + snum.str(), 100, -3, 3 ) );
    PrimVPhi_.push_back( new TH1F( PrimVPhi_Name + snum.str(), PrimVPhi_Title + snum.str(), 100, -3.15, 3.15 ) );

    SecVNum_.push_back( new TH1F( SecVNum_Name + snum.str(), SecVNum_Title + snum.str(), 30, 0, 30 ) );
    SecVPt_.push_back( new TH1F( SecVPt_Name + snum.str(), SecVPt_Title + snum.str(), 100, 0, 1000 ) );
    SecVEta_.push_back( new TH1F( SecVEta_Name + snum.str(), SecVEta_Title + snum.str(), 100, -3, 3 ) );
    SecVPhi_.push_back( new TH1F( SecVPhi_Name + snum.str(), SecVPhi_Title + snum.str(), 100, -3.15, 3.15 ) );

    AllSecVNum_.push_back( new TH1F( AllSecVNum_Name + snum.str(), AllSecVNum_Title + snum.str(), 30, 0, 30 ) );
    AllSecVPt_.push_back( new TH1F( AllSecVPt_Name + snum.str(), AllSecVPt_Title + snum.str(), 100, 0, 1000 ) );
    AllSecVEta_.push_back( new TH1F( AllSecVEta_Name + snum.str(), AllSecVEta_Title + snum.str(), 100, -3, 3 ) );
    AllSecVPhi_.push_back( new TH1F( AllSecVPhi_Name + snum.str(), AllSecVPhi_Title + snum.str(), 100, -3.15, 3.15 ) );

    // Empty the ostringstream
    snum.str("");
  }
  // ---------------------------

  // Go back to the main file
  OutputFile->cd();

  DPhimin_ = new TH1F( "DPhimin", "Minimum distance in (R,Phi) between MEt and closest jet", 100, 0, 3.15 );

  Eff_ = 0;
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


  // should get rid of this try/catch?
  try {
    edm::InputTag L1CJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eCentralJetsSource");
    edm::InputTag L1FJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eForwardJetsSource");
    edm::InputTag L1TauJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eTauJetsSource");

    iEvent.getByLabel(L1CJetLabel, l1eCenJets);
    iEvent.getByLabel(L1FJetLabel, l1eForJets);
    iEvent.getByLabel(L1TauJetLabel, l1eTauJets);

  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // Level 1 trigger
  // ---------------

  // All the jets together fot the L1Trigger
  vector<SimpleJet> vec_TriggerJet;
  for ( L1JetParticleCollection::const_iterator tcj = l1eCenJets->begin(); tcj != l1eCenJets->end(); ++tcj ) {
    vec_TriggerJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
  }
  for ( L1JetParticleCollection::const_iterator tfj = l1eForJets->begin(); tfj != l1eForJets->end(); ++tfj ) {
    vec_TriggerJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
  }
  // Tau jets
  for ( L1JetParticleCollection::const_iterator ttj = l1eTauJets->begin(); ttj != l1eTauJets->end(); ++ttj ) {
    vec_TriggerJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
  }

  L1Trigger.Fill( vec_TriggerJet );

#ifdef DEBUG
  std::cout << "Number of L1Jets = " << vec_TriggerJet.size() << std::endl;
  std::cout << "L1Trigger response = " << L1Trigger.Response() << std::endl;
#endif

  // Count the trigger efficiency
  if ( L1Trigger.Response() ) {
    ++Eff_;
  }


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

  // Pixel jets
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
    Vertex_Num_[numdz]->Fill( numVertex );

#ifdef DEBUG
    std::cout << "Loop on primary vertex" << std::endl;
#endif
    if ( numVertex > 0 ) {
      //Primary vertex characteristics
      Vertex<PixelJet> prim_vertex( vec_vertex_aptr->back() );

      PrimVNum_[numdz]->Fill( prim_vertex.number() );
      PrimVPt_[numdz]->Fill( prim_vertex.pt() );
      PrimVEta_[numdz]->Fill( prim_vertex.eta() );
      PrimVPhi_[numdz]->Fill( prim_vertex.phi() );

      // Loop on the remaining verteces and fill their characteristics
      if ( numVertex > 1 ) {
        // End returns just just after the last element, to take the one before the last go two times back
        vector<Vertex<PixelJet> >::const_iterator svec_it = (vec_vertex_aptr->end()-2);
        SecVNum_[numdz]->Fill( svec_it->number() );
        SecVPt_[numdz]->Fill( svec_it->pt() );
        SecVEta_[numdz]->Fill( svec_it->eta() );
        SecVPhi_[numdz]->Fill( svec_it->phi() );

        vector<Vertex<PixelJet> >::const_iterator all_svec_it = (vec_vertex_aptr->begin());
        for ( ; all_svec_it != vec_vertex_aptr->end()-1; ++all_svec_it ) {
          AllSecVNum_[numdz]->Fill( all_svec_it->number() );
          AllSecVPt_[numdz]->Fill( all_svec_it->pt() );
          AllSecVEta_[numdz]->Fill( all_svec_it->eta() );
          AllSecVPhi_[numdz]->Fill( all_svec_it->phi() );
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
          Vertex_Dz_[numdz]->Fill( fabs(svec_z_it->z() - (svec_z_it-1)->z()) );      
          Multi_Vertex_Dz_->Fill( fabs(svec_z_it->z() - (svec_z_it-1)->z()), numdz );
        }
      }

      // Evaluate the distance between primary and secondary vertex
      if ( numVertex > 1 ) {
        // Take the last and the one before
        vector<Vertex<PixelJet> >::const_iterator primV_secondV_z_it = (vec_vertex_aptr->end()-1);
        Prim_Second_Vertex_Dz_[numdz]->Fill( fabs(primV_secondV_z_it->z() - (primV_secondV_z_it-1)->z()) );      
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

  Vertex_Dz_Mean_->GetXaxis()->CenterLabels();
  Vertex_Dz_Mean_->GetXaxis()->SetNdivisions(Vertex_Dz_Mean_->GetSize()-2, false);
  Vertex_Num_Mean_->GetXaxis()->CenterLabels();
  Vertex_Num_Mean_->GetXaxis()->SetNdivisions(Vertex_Num_Mean_->GetSize()-2, false);
  Prim_Second_Vertex_Dz_Mean_->GetXaxis()->CenterLabels();
  Prim_Second_Vertex_Dz_Mean_->GetXaxis()->SetNdivisions(Prim_Second_Vertex_Dz_Mean_->GetSize()-2, false);

  PrimVNum_Mean_->GetXaxis()->CenterLabels();
  PrimVPt_Mean_->GetXaxis()->CenterLabels();
  PrimVNum_Mean_->GetXaxis()->SetNdivisions(PrimVPt_Mean_->GetSize()-2, false);
  PrimVEta_Mean_->GetXaxis()->CenterLabels();
  PrimVEta_Mean_->GetXaxis()->SetNdivisions(PrimVEta_Mean_->GetSize()-2, false);
  PrimVPhi_Mean_->GetXaxis()->CenterLabels();
  PrimVPhi_Mean_->GetXaxis()->SetNdivisions(PrimVPhi_Mean_->GetSize()-2, false);

  SecVNum_Mean_->GetXaxis()->CenterLabels();
  SecVNum_Mean_->GetXaxis()->SetNdivisions(SecVNum_Mean_->GetSize()-2, false);
  SecVPt_Mean_->GetXaxis()->CenterLabels();
  SecVPt_Mean_->GetXaxis()->SetNdivisions(SecVNum_Mean_->GetSize()-2, false);
  SecVEta_Mean_->GetXaxis()->CenterLabels();
  SecVEta_Mean_->GetXaxis()->SetNdivisions(SecVNum_Mean_->GetSize()-2, false);
  SecVPhi_Mean_->GetXaxis()->CenterLabels();
  SecVPhi_Mean_->GetXaxis()->SetNdivisions(SecVNum_Mean_->GetSize()-2, false);

  AllSecVNum_Mean_->GetXaxis()->CenterLabels();
  AllSecVNum_Mean_->GetXaxis()->SetNdivisions(AllSecVNum_Mean_->GetSize()-2, false);
  AllSecVPt_Mean_->GetXaxis()->CenterLabels();
  AllSecVPt_Mean_->GetXaxis()->SetNdivisions(AllSecVNum_Mean_->GetSize()-2, false);
  AllSecVEta_Mean_->GetXaxis()->CenterLabels();
  AllSecVEta_Mean_->GetXaxis()->SetNdivisions(AllSecVNum_Mean_->GetSize()-2, false);
  AllSecVPhi_Mean_->GetXaxis()->CenterLabels();
  AllSecVPhi_Mean_->GetXaxis()->SetNdivisions(AllSecVNum_Mean_->GetSize()-2, false);

  // Take the means of the Vertex_dz histograms
  for( int numdz = 0; numdz < bins_; ++numdz ) {

    Vertex_Dz_Mean_->SetBinContent( numdz+1, Vertex_Dz_[numdz]->GetMean() );
    Vertex_Dz_Mean_->SetBinError( numdz+1, Vertex_Dz_[numdz]->GetMeanError() );
    Vertex_Num_Mean_->SetBinContent( numdz+1, Vertex_Num_[numdz]->GetMean() );;
    Vertex_Num_Mean_->SetBinError( numdz+1, Vertex_Num_[numdz]->GetMeanError() );;
    Prim_Second_Vertex_Dz_Mean_->SetBinContent( numdz+1, Prim_Second_Vertex_Dz_[numdz]->GetMean() );
    Prim_Second_Vertex_Dz_Mean_->SetBinError( numdz+1, Prim_Second_Vertex_Dz_[numdz]->GetMeanError() );

    PrimVNum_Mean_->SetBinContent( numdz+1, PrimVNum_[numdz]->GetMean() );
    PrimVNum_Mean_->SetBinError( numdz+1, PrimVNum_[numdz]->GetMeanError() );
    PrimVPt_Mean_->SetBinContent( numdz+1, PrimVPt_[numdz]->GetMean() );
    PrimVPt_Mean_->SetBinError( numdz+1, PrimVPt_[numdz]->GetMeanError() );
    PrimVEta_Mean_->SetBinContent( numdz+1, PrimVEta_[numdz]->GetMean() );
    PrimVEta_Mean_->SetBinError( numdz+1, PrimVEta_[numdz]->GetMeanError() );
    PrimVPhi_Mean_->SetBinContent( numdz+1, PrimVPhi_[numdz]->GetMean() );
    PrimVPhi_Mean_->SetBinError( numdz+1, PrimVPhi_[numdz]->GetMeanError() );

    SecVNum_Mean_->SetBinContent( numdz+1, SecVNum_[numdz]->GetMean() );
    SecVNum_Mean_->SetBinError( numdz+1, SecVNum_[numdz]->GetMeanError() );
    SecVPt_Mean_->SetBinContent( numdz+1, SecVPt_[numdz]->GetMean() );
    SecVPt_Mean_->SetBinError( numdz+1, SecVPt_[numdz]->GetMeanError() );
    SecVEta_Mean_->SetBinContent( numdz+1, SecVEta_[numdz]->GetMean() );
    SecVEta_Mean_->SetBinError( numdz+1, SecVEta_[numdz]->GetMeanError() );
    SecVPhi_Mean_->SetBinContent( numdz+1, SecVPhi_[numdz]->GetMean() );
    SecVPhi_Mean_->SetBinError( numdz+1, SecVPhi_[numdz]->GetMeanError() );

    AllSecVNum_Mean_->SetBinContent( numdz+1, AllSecVNum_[numdz]->GetMean() );
    AllSecVNum_Mean_->SetBinError( numdz+1, AllSecVNum_[numdz]->GetMeanError() );
    AllSecVPt_Mean_->SetBinContent( numdz+1, AllSecVPt_[numdz]->GetMean() );
    AllSecVPt_Mean_->SetBinError( numdz+1, AllSecVPt_[numdz]->GetMeanError() );
    AllSecVEta_Mean_->SetBinContent( numdz+1, AllSecVEta_[numdz]->GetMean() );
    AllSecVEta_Mean_->SetBinError( numdz+1, AllSecVEta_[numdz]->GetMeanError() );
    AllSecVPhi_Mean_->SetBinContent( numdz+1, AllSecVPhi_[numdz]->GetMean() );
    AllSecVPhi_Mean_->SetBinError( numdz+1, AllSecVPhi_[numdz]->GetMeanError() );
  }

  ofstream Effoutputfile( OutputEffFileName.c_str() );

  Effoutputfile << float(Eff_)/float(eventcounter_);
  Effoutputfile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrigPixelAnalyzer);
