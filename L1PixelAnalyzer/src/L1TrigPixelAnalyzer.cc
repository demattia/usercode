//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/L1TrigPixelAnalyzer.h"

#include "/data/demattia/PJVERTEX_CMSSW/Classes/SimpleJet/SimpleJet.h"
#include "/data/demattia/PJVERTEX_CMSSW/Classes/Associator/Associator.h"
#include "/data/demattia/PJVERTEX_CMSSW/Classes/DeltaPhi/DeltaPhi.h"
#include "/data/demattia/PJVERTEX_CMSSW/Classes/L1PixelTrig/L1PixelTrig.h"

// For the offline jets and corrections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// For file output
#include <fstream>

#include <cmath>

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
  MEt_CorrIC5_Pt_ = new TH1F( "MEt_CorrIC5_Pt", "MEt corrected with IC5 jets", 100, 0, 100 );
  MEt_CorrIC5_Phi_ = new TH1F( "MEt_CorrIC5_Phi", "MEt Phi corrected with IC5 jets", 100, -3.15, 3.15 );
  MEt_CorrIC5_SumEt_ = new TH1F( "MEt_CorrIC5_SumEt", "SumEt corrected with IC5 jets", 100, 0, 2000 );
  MEt_CorrIC5_mEtSig_ = new TH1F( "MEt_CorrIC5_mEtSig", "MEt significance corrected with IC5 jets", 100, 0, 10 );
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
  std::cout << "First of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
  const Candidate* mother_1 = &*Mothers;
  ++Mothers;
  std::cout << "Second of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
  const Candidate* mother_2 = &*Mothers;
  ++Mothers;
  // Take the other mother only if it is a Higgs
  const Candidate* mother_3 = 0;
  if ( Mothers->pdgId() == 25 ) {
    std::cout << "Third of the mothers has pdgId = " << Mothers->pdgId() << std::endl;
    mother_3 = &*Mothers;
  }

  std::vector<const Candidate*> vec_Partons;
  int counter = 0;
  ofstream Partons( "Partons.txt", ios::app );
  Partons << std::endl;
  for( CandidateCollection::const_iterator MCp = MCpartons->begin()+8; MCp != MCpartons->end(); ++MCp, ++counter ) {
    const Candidate* temp_mother = MCp->mother();
    if ( temp_mother != 0 ) {

      // Store the status = 3 partons
      if ( MCp->status() == 3 ) {
        vec_Partons.push_back( &*MCp );
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

  // Pixel trigger requiring at least 5 pixel jets coming from the primary vertex
  // (constructed from pixel jets, and taken as the vertex with the highest ptsum)
  L1PixelTrig<PixelJet> PJtrig(2);

  PJtrig.Fill( pixeljets );

  std::cout << "Pixel trigger response = " << PJtrig.Response() << std::endl;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void L1TrigPixelAnalyzer::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void L1TrigPixelAnalyzer::endJob() {

  ofstream Effoutputfile( OutputEffFileName.c_str() );

  Effoutputfile << float(Eff_)/float(eventcounter_);
  Effoutputfile.close();

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrigPixelAnalyzer);
