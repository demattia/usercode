//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/L1TrigPixelAnalyzer.h"

#include "/data/demattia/PJVERTEX_CMSSW/Classes/SimpleJet/SimpleJet.h"

// For the offline jets and corrections
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// For file output
#include <ofstream.h>

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
  JetNumber_IC5_ = new TH1F( "JetNumber_IC5", "Number of IC5 jets", 20, 0, 20 );

  Eff_ = 0;
  eventcounter_ = 0;
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

  JetNumber_IC5_->Fill( vec_TriggerJet.size() );

  // Count the trigger efficiency
  if ( L1Trigger.Response() ) {
    ++Eff_;
  }

  // HiVariables
  // -----------

  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel( CaloJetAlgorithm, caloJets );

  vector<SimpleJet> vec_calojet; 
  // Correct offline jets on the fly
  const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService, iSetup);
  for( reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal ) {
    double scale = corrector->correction( *cal );
    double corPt = scale*cal->pt();
    vec_calojet.push_back( SimpleJet( corPt, cal->eta(), cal->phi() ) );
    uncorr_JetPt_IC5_->Fill( cal->pt() );
    corr_JetPt_IC5_->Fill( corPt );
  }

  HiVar->Fill( vec_calojet );

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
