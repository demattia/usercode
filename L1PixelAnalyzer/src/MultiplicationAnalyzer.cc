// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/MultiplicationAnalyzer.h"

// Classes to be accessed
// ----------------------
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
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"

// ModJet class used in the multiplication
#include "AnalysisExamples/AnalysisObjects/interface/ModOfflineJet.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

// Constants, enums and typedefs
// -----------------------------

// Static data member definitions
// ------------------------------

L1Trig MultiplicationAnalyzer::L1Trigger;

// Constructors and destructor
// ---------------------------
MultiplicationAnalyzer::MultiplicationAnalyzer(const edm::ParameterSet& iConfig) :
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
  genJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GenJets" ) ),
  minMultiplicationEt_( iConfig.getUntrackedParameter<double>( "MinMultiplicationEt" ) ),
  mEtAlpha_( iConfig.getUntrackedParameter<double>( "MEtAlpha" ) ),
  Ntot_( iConfig.getUntrackedParameter<int>( "Ntot" ) )
{

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;

  // Multiplier, pass the minimum et of jets to be changed and alpha factor for MEt
  multiplier = new Multiplier( 0.3, minMultiplicationEt_, mEtAlpha_ );

  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd();

  // White background for the canvases
  // ---------------------------------
  gROOT->SetStyle("Plain");

  // Control histograms
  NJET = new TH1F ( "NJET", "Njet", 12, -0.5, 11.5 ); 
  JETET = new TH1F ( "JETET", "Jet Et", 100, 0, 250 ); 
  JETET_5j = new TH1F ( "JETET_5j", "Jet Et for NgoodJet >= 5j", 100, 0, 250 ); 
  JETET_5j_s3 = new TH1F ( "JETET_5j_s3", "Jet Et for NgoodJet >= 5j smet>3", 100, 0, 250 ); 
  MET = new TH1F ( "MET", "MEt", 100, 0, 250 ); 
  MET_5j = new TH1F ( "MET_5j", "MEt >= 5j", 100, 0, 250 ); 
  MET_5j_s3 = new TH1F ( "MET_5j_s3", "MEt >= 5j smet>3", 100, 0, 250 ); 
  METX = new TH1F ( "METX", "MEtx", 100, -150, 150 ); 
  METX_5j = new TH1F ( "METX_5j", "MEtx >= 5j", 100, -150, 150 ); 
  METX_5j_s3 = new TH1F ( "METX_5j_s3", "MEtx >= 5j smet>3", 100, -150, 150 ); 
  SET = new TH1F ( "SET", "SEt",100, 0, 2000 ); 
  SET_5j = new TH1F ( "SET_5j", "SEt >= 5j",100, 0, 2000 ); 
  SET_5j_s3 = new TH1F ( "SET_5j_s3", "SEt >= 5j smet>3",100, 0, 2000 ); 
  SMET = new TH1F ( "SMET", "sMEt", 100, 0, 5 ); 
  SMET_5j = new TH1F ( "SMET_5j", "sMEt >= 5j", 100, 0, 5 ); 
  SMET_5j_s3 = new TH1F ( "SMET_5j_s3", "sMEt >= 5j smet>3", 100, 0, 5 ); 
  NJETC = new TH1F ( "NJETC", "Njetc", 12, -0.5, 11.5 ); 
  METC = new TH1F ( "METC", "MEtc", 100, 0, 250 ); 
  METC_5j = new TH1F ( "METC_5j", "MEtc >= 5j", 100, 0, 250 ); 
  METC_5j_s3 = new TH1F ( "METC_5j_s3", "MEtc >= 5j smet>3", 100, 0, 250 ); 
  JETETC = new TH1F ( "JETETC", "changed Jet Et", 100, 0, 250 ); 
  JETETC_5j = new TH1F ( "JETETC_5j", "changed Jet Et for NgoodJet >= 5j", 100, 0, 250 ); 
  JETETC_5j_s3 = new TH1F ( "JETETC_5j_s3", "changed Jet Et for NgoodJet >= 5j smet>3", 100, 0, 250 ); 
  METXC = new TH1F ( "METXC", "MEtxC", 100, -150, 150 ); 
  METXC_5j = new TH1F ( "METXC_5j", "MEtxC >= 5j", 100, -150, 150 ); 
  METXC_5j_s3 = new TH1F ( "METXC_5j_s3", "MEtxC >= 5j smet>3", 100, -150, 150 ); 
  SETC = new TH1F ( "SETC", "SEtc",100, 0, 2000 ); 
  SETC_5j = new TH1F ( "SETC_5j", "SEtc >= 5j",100, 0, 2000 ); 
  SETC_5j_s3 = new TH1F ( "SETC_5j_s3", "SEtc >= 5j smet>3",100, 0, 2000 ); 
  SMETC = new TH1F ( "SMETC", "sMEtc", 100, 0, 5 ); 
  SMETC_5j = new TH1F ( "SMETC_5j", "sMEtc >= 5j", 100, 0, 5 ); 
  SMETC_5j_s3 = new TH1F ( "SMETC_5j_s3", "sMEtc >= 5j smet>3", 100, 0, 5 ); 

  DPtmod = new TH1F ( "DPtmod", "DPtmod", 100, -3, 1 );
  DSET_DMET = new TH2F ( "DSET_DMET", "DSET_DMET", 100, -200, 200, 100, -100, 100 );
  DSET = new TH1F ( "DSET", "DSET", 100, -200, 200 );
  DMET = new TH1F ( "DMET", "DMET", 100, -100, 100 );

  ChangedNum = new TH1F ( "ChangedNum", "Number of changed jets", 12, -0.5, 11.5 ); 
  ChangedRatio = new TH1F ( "ChangedRatio", "Number of jets changed/goodJets", 102, -0.01, 1.01 ); 

  KolmogorovTestMEt = new TH1F ( "KolmogorovTestMEt", "Kolmogorov test for the MEt", 50, -0, 1 );
  KolmogorovTestSumEt = new TH1F ( "KolmogorovTestSumEt", "Kolmogorov test for the SumEt", 50, -0, 1 );
  for( int i=0; i<50; ++i ) {
    stringstream id;
    id << i+1;
    stringstream idFloat;
    idFloat << float(i+1)*0.1;
    KolMEt[i] = new TH1F( TString("MEtForAlpha_" + id.str()), TString("MEt for alpha = " + idFloat.str()), 100, 0 ,250 );
    KolSumEt[i] = new TH1F( TString("SumEtForAlpha_" + id.str()), TString("SumEt for alpha = " + idFloat.str()), 100, 0 ,2000 );
  }
}

MultiplicationAnalyzer::~MultiplicationAnalyzer()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------

  OutputFile->Write();
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void MultiplicationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace anaobj;

  // L1 Calo
  // -------
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  // Should we get rid of this try/catch?
  // ------------------------------------
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
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;

  // SimpleElectrons
  // ---------------
  edm::Handle < SimpleElectronCollection > simpleElectrons;

  // SimpleTaus
  // ----------
  edm::Handle < SimpleTauCollection > simpleTaus;

  // Summary
  // -------
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


  // Offline analysis
  // ----------------

  // MEt
  // ---
  edm::Handle<OfflineMEt> caloMET;
  iEvent.getByLabel( offlineMEtLabel_, caloMET );

  // HiVariables
  // -----------
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Take genParticleCandidates
  // --------------------------
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  // Take the mothers
  // ----------------
  MCParticleCollection::const_iterator MCp = MCpartons->begin();

  // Pixel jets
  // ----------
  edm::Handle<SimplePixelJetCollection> pixelJetsHandle;
  iEvent.getByLabel( simplePixelJetLabel_, pixelJetsHandle );
  const SimplePixelJetCollection pixelJets = *(pixelJetsHandle.product());

  // GenJets
  // -------
  edm::Handle < GenJetCollection > genJets;
  iEvent.getByLabel( genJetLabel_, genJets );

  // -------------------
  // end initializations
  // -------------------

  //  Multiplier multiplier( *l1eCenJets, *l1eTauJets, *l1eForJets, *caloJets );
  vector<ModOfflineJet> * modOfflineJetPtr = multiplier->initialize( *caloJets, *caloMET );
  OfflineMEt * offlineMEt = multiplier->offlineMEt();

  if ( Ntot_ > 1 ) multiplier->fillGen( *genJets );

  int count = 0;
  int numChanged = 0;
  for ( int N=0; N<Ntot_; ++N, ++count ) {
    // Change the et of offlineJets associated to genJets
    if ( N > 0 ) {
      multiplier->setAlpha( mEtAlpha_ );
      multiplier->generate();
      numChanged = multiplier->numChanged();
    }
    // The event has been regenerated, take the new values
    float mEtSig = offlineMEt->mEtSig();

    int numGoodJets = 0;
    vector<ModOfflineJet>::const_iterator good_it = modOfflineJetPtr->begin();
    for( ; good_it != modOfflineJetPtr->end(); ++good_it ) {
      float jetEt = good_it->et();
      if ( N == 0 ) {
        JETET->Fill( jetEt );
        if ( jetEt > 25. && fabs(good_it->eta()) < 3. ) {
          JETET_5j->Fill( jetEt );
          if ( mEtSig > 3 ) JETET_5j_s3->Fill( jetEt );
        }
      }
      // Evalute the number of goodJets
      if ( jetEt > 25. && fabs(good_it->eta()) < 3. ) {
        ++numGoodJets;
      }

      // Only for the second pass
      if ( N > 0 && numChanged != 0 ) {
        // Evaluate differences
        float diff = good_it->origEt() - good_it->et();
        if ( diff != 0 ) DPtmod->Fill(diff/(good_it->origEt()));

        JETETC->Fill( jetEt );
        if ( jetEt > 25. && fabs(good_it->eta()) < 3. ) {
          JETETC_5j->Fill( jetEt );
          if ( mEtSig > 3 ) JETETC_5j_s3->Fill( jetEt );
        }
      }
    }

    // Fill the ratio of changed/goodJets
    ChangedNum->Fill( numChanged );
    ChangedRatio->Fill( float(numChanged)/float(numGoodJets) );

    // First pass, no modifications to the jets and MEt are done
    if ( N == 0 ) {
      NJET->Fill(numGoodJets);
      float MEtEt = offlineMEt->et();
      float MEtPhi = offlineMEt->phi();
      MET->Fill(MEtEt);
      METX->Fill(MEtEt*cos(MEtPhi));
      METX->Fill(MEtEt*sin(MEtPhi));
      SET->Fill(offlineMEt->sumEt());
      SMET->Fill(mEtSig);
      if ( numGoodJets >= 5 ) {
        MET_5j->Fill(MEtEt);
        METX_5j->Fill(MEtEt*cos(MEtPhi));
        METX_5j->Fill(MEtEt*sin(MEtPhi));
        SET_5j->Fill(offlineMEt->sumEt());
        SMET_5j->Fill(mEtSig);
        if ( mEtSig > 3 ) {
        MET_5j_s3->Fill(MEtEt);
        METX_5j_s3->Fill(MEtEt*cos(MEtPhi));
        METX_5j_s3->Fill(MEtEt*sin(MEtPhi));
        SET_5j_s3->Fill(offlineMEt->sumEt());
        SMET_5j_s3->Fill(mEtSig);
        } // if met_smet>3
      } // if Ngj>=5
    }

    // From second pass, new et and MEt
    if ( N > 0 && numChanged != 0 ) {
      NJETC->Fill(numGoodJets);
      float MEtEt = offlineMEt->et();
      float MEtPhi = offlineMEt->phi();
      METC->Fill(MEtEt);
      METXC->Fill(MEtEt*cos(MEtPhi));
      METXC->Fill(MEtEt*sin(MEtPhi));
      SETC->Fill(offlineMEt->sumEt());
      SMETC->Fill(mEtSig);
      if ( numGoodJets >= 5 ) {
        METC_5j->Fill(MEtEt);
        METXC_5j->Fill(MEtEt*cos(MEtPhi));
        METXC_5j->Fill(MEtEt*sin(MEtPhi));
        SETC_5j->Fill(offlineMEt->sumEt());
        SMETC_5j->Fill(mEtSig);
        if ( mEtSig > 3 ) {
        METC_5j_s3->Fill(MEtEt);
        METXC_5j_s3->Fill(MEtEt*cos(MEtPhi));
        METXC_5j_s3->Fill(MEtEt*sin(MEtPhi));
        SETC_5j_s3->Fill(offlineMEt->sumEt());
        SMETC_5j_s3->Fill(mEtSig);
        } // if met_smet>3
      } // if Ngj>=5

      float sumEtDiff = offlineMEt->sumEt() - caloMET->sumEt();
      float mEtDiff = MEtEt - caloMET->et();
      DSET_DMET->Fill( sumEtDiff, mEtDiff );
      DSET->Fill( sumEtDiff );
      DMET->Fill( mEtDiff );

    } // end from second pass

//     vector<ModOfflineJet>::const_iterator off_it = modOfflineJetPtr->begin();
//     int jetCount = 0;
//     for( ; off_it != modOfflineJetPtr->end(); ++off_it, ++jetCount ) {
//       cout << "Generated event = " << count << " Et["<<jetCount<<"] = " << off_it->et() << endl;
//     }

  }


  // Multiple passes to evaluate alpha factor:
  for ( int iAlpha=0; iAlpha<50; ++iAlpha ) {
    double alpha = float(iAlpha+1)*0.02;
    // Change the alpha factor and regenerate the event
    multiplier->setAlpha( alpha );
    multiplier->generate();
    KolMEt[iAlpha]->Fill( offlineMEt->et() );
    KolSumEt[iAlpha]->Fill( offlineMEt->sumEt() );
  }


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif

}

// ------------ method called once each job just before starting event loop  ------------
// --------------------------------------------------------------------------------------
void MultiplicationAnalyzer::beginJob(const edm::EventSetup&) {
}


// ------------ method called once each job just after ending the event loop  ------------
// ---------------------------------------------------------------------------------------
void MultiplicationAnalyzer::endJob() {

  // Perform Kolmogorov test on the MEt and SumEt histograms
  for( int i=0; i<50; ++i ) {
    KolmogorovTestMEt->SetBinContent( i+1, MET->KolmogorovTest( KolMEt[i] ) );
    KolmogorovTestSumEt->SetBinContent( i+1, SET->KolmogorovTest( KolSumEt[i] ) );
  }
  KolmogorovTestMEt->SetLabelSize(0.05);
  KolmogorovTestMEt->SetTitleSize(0.05);
  KolmogorovTestMEt->GetXaxis()->SetLabelSize(0.05);
  KolmogorovTestMEt->GetYaxis()->SetLabelSize(0.05);
  KolmogorovTestMEt->SetXTitle("alpha");

  KolmogorovTestSumEt->SetLabelSize(0.05);
  KolmogorovTestSumEt->SetTitleSize(0.05);
  KolmogorovTestSumEt->GetXaxis()->SetLabelSize(0.05);
  KolmogorovTestSumEt->GetYaxis()->SetLabelSize(0.05);
  KolmogorovTestSumEt->SetXTitle("alpha");

  // Ok now display comparison between original and bootstrapped histograms
  // ----------------------------------------------------------------------
  TCanvas * JetResTest1 = new TCanvas( "JetResTest1", "Jet res test plots", 400, 800 );
  JetResTest1->SetFillColor(kNone);
  JetResTest1->SetBorderSize(2);
  JetResTest1->SetBorderMode(0);
  JetResTest1->Divide(2,3);

  JetResTest1->cd(1);
  NJET->SetLabelSize(0.05);
  NJET->SetTitleSize(0.05);
  NJET->GetXaxis()->SetLabelSize(0.05);
  NJET->GetYaxis()->SetLabelSize(0.05);
  NJET->SetXTitle("Number of jets");
  NJET->SetLineColor(kRed);
  NJET->SetNormFactor(1);
  NJET->Draw();
  NJETC->SetLabelSize(0.05);
  NJETC->SetTitleSize(0.05);
  NJETC->GetXaxis()->SetLabelSize(0.05);
  NJETC->GetYaxis()->SetLabelSize(0.05);
  NJETC->SetXTitle("Number of jets");
  NJETC->SetLineColor(kBlue);
  NJETC->SetNormFactor(1);
  NJETC->Draw("SAME");
  JetResTest1->cd(2);
  MET->SetLabelSize(0.05);
  MET->SetTitleSize(0.05);
  MET->GetXaxis()->SetLabelSize(0.05);
  MET->GetYaxis()->SetLabelSize(0.05);
  MET->SetXTitle("Missing Et (GeV)");
  MET->SetLineColor(kRed);
  MET->SetNormFactor(1);
  MET->Draw();
  METC->SetLabelSize(0.05);
  METC->SetTitleSize(0.05);
  METC->GetXaxis()->SetLabelSize(0.05);
  METC->GetYaxis()->SetLabelSize(0.05);
  METC->SetXTitle("Missing Et (GeV)");
  METC->SetLineColor(kBlue);
  METC->SetNormFactor(1);
  METC->Draw("SAME");
  JetResTest1->cd(3);
  METX->SetLabelSize(0.05);
  METX->SetTitleSize(0.05);
  METX->GetXaxis()->SetLabelSize(0.05);
  METX->GetYaxis()->SetLabelSize(0.05);
  METX->SetXTitle("Missing Et x comp(GeV)");
  METX->SetLineColor(kRed);
  METX->SetNormFactor(1);
  METX->Draw();
  METXC->SetLabelSize(0.05);
  METXC->SetTitleSize(0.05);
  METXC->GetXaxis()->SetLabelSize(0.05);
  METXC->GetYaxis()->SetLabelSize(0.05);
  METXC->SetXTitle("Missing Et x comp (GeV)");
  METXC->SetLineColor(kBlue);
  METXC->SetNormFactor(1);
  METXC->Draw("SAME");
  JetResTest1->cd(4);
  SET->SetLabelSize(0.05);
  SET->SetTitleSize(0.05);
  SET->GetXaxis()->SetLabelSize(0.05);
  SET->GetYaxis()->SetLabelSize(0.05);
  SET->SetXTitle("Sum Et (GeV)");
  SET->SetLineColor(kRed);
  SET->SetNormFactor(1);
  SET->Draw();
  SETC->SetLabelSize(0.05);
  SETC->SetTitleSize(0.05);
  SETC->GetXaxis()->SetLabelSize(0.05);
  SETC->GetYaxis()->SetLabelSize(0.05);
  SETC->SetXTitle("Sum Et (GeV)");
  SETC->SetLineColor(kBlue);
  SETC->SetNormFactor(1);
  SETC->Draw("SAME");
  JetResTest1->cd(5);
  SMET->SetLabelSize(0.05);
  SMET->SetTitleSize(0.05);
  SMET->GetXaxis()->SetLabelSize(0.05);
  SMET->GetYaxis()->SetLabelSize(0.05);
  SMET->SetXTitle("Met significance");
  SMET->SetLineColor(kRed);
  SMET->SetNormFactor(1);
  SMET->Draw();
  SMETC->SetLabelSize(0.05);
  SMETC->SetTitleSize(0.05);
  SMETC->GetXaxis()->SetLabelSize(0.05);
  SMETC->GetYaxis()->SetLabelSize(0.05);
  SMETC->SetXTitle("Met Significance");
  SMETC->SetLineColor(kBlue);
  SMETC->SetNormFactor(1);
  SMETC->Draw("SAME");
  JetResTest1->cd(6);
  JETET->SetLabelSize(0.05);
  JETET->SetTitleSize(0.05);
  JETET->GetXaxis()->SetLabelSize(0.05);
  JETET->GetYaxis()->SetLabelSize(0.05);
  JETET->SetXTitle("Jet Et (GeV)");
  JETET->SetLineColor(kRed);
  JETET->SetNormFactor(1);
  JETET->Draw();
  JETETC->SetLabelSize(0.05);
  JETETC->SetTitleSize(0.05);
  JETETC->GetXaxis()->SetLabelSize(0.05);
  JETETC->GetYaxis()->SetLabelSize(0.05);
  JETETC->SetXTitle("Jet Et (GeV)");
  JETETC->SetLineColor(kBlue);
  JETETC->SetNormFactor(1);
  JETETC->Draw("SAME");
//  JetResTest1->Update();
  OutputFile->cd();
  JetResTest1->Write();

  TCanvas * JetResTest2 = new TCanvas( "JetResTest2", "Jet res test plots - Nj>=4", 400, 800 );
  JetResTest2->SetFillColor(kNone);
  JetResTest2->SetBorderSize(2);
  JetResTest2->SetBorderMode(0);
  JetResTest2->Divide(2,3);

  JetResTest2->cd(2);
  MET_5j->SetLabelSize(0.05);
  MET_5j->SetTitleSize(0.05);
  MET_5j->GetXaxis()->SetLabelSize(0.05);
  MET_5j->GetYaxis()->SetLabelSize(0.05);
  MET_5j->SetXTitle("Missing Et");
  MET_5j->SetLineColor(kRed);
  MET_5j->SetNormFactor(1);
  MET_5j->Draw();
  METC_5j->SetLabelSize(0.05);
  METC_5j->SetTitleSize(0.05);
  METC_5j->GetXaxis()->SetLabelSize(0.05);
  METC_5j->GetYaxis()->SetLabelSize(0.05);
  METC_5j->SetXTitle("Missing Et");
  METC_5j->SetLineColor(kBlue);
  METC_5j->SetNormFactor(1);
  METC_5j->Draw("SAME");
  JetResTest2->cd(3);
  METX_5j->SetLabelSize(0.05);
  METX_5j->SetTitleSize(0.05);
  METX_5j->GetXaxis()->SetLabelSize(0.05);
  METX_5j->GetYaxis()->SetLabelSize(0.05);
  METX_5j->SetXTitle("Missing Et x comp (GeV)");
  METX_5j->SetLineColor(kRed);
  METX_5j->SetNormFactor(1);
  METX_5j->Draw();
  METXC_5j->SetLabelSize(0.05);
  METXC_5j->SetTitleSize(0.05);
  METXC_5j->GetXaxis()->SetLabelSize(0.05);
  METXC_5j->GetYaxis()->SetLabelSize(0.05);
  METXC_5j->SetXTitle("Missing Et x comp (GeV)");
  METXC_5j->SetLineColor(kBlue);
  METXC_5j->SetNormFactor(1);
  METXC_5j->Draw("SAME");
  JetResTest2->cd(4);
  SET_5j->SetLabelSize(0.05);
  SET_5j->SetTitleSize(0.05);
  SET_5j->GetXaxis()->SetLabelSize(0.05);
  SET_5j->GetYaxis()->SetLabelSize(0.05);
  SET_5j->SetXTitle("Sum Et (GeV)");
  SET_5j->SetLineColor(kRed);
  SET_5j->SetNormFactor(1);
  SET_5j->Draw();
  SETC_5j->SetLabelSize(0.05);
  SETC_5j->SetTitleSize(0.05);
  SETC_5j->GetXaxis()->SetLabelSize(0.05);
  SETC_5j->GetYaxis()->SetLabelSize(0.05);
  SETC_5j->SetXTitle("Sum Et (GeV)");
  SETC_5j->SetLineColor(kBlue);
  SETC_5j->SetNormFactor(1);
  SETC_5j->Draw("SAME");
  JetResTest2->cd(5);
  SMET_5j->SetLabelSize(0.05);
  SMET_5j->SetTitleSize(0.05);
  SMET_5j->GetXaxis()->SetLabelSize(0.05);
  SMET_5j->GetYaxis()->SetLabelSize(0.05);
  SMET_5j->SetXTitle("Missing Et significance");
  SMET_5j->SetLineColor(kRed);
  SMET_5j->SetNormFactor(1);
  SMET_5j->Draw();
  SMETC_5j->SetLabelSize(0.05);
  SMETC_5j->SetTitleSize(0.05);
  SMETC_5j->GetXaxis()->SetLabelSize(0.05);
  SMETC_5j->GetYaxis()->SetLabelSize(0.05);
  SMETC_5j->SetXTitle("Missing Et significance");
  SMETC_5j->SetLineColor(kBlue);
  SMETC_5j->SetNormFactor(1);
  SMETC_5j->Draw("SAME");
//  JetResTest2->Update();
  OutputFile->cd();
  JetResTest2->Write();

  TCanvas * JetResTest3 = new TCanvas( "JetResTest3", "Jet res test plots - Nj>=4 s>3", 400, 800 );
  JetResTest3->SetFillColor(kNone);
  JetResTest3->SetBorderSize(2);
  JetResTest3->SetBorderMode(0);
  JetResTest3->Divide(2,3);

  JetResTest3->cd(2);
  MET_5j_s3->SetLabelSize(0.05);
  MET_5j_s3->SetTitleSize(0.05);
  MET_5j_s3->GetXaxis()->SetLabelSize(0.05);
  MET_5j_s3->GetYaxis()->SetLabelSize(0.05);
  MET_5j_s3->SetXTitle("Missing Et");
  MET_5j_s3->SetLineColor(kRed);
  MET_5j_s3->SetNormFactor(1);
  MET_5j_s3->Draw();
  METC_5j_s3->SetLabelSize(0.05);
  METC_5j_s3->SetTitleSize(0.05);
  METC_5j_s3->GetXaxis()->SetLabelSize(0.05);
  METC_5j_s3->GetYaxis()->SetLabelSize(0.05);
  METC_5j_s3->SetXTitle("Missing Et");
  METC_5j_s3->SetLineColor(kBlue);
  METC_5j_s3->SetNormFactor(1);
  METC_5j_s3->Draw("SAME");
  JetResTest3->cd(3);
  METX_5j_s3->SetLabelSize(0.05);
  METX_5j_s3->SetTitleSize(0.05);
  METX_5j_s3->GetXaxis()->SetLabelSize(0.05);
  METX_5j_s3->GetYaxis()->SetLabelSize(0.05);
  METX_5j_s3->SetXTitle("Missing Et x comp (GeV)");
  METX_5j_s3->SetLineColor(kRed);
  METX_5j_s3->SetNormFactor(1);
  METX_5j_s3->Draw();
  METXC_5j_s3->SetLabelSize(0.05);
  METXC_5j_s3->SetTitleSize(0.05);
  METXC_5j_s3->GetXaxis()->SetLabelSize(0.05);
  METXC_5j_s3->GetYaxis()->SetLabelSize(0.05);
  METXC_5j_s3->SetXTitle("Missing Et x comp (GeV)");
  METXC_5j_s3->SetLineColor(kBlue);
  METXC_5j_s3->SetNormFactor(1);
  METXC_5j_s3->Draw("SAME");
  JetResTest3->cd(4);
  SET_5j_s3->SetLabelSize(0.05);
  SET_5j_s3->SetTitleSize(0.05);
  SET_5j_s3->GetXaxis()->SetLabelSize(0.05);
  SET_5j_s3->GetYaxis()->SetLabelSize(0.05);
  SET_5j_s3->SetXTitle("Sum Et (GeV)");
  SET_5j_s3->SetLineColor(kRed);
  SET_5j_s3->SetNormFactor(1);
  SET_5j_s3->Draw();
  SETC_5j_s3->SetLabelSize(0.05);
  SETC_5j_s3->SetTitleSize(0.05);
  SETC_5j_s3->GetXaxis()->SetLabelSize(0.05);
  SETC_5j_s3->GetYaxis()->SetLabelSize(0.05);
  SETC_5j_s3->SetXTitle("Sum Et (GeV)");
  SETC_5j_s3->SetLineColor(kBlue);
  SETC_5j_s3->SetNormFactor(1);
  SETC_5j_s3->Draw("SAME");
  JetResTest3->cd(5);
  SMET_5j_s3->SetLabelSize(0.05);
  SMET_5j_s3->SetTitleSize(0.05);
  SMET_5j_s3->GetXaxis()->SetLabelSize(0.05);
  SMET_5j_s3->GetYaxis()->SetLabelSize(0.05);
  SMET_5j_s3->SetXTitle("Missing Et significance");
  SMET_5j_s3->SetLineColor(kRed);
  SMET_5j_s3->SetNormFactor(1);
  SMET_5j_s3->Draw();
  SMETC_5j_s3->SetLabelSize(0.05);
  SMETC_5j_s3->SetTitleSize(0.05);
  SMETC_5j_s3->GetXaxis()->SetLabelSize(0.05);
  SMETC_5j_s3->GetYaxis()->SetLabelSize(0.05);
  SMETC_5j_s3->SetXTitle("Missing Et significance");
  SMETC_5j_s3->SetLineColor(kBlue);
  SMETC_5j_s3->SetNormFactor(1);
  SMETC_5j_s3->Draw("SAME");
  JetResTest3->Update();
  OutputFile->cd();
  JetResTest3->Write();

  TCanvas * DPt = new TCanvas( "DPt", "DPt before/after correction", 500, 500 );
  DPt->SetFillColor(kNone);
  DPt->SetBorderSize(2);
  DPt->SetBorderMode(0);
  DPt->Divide(2,2);

  DPt->cd(1);
  DPtmod->SetLabelSize(0.05);
  DPtmod->SetTitleSize(0.05);
  DPtmod->GetXaxis()->SetLabelSize(0.05);
  DPtmod->GetYaxis()->SetLabelSize(0.05);
  DPtmod->SetXTitle("DeltaPt/Pt");
  DPtmod->SetLineColor(kBlue);
  DPtmod->SetNormFactor(1);
  DPtmod->Draw();
  DPt->cd(2);
  DMET->Draw();
  DPt->cd(3);
  DSET->Draw();
  DPt->cd(4);
  DSET_DMET->Draw();
//  DPt->Print ("JetResTest6_5.eps");
//  DPt->Update();
  OutputFile->cd();
  DPt->Write();

  // Write out histograms
  // --------------------
//  TFile * outfile = new TFile ("JetResTest6.root","RECREATE");
//  outfile->cd();
//   NJET->Write();
//   SET->Write();
//   SET_5j->Write();
//   SET_5j_s3->Write();
//   MET->Write();
//   MET_5j->Write();
//   MET_5j_s3->Write();
//   SMET->Write();
//   SMET_5j->Write();
//   SMET_5j_s3->Write();
//   NJETC->Write();
//   METC->Write();
//   METC_5j->Write();
//   METC_5j_s3->Write();
//   SETC->Write();
//   SETC_5j->Write();
//   SETC_5j_s3->Write();
//   SMETC->Write();
//   SMETC_5j->Write();
//   SMETC_5j_s3->Write();

//  outfile->Close();
}

// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(MultiplicationAnalyzer);
