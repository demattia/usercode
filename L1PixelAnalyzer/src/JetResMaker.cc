// #define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/JetResMaker.h"

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

L1Trig JetResMaker::L1Trigger;

// Constructors and destructor
// ---------------------------
JetResMaker::JetResMaker(const edm::ParameterSet& iConfig) :
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
  genJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GenJets" ) )
{

  // Now do what ever initialization is needed
  // -----------------------------------------
  eventcounter_=0;
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","L1TrigPixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  OutputFile->cd();

  const float drmin2 = 0.09;                 // matching cone jets/partons
  const float drmin2ext = drmin2*sqrt(2);    // outside ring for background subtraction

  // Create the histograms
  // ---------------------
  char name1[20];
  char name2[20];
  for ( int i=0; i<40; ++i ) {
    sprintf ( name1, "Ptres%d", i );
    Ptres[i] = new TH1F ( name1, name1, 100, -3., 1. );
    sprintf ( name2, "Dr%d", i );
    Dr[i] = new TH1F ( name2, name2, 100, 0., sqrt(drmin2ext) );
  }
  char name3[20];
  for ( int j=0; j<4; ++j ) {
    sprintf ( name3, "Ave_vs_pt%d", j );
    Ave_vs_pt[j] = new TH1F ( name3, name3, 10, 0., 200. );
  }
  char name4[20];
  for ( int j=0; j<4; ++j ) {
    sprintf ( name4, "Res_vs_pt%d", j );
    Res_vs_pt[j] = new TH1F ( name4, name4, 10, 0., 200. );
  }

  // Debug histograms
  // ----------------
  Ipt = new TH1F ( "Ipt", "Pt spectrum of associated partons", 10, 0., 200. );
  Ieta= new TH1F ( "Ieta", "eta spectrum of associated partons", 4, 0., 4.8 ); 
  Totptleft = new TH1F ( "Totptleft", "Total Pt of unassoc partons", 50, 0., 500. );
  Totptass = new TH1F ( "Totptass", "Total Pt of assoc partons", 50, 0., 500. );
  Nleft = new TH1F ( "Nleft", "Number of partons left out", 200, 0., 200. );
  Nass = new TH1F ( "Nass", "Number of partons associated", 10, 0., 10. );
  Njtot = new TH1F ( "Njtot", "Total number of jets", 20, 0., 20. );
  Njgood= new TH1F ( "Njgood", "Number of good jets", 10, 0., 10. );

  // White background for the canvases
  // ---------------------------------
  gROOT->SetStyle("Plain");

}

JetResMaker::~JetResMaker()
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
void JetResMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  // Constants
  // ---------
  const float pi = 3.1415926;
  const float drmin = 0.3;                   // matching cone jets/partons
  const float drmin2 = 0.09;                 // matching cone jets/partons
  const float drmin2ext = drmin2*sqrt(2);    // outside ring for background subtraction
  const float seed_ass_ptmin = 10;           // minimum pt of clusters associated to jets
  const float jet_ass_ptmin = 10;            // minimum pt of jets associated to clusters
  const float l1jet_ass_ptmin = 10;          // minimum pt of l1 jets associated to jets
//  const float max_eta = 2.4;                 // maximum eta of offlineJets and genJets associated

  // Create new collections with jets satisfying the thresholds

  // L1Jets (central and tau)
  BaseJetCollection goodL1Jets;
  BaseJetCollection::const_iterator l1Jet_it = l1eCenJets->begin();
  for( ; l1Jet_it != l1eCenJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  for( l1Jet_it = l1eTauJets->begin(); l1Jet_it != l1eTauJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  for( l1Jet_it = l1eForJets->begin(); l1Jet_it != l1eForJets->end(); ++l1Jet_it ) {
    if ( l1Jet_it->et() > l1jet_ass_ptmin ) {
      goodL1Jets.push_back(*l1Jet_it);
    }
  }
  sort( goodL1Jets.rbegin(), goodL1Jets.rend(), EtSort<BaseJet>() );

  // offlineJets
  OfflineJetCollection goodOfflineJets;
  OfflineJetCollection::const_iterator offlineJet_it = caloJets->begin();
  for( ; offlineJet_it != caloJets->end(); ++offlineJet_it ) {
    //    if ( offlineJet_it->et() > jet_ass_ptmin && fabs(offlineJet_it->eta()) < max_eta ) {
    if ( offlineJet_it->et() > jet_ass_ptmin ) {
      goodOfflineJets.push_back(*offlineJet_it);
    }
  }
  sort( goodOfflineJets.rbegin(), goodOfflineJets.rend(), EtSort<OfflineJet>() );

  // GenJets
  GenJetCollection goodGenJets;
  GenJetCollection::const_iterator genJet_it = genJets->begin();
  for( ; genJet_it != genJets->end(); ++genJet_it ) {
    //    if ( genJet_it->et() > seed_ass_ptmin && fabs(genJet_it->eta()) < max_eta ) {
    if ( genJet_it->et() > seed_ass_ptmin ) {
      goodGenJets.push_back(*genJet_it);
    }
  }
  sort( goodGenJets.rbegin(), goodGenJets.rend(), EtSort<GenJet>() );

  // ATTENTION
  // If the vectors passed for the association are altered, the map is
  // invalidated (contains pointers to the elements of the vectors).
  // ---------

  // Associate L1Jets with offlineJets
  AssociatorEt<BaseJet, OfflineJet> l1OffAssociator( drmin );
  auto_ptr<map<const BaseJet*, const OfflineJet*> > l1JetOfflineJetMap( l1OffAssociator.Associate( goodL1Jets, goodOfflineJets ) );

  // Fill the vector of offlineJets associated to L1Jets
  OfflineJetCollection goodAssocOfflineJets;
  map<const BaseJet*, const OfflineJet*>::const_iterator l1OffMap_it = l1JetOfflineJetMap->begin();
  for( ; l1OffMap_it != l1JetOfflineJetMap->end(); ++l1OffMap_it ) {
    goodAssocOfflineJets.push_back( *(l1OffMap_it->second) );
    cout << "goodAssocOfflineJet et = " << l1OffMap_it->first->et() << endl;
  }

  // Associate this new collection the goodGenJets
  AssociatorEt<GenJet, OfflineJet> genOffAssociator( drmin );
  auto_ptr<map<const GenJet*, const OfflineJet*> > genJetOfflineJetMap( genOffAssociator.Associate( goodGenJets, goodAssocOfflineJets ) );

  // Loop on the map and fill the histograms
  map<const GenJet*, const OfflineJet*>::const_iterator genOffMap_it = genJetOfflineJetMap->begin();
  for( ; genOffMap_it != genJetOfflineJetMap->end(); ++genOffMap_it ) {
    const GenJet * genJetPtr = genOffMap_it->first;
    const OfflineJet * offlineJetPtr = genOffMap_it->second;
    int ipt = (int)(genJetPtr->et()/20.);
    if ( ipt>9 ) ipt=9; 
    int ieta = (int)(fabs(genJetPtr->eta())/1.2);
    if ( ieta>3 ) ieta=3;
    float genJetEt = genJetPtr->et();
    Ptres[ipt*4+ieta]->Fill( (genJetEt - offlineJetPtr->et())/genJetEt );
    float genJetEta = genJetPtr->eta();
    float deta = genJetEta - offlineJetPtr->eta();
    float dphi = pi - fabs(fabs(genJetPtr->phi()-offlineJetPtr->phi())-pi);
    float dr2 = deta*deta+dphi*dphi;
    Dr[ipt*4+ieta]->Fill( sqrt(dr2) );
    Ipt->Fill(genJetEt);
    Ieta->Fill(fabs(genJetEta));
  }


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif

}

// ------------ method called once each job just before starting event loop  ------------
// --------------------------------------------------------------------------------------
void JetResMaker::beginJob(const edm::EventSetup&) {
}


// ------------ method called once each job just after ending the event loop  ------------
// ---------------------------------------------------------------------------------------
void JetResMaker::endJob() {

  // Ok now display some stuff
  // -------------------------
  TF1 * ResFun = new TF1 ( "ResFun", "gaus", -2, 2 );
  ResFun->SetLineWidth(1);
  ResFun->SetLineColor(kBlue);
  TCanvas * JetRes = new TCanvas( "JetRes", "Jet resolution plots", 900, 900 );
  JetRes->SetFillColor(kNone);
  JetRes->SetBorderSize(2);
  JetRes->SetBorderMode(0);
  JetRes->Divide(4,10);
  for ( int i=0; i<10; i++ ) {
    for ( int j=0; j<4; j++ ) {
      JetRes->cd(i*4+j+1);
      Ptres[i*4+j]->SetLabelSize(0.05);
      Ptres[i*4+j]->SetTitleSize(0.05);
      Ptres[i*4+j]->GetXaxis()->SetLabelSize(0.05);
      Ptres[i*4+j]->GetYaxis()->SetLabelSize(0.05);
      Ptres[i*4+j]->SetXTitle("DPt/Pt");
      Ptres[i*4+j]->SetLineColor(kRed);
      Ptres[i*4+j]->Fit("ResFun","V","PE", -1, 1);
      Ave_vs_pt[j]->SetBinContent(i+1,ResFun->GetParameter(1));
      Ave_vs_pt[j]->SetBinError(i+1,ResFun->GetParError(1));
      Res_vs_pt[j]->SetBinContent(i+1,ResFun->GetParameter(2));
      Res_vs_pt[j]->SetBinError(i+1,ResFun->GetParError(2));
    }
  }
  JetRes->Print("JetRes3_dptfits_5c.eps");

  TCanvas * JetResDr = new TCanvas( "JetResDr", "Jet resolution Dr plots", 900, 900 );
  JetResDr->SetFillColor(kNone);
  JetResDr->SetBorderSize(2);
  JetResDr->SetBorderMode(0);
  JetResDr->Divide(4,10);
  for ( int i=0; i<10; i++ ) {
    for ( int j=0; j<4; j++ ) {
      JetResDr->cd(i*4+j+1);
      Dr[i*4+j]->SetLabelSize(0.05);
      Dr[i*4+j]->SetTitleSize(0.05);
      Dr[i*4+j]->GetXaxis()->SetLabelSize(0.05);
      Dr[i*4+j]->GetYaxis()->SetLabelSize(0.05);
      Dr[i*4+j]->SetXTitle("Delta R");
      Dr[i*4+j]->SetLineColor(kRed);
      Dr[i*4+j]->SetFillStyle(1);
      Dr[i*4+j]->SetFillColor(kRed);
      Dr[i*4+j]->Draw();

    }
  }
  JetResDr->Print("JetRes3_dr_5c.eps");
  JetResDr->Update();

  TCanvas * JetRes3a = new TCanvas( "JetRes3a", "Averaeg DPt/Pt vs Pt", 900, 600 );
  JetRes3a->SetFillColor(kNone);
  JetRes3a->SetBorderSize(2);
  JetRes3a->SetBorderMode(0);
  JetRes3a->Divide(2,2);
  for ( int j=0; j<4; j++ ) {
    JetRes3a->cd(j+1);
    Ave_vs_pt[j]->SetLabelSize(0.05);
    Ave_vs_pt[j]->SetTitleSize(0.05);
    Ave_vs_pt[j]->GetXaxis()->SetLabelSize(0.05);
    Ave_vs_pt[j]->GetYaxis()->SetLabelSize(0.05);
    Ave_vs_pt[j]->SetXTitle("Pt (GeV)");
    Ave_vs_pt[j]->SetYTitle("Average DPt/Pt");
    Ave_vs_pt[j]->SetLineColor(kRed);
    Ave_vs_pt[j]->Draw("PE");
  }
  JetRes3a->Print("JetRes3_ave_vs_pt_5c.eps");
  JetRes3a->Update();

  TCanvas * JetRes3 = new TCanvas( "JetRes3", "Resolution vs Pt", 900, 600 );
  JetRes3->SetFillColor(kNone);
  JetRes3->SetBorderSize(2);
  JetRes3->SetBorderMode(0);
  JetRes3->Divide(2,2);
  for ( int j=0; j<4; j++ ) {
    JetRes3->cd(j+1);
    Res_vs_pt[j]->SetLabelSize(0.05);
    Res_vs_pt[j]->SetTitleSize(0.05);
    Res_vs_pt[j]->GetXaxis()->SetLabelSize(0.05);
    Res_vs_pt[j]->GetYaxis()->SetLabelSize(0.05);
    Res_vs_pt[j]->SetXTitle("Pt (GeV)");
    Res_vs_pt[j]->SetYTitle("Resolution %");
    Res_vs_pt[j]->SetLineColor(kRed);
    Res_vs_pt[j]->Draw("PE");
  }
  JetRes3->Print("JetRes3_res_vs_pt_5c.eps");
  JetRes3->Update();

  TCanvas * Debug = new TCanvas( "Debug", "Debug histograms", 600, 600 );
  Debug->SetFillColor(kNone);
  Debug->SetBorderSize(2);
  Debug->SetBorderMode(0);
  Debug->Divide(2,4);
  Debug->cd(1);
  Ipt->SetLabelSize(0.05);
  Ipt->SetTitleSize(0.05);
  Ipt->GetXaxis()->SetLabelSize(0.05);
  Ipt->GetYaxis()->SetLabelSize(0.05);
  Ipt->SetXTitle("Pt (GeV)");
  Ipt->SetLineColor(kRed);
  Ipt->Draw("PE");
  Debug->cd(2);
  Ieta->SetLabelSize(0.05);
  Ieta->SetTitleSize(0.05);
  Ieta->GetXaxis()->SetLabelSize(0.05);
  Ieta->GetYaxis()->SetLabelSize(0.05);
  Ieta->SetXTitle("Ieta");
  Ieta->SetLineColor(kRed);
  Ieta->Draw("PE");
  Debug->cd(3);
  Totptass->SetLabelSize(0.05);
  Totptass->SetTitleSize(0.05);
  Totptass->GetXaxis()->SetLabelSize(0.05);
  Totptass->GetYaxis()->SetLabelSize(0.05);
  Totptass->SetXTitle("Pt (GeV)");
  Totptass->SetLineColor(kRed);
  Totptass->Draw("PE");
  Debug->cd(4);
  Nass->SetLabelSize(0.05);
  Nass->SetTitleSize(0.05);
  Nass->GetXaxis()->SetLabelSize(0.05);
  Nass->GetYaxis()->SetLabelSize(0.05);
  Nass->SetXTitle("partons ass");
  Nass->SetLineColor(kRed);
  Nass->Draw("PE");
  Debug->cd(5);
  Totptleft->SetLabelSize(0.05);
  Totptleft->SetTitleSize(0.05);
  Totptleft->GetXaxis()->SetLabelSize(0.05);
  Totptleft->GetYaxis()->SetLabelSize(0.05);
  Totptleft->SetXTitle("Pt left");
  Totptleft->SetLineColor(kRed);
  Totptleft->Draw("PE");
  Debug->cd(6);
  Nleft->SetLabelSize(0.05);
  Nleft->SetTitleSize(0.05);
  Nleft->GetXaxis()->SetLabelSize(0.05);
  Nleft->GetYaxis()->SetLabelSize(0.05);
  Nleft->SetXTitle("Partons left");
  Nleft->SetLineColor(kRed);
  Nleft->Draw("PE");
  Debug->cd(7);
  Njtot->SetLabelSize(0.05);
  Njtot->SetTitleSize(0.05);
  Njtot->GetXaxis()->SetLabelSize(0.05);
  Njtot->GetYaxis()->SetLabelSize(0.05);
  Njtot->SetXTitle("N all jets");
  Njtot->SetLineColor(kRed);
  Njtot->Draw("PE");
  Debug->cd(8);
  Njgood->SetLabelSize(0.05);
  Njgood->SetTitleSize(0.05);
  Njgood->GetXaxis()->SetLabelSize(0.05);
  Njgood->GetYaxis()->SetLabelSize(0.05);
  Njgood->SetXTitle("N good jets");
  Njgood->SetLineColor(kRed);
  Njgood->Draw();
  Debug->Update();

  // Write out histograms
  // --------------------
  TFile * outfile;
  outfile = new TFile ("JetRes3.hist","RECREATE");
  outfile->cd();
  for ( int i=0; i<40; i++ ) {
    Ptres[i]->Write();
    Dr[i]->Write();
  }
  for ( int j=0; j<4; j++ ) {
    Ave_vs_pt[j]->Write();    
  }
  for ( int j=0; j<4; j++ ) {
    Res_vs_pt[j]->Write();    
  }
  Ipt->Write();
  Ieta->Write();
  Totptass->Write();
  Nass->Write();
  Totptleft->Write();
  Nleft->Write();
  outfile->Close();

}

// Define this as a plug-in
// ------------------------
DEFINE_FWK_MODULE(JetResMaker);
