//#define DEBUG

#include "AnalysisExamples/L1PixelAnalyzer/interface/L1PixelAnalyzer.h"

#include "AnalysisExamples/AnalysisClasses/interface/SimpleJet.h"
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaPhi.h"
#include "AnalysisExamples/AnalysisClasses/interface/L1PixelTrig.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1PixelAnalyzer::L1PixelAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig )
{
  //now do what ever initialization is needed

  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,"RECREATE","L1PixelAnaOutput");
  // The file must be opened first, so that becomes the default position for all the histograms
  OutputFile->cd();

  int binning = 100;
  float first_bin_et = 0.;
  float last_bin_et = 200.;
  float first_bin_eta = -3.;
  float last_bin_eta = 3.;
  float first_bin_phi = -3.15;
  float last_bin_phi = 3.15;

  float first_bin_tke = 0.;
  float last_bin_tke = 100.;
  float first_bin_zip = -5.;
  float last_bin_zip = 5.;
  float first_bin_tip = -0.5;
  float last_bin_tip = 0.5;
  float first_bin_chi2 = 0.;
  float last_bin_chi2 = 50.;
  float first_bin_vxy = -0.5;
  float last_bin_vxy = 0.5;
  float first_bin_vz = -10.;
  float last_bin_vz = 10.;
  float first_bin_xy = -11.;
  float last_bin_xy = 11.;
  float first_bin_z = -30.;
  float last_bin_z = 30.;
  float first_bin_pt = 0.;
  float last_bin_pt = 50.;

  float first_bin_numtk = 0.;
  float last_bin_numtk = 0.;

  float first_bin_eta_ext = -6.;
  float last_bin_eta_ext = 6.;

  L1ExtraCenJetsEt_ = new TH1F( "L1ExtraCenJetsEt", "L1 Central Jets Et", binning, first_bin_et, last_bin_et );
  L1ExtraCenJetsEta_ = new TH1F( "L1ExtraCenJetsEta", "L1 Central Jets Eta", binning, first_bin_eta, last_bin_eta );
  L1ExtraCenJetsPhi_ = new TH1F( "L1ExtraCenJetsPhi", "L1 Central Jets Phi", binning, first_bin_phi, last_bin_phi );

  L1ExtraForJetsEt_ = new TH1F( "L1ExtraForJetsEt", "L1 Forward Jets Et", binning, first_bin_et, last_bin_et );
  L1ExtraForJetsEta_ = new TH1F( "L1ExtraForJetsEta", "L1 Forward Jets Eta", binning, first_bin_eta, last_bin_eta );
  L1ExtraForJetsPhi_ = new TH1F( "L1ExtraForJetsPhi", "L1 Forward Jets Phi", binning, first_bin_phi, last_bin_phi );

  L1ExtraTauJetsEt_ = new TH1F( "L1ExtraTauJetsEt", "L1 Tau Jets Et", binning, first_bin_et, last_bin_et );
  L1ExtraTauJetsEta_ = new TH1F( "L1ExtraTauJetsEta", "L1 Tau Jets Eta", binning, first_bin_eta, last_bin_eta );
  L1ExtraTauJetsPhi_ = new TH1F( "L1ExtraTauJetsPhi", "L1 Tau Jets Phi", binning, first_bin_phi, last_bin_phi );

  L1ExtraIsoEmEt_ = new TH1F( "L1ExtraIsoEmEt", "L1 Iso EGamma Et", binning, first_bin_et, last_bin_et );
  L1ExtraIsoEmEta_ = new TH1F( "L1ExtraIsoEmEta", "L1 Iso EGamma Eta", binning, first_bin_eta, last_bin_eta );
  L1ExtraIsoEmPhi_ = new TH1F( "L1ExtraIsoEmPhi", "L1 Iso EGamma Phi", binning, first_bin_phi, last_bin_phi );

  L1ExtraNonIsoEmEt_ = new TH1F( "L1ExtraNonIsoEmEt", "L1 Non Iso EGamma Et", binning, first_bin_et, last_bin_et );
  L1ExtraNonIsoEmEta_ = new TH1F( "L1ExtraNonIsoEmEta", "L1 Non Iso EGamma Eta", binning, first_bin_eta, last_bin_eta );
  L1ExtraNonIsoEmPhi_ = new TH1F( "L1ExtraNonIsoEmPhi", "L1 Non Iso EGamma Phi", binning, first_bin_phi, last_bin_phi );

  L1ExtraEtMiss_ = new TH1F( "L1ExtraEtMiss", "L1 Missing Et", binning, first_bin_et, last_bin_et );
  L1ExtraEtMissPhi_ = new TH1F( "L1ExtraEtMissPhi", "L1 Missing Et Phi", binning, first_bin_phi, last_bin_phi );
  L1ExtraEtTotal_ = new TH1F( "L1ExtraEtTotal", "L1 Total Et", binning, first_bin_et, last_bin_et );
  L1ExtraEtHad_ = new TH1F( "L1ExtraEtHad", "L1 EtHad", binning, first_bin_et, last_bin_et );

  // L1Pixel
  // -------
  PixelTrack_P_X_ = new TH1F( "PixelTrack_P_X", "PixelTrack momentum along X", binning, first_bin_tke, last_bin_tke );
  PixelTrack_P_Y_ = new TH1F( "PixelTrack_P_Y", "PixelTrack momentum along Y", binning, first_bin_tke, last_bin_tke );
  PixelTrack_P_Z_ = new TH1F( "PixelTrack_P_Z", "PixelTrack momentum along Z", binning, first_bin_tke, last_bin_tke );
  PixelTrack_Pt_ = new TH1F( "PixelTrack_Pt", "PixelTrack Pt", binning, first_bin_pt, last_bin_pt );
  PixelTrack_Eta_ = new TH1F( "PixelTrack_Eta", "PixelTrack Eta", binning, first_bin_eta, last_bin_eta );
  PixelTrack_Phi_ = new TH1F( "PixelTrack_Phi", "PixelTrack Phi", binning, first_bin_phi, last_bin_phi );
  PixelTrack_Charge_ = new TH1F( "PixelTrack_Charge", "PixelTrack Charge", binning, -1.1, 1.1 );
  PixelTrack_Chi2_ = new TH1F( "PixelTrack_Chi2", "PixelTrack Chi2", binning, first_bin_chi2, last_bin_chi2 );
  PixelTrack_Zip_ = new TH1F( "PixelTrack_Zip", "PixelTrack Zip", binning, first_bin_zip, last_bin_zip );
  PixelTrack_Tip_ = new TH1F( "PixelTrack_Tip", "PixelTrack Tip", binning, first_bin_tip, last_bin_tip );
  PixelTrack_Vertex_X_ = new TH1F( "PixelTrack_Vertex_X", "PixelTrack Vertex X", binning, first_bin_vxy, last_bin_vxy );
  PixelTrack_Vertex_Y_ = new TH1F( "PixelTrack_Vertex_Y", "PixelTrack Vertex Y", binning, first_bin_vxy, last_bin_vxy );
  PixelTrack_Vertex_Z_ = new TH1F( "PixelTrack_Vertex_Z", "PixelTrack Vertex Z", binning, first_bin_vz, last_bin_vz );

  PixelHit_X_ = new TH1F( "PixelHit_X", "PixelHit X ", binning, first_bin_xy, last_bin_xy );
  PixelHit_Y_ = new TH1F( "PixelHit_Y", "PixelHit Y", binning, first_bin_xy, last_bin_xy );
  PixelHit_Z_ = new TH1F( "PixelHit_Z", "PixelHit Z", binning, first_bin_z, last_bin_z );
  PixelHit_XY_ = new TH2F( "PixelHit_XY", "PixelHit XY", binning, first_bin_xy, last_bin_xy, binning, first_bin_xy, last_bin_xy );

  // PixelJet
  PixelJet_Pt_ = new TH1F( "PixelJet_Pt", "PixelJet Pt", binning, first_bin_pt, last_bin_pt );
  PixelJet_Eta_ = new TH1F( "PixelJet_Eta", "PixelJet Eta", binning, first_bin_eta, last_bin_eta );
  PixelJet_Phi_ = new TH1F( "PixelJet_Phi", "PixelJet Phi", binning, first_bin_phi, last_bin_phi );
  PixelJet_NumTk_ = new TH1F( "PixelJet_NumTk", "PixelJet number of tracks", binning , first_bin_numtk, last_bin_numtk );
  PixelJet_Vertex_Z_ = new TH1F( "PixelJet_Vertex_Z", "PixelJet Vertex Z", binning, first_bin_vz, last_bin_vz );

  // GenJet
  GenJet_Pt_ = new TH1F( "GenJet_Pt", "GenJet Pt", binning, first_bin_et, last_bin_et );
  GenJet_Eta_ = new TH1F( "GenJet_Eta", "GenJet Eta", binning, first_bin_eta_ext, last_bin_eta_ext );
  GenJet_Phi_ = new TH1F( "GenJet_Phi", "GenJet Phi", binning, first_bin_phi, last_bin_phi );

  //PixelJet vs GenJets Pt
  PJ_PtRes_ = new TProfile( "PJ_PtRes", "PixelJet vs GenJet Pt", binning, first_bin_et, last_bin_et, first_bin_et, last_bin_et );
  PJ_L1J_PtRes_ = new TProfile( "PJ_L1J_PtRes", "PixelJet vs L1Jet Pt", binning, first_bin_et, last_bin_et, first_bin_et, last_bin_et );

  eventcounter = 0;
}


L1PixelAnalyzer::~L1PixelAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  OutputFile->Write();

}

//
// member functions
//

// ------------ method called to for each event  ------------
void
L1PixelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace l1extra;

  // L1 Calo
  edm::Handle < L1EmParticleCollection > l1eIsoEm;
  edm::Handle < L1EmParticleCollection > l1eNonIsoEm;
  edm::Handle < L1JetParticleCollection > l1eCenJets;
  edm::Handle < L1JetParticleCollection > l1eForJets;
  edm::Handle < L1JetParticleCollection > l1eTauJets;
  edm::Handle < L1EtMissParticle > l1eEtMiss;

  // L1 Pixel tracks
  edm::Handle<reco::TrackCollection> pixeltracks;
  edm::Handle<SiPixelRecHitCollection> pixelhits;

  // Pixel jets
  edm::Handle<PixelJetCollection> pixeljets;

  // GenJets
  edm::Handle<reco::GenJetCollection> genJets;

  // should get rid of this try/catch?
  try {
    edm::InputTag L1CJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eCentralJetsSource");
    edm::InputTag L1FJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eForwardJetsSource");
    edm::InputTag L1TauJetLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eTauJetsSource");
    edm::InputTag L1EtMissLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eEtMissSource");
    edm::InputTag L1IsoEGammaLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eIsoEmSource");
    edm::InputTag L1NonIsoEGammaLabel = conf_.getUntrackedParameter<edm::InputTag>("l1eNonIsoEmSource");

    iEvent.getByLabel(L1IsoEGammaLabel, l1eIsoEm);
    iEvent.getByLabel(L1NonIsoEGammaLabel, l1eNonIsoEm);
    iEvent.getByLabel(L1CJetLabel, l1eCenJets);
    iEvent.getByLabel(L1FJetLabel, l1eForJets);
    iEvent.getByLabel(L1TauJetLabel, l1eTauJets);
    iEvent.getByLabel(L1EtMissLabel, l1eEtMiss);

    // L1 Pixel
    // --------
    std::string PixelTracksLabel = conf_.getUntrackedParameter<std::string>("PixelTrackSource");
    std::string PixelHitsLabel = conf_.getUntrackedParameter<std::string>("PixelHitSource");
    iEvent.getByLabel(PixelTracksLabel, pixeltracks);
    iEvent.getByLabel(PixelHitsLabel, pixelhits);

    // PixelJets
    edm::InputTag PixelJetsLabel = conf_.getUntrackedParameter<edm::InputTag>("PixelJetSource");
    iEvent.getByLabel(PixelJetsLabel, pixeljets);

    // GenJets
    std::string GenJetsLabel = conf_.getUntrackedParameter<std::string>("GenJetSource");
    iEvent.getByLabel( GenJetsLabel, genJets );
  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  eventcounter++;
  std::cout << "Event number " << eventcounter << std::endl;

  // Fill the L1Extra histograms

  // Central jets
#ifdef DEBUG
  int cjet_counter = 0;
  std::cout << "number of cjets = " << l1eCenJets->size() << std::endl;
#endif
  for ( L1JetParticleCollection::const_iterator cj = l1eCenJets->begin(); cj != l1eCenJets->end(); ++cj ) {
#ifdef DEBUG
    std::cout << "cj->phi("<<cjet_counter<<") = " << cj->phi() << std::endl;
    std::cout << "cj->eta("<<cjet_counter<<") = " << cj->eta() << std::endl;
    std::cout << "cj->et("<<cjet_counter<<") = " << cj->et() << std::endl;
    ++cjet_counter;
#endif
    L1ExtraCenJetsEt_->Fill( cj->et() );
    L1ExtraCenJetsEta_->Fill( cj->eta() );
    L1ExtraCenJetsPhi_->Fill( cj->phi() );
  }

  // Forward jets
  for (L1JetParticleCollection::const_iterator fj = l1eForJets->begin(); fj != l1eForJets->end(); ++fj) {
    L1ExtraForJetsEt_->Fill( fj->et() );
    L1ExtraForJetsEta_->Fill( fj->eta() );
    L1ExtraForJetsPhi_->Fill( fj->phi() );
  }

  // Tau jets
  for (L1JetParticleCollection::const_iterator tj = l1eTauJets->begin(); tj != l1eTauJets->end(); ++tj) {
    L1ExtraTauJetsEt_->Fill( tj->et() );
    L1ExtraTauJetsEta_->Fill( tj->eta() );
    L1ExtraTauJetsPhi_->Fill( tj->phi() );
  }

  // Isolated EM
  for (L1EmParticleCollection::const_iterator ie = l1eIsoEm->begin(); ie != l1eIsoEm->end(); ++ie) {
    L1ExtraIsoEmEt_->Fill( ie->et() );
    L1ExtraIsoEmEta_->Fill( ie->eta() );
    L1ExtraIsoEmPhi_->Fill( ie->phi() );
  }

  // Non-isolated EM
  for (L1EmParticleCollection::const_iterator ne = l1eNonIsoEm->begin(); ne != l1eNonIsoEm->end(); ++ne) {
    L1ExtraNonIsoEmEt_->Fill( ne->et() );
    L1ExtraNonIsoEmEta_->Fill( ne->eta() );
    L1ExtraNonIsoEmPhi_->Fill( ne->phi() );
  }

  // Energy sums
  L1ExtraEtMiss_->Fill( l1eEtMiss->et());
  L1ExtraEtMissPhi_->Fill( l1eEtMiss->phi() );
  L1ExtraEtTotal_->Fill( l1eEtMiss->etTotal() );
  L1ExtraEtHad_->Fill( l1eEtMiss->etHad() );



  using namespace std;
  // L1 Pixel tracks
  // ---------------
  typedef reco::TrackCollection::const_iterator IT;
  const reco::TrackCollection tracks = *(pixeltracks.product());
#ifdef DEBUG
  std::cout << "Number of pixel-tracks: "<< tracks.size() << " pixel-tracks" << std::endl;
#endif
  for (IT it=tracks.begin(); it!=tracks.end(); ++it) {

    PixelTrack_P_X_->Fill(it->momentum().x());
    PixelTrack_P_Y_->Fill(it->momentum().y());
    PixelTrack_P_Z_->Fill(it->momentum().z());
    PixelTrack_Pt_->Fill(it->pt());
    PixelTrack_Eta_->Fill(it->eta());
    PixelTrack_Phi_->Fill(it->phi());
    PixelTrack_Charge_->Fill(it->charge());
    PixelTrack_Chi2_->Fill(it->chi2());
    PixelTrack_Zip_->Fill(it->dz());
    PixelTrack_Tip_->Fill(it->d0());
    PixelTrack_Vertex_X_->Fill(it->vertex().x());
    PixelTrack_Vertex_Y_->Fill(it->vertex().y());
    PixelTrack_Vertex_Z_->Fill(it->vertex().z());
  }

  // Taken from CMCMSSW/Validation/RecoTrack/plugins/SiPixelTrackingRecHitsValid.cc
  edm::ESHandle<TrackerGeometry> pDD;
  iSetup.get<TrackerDigiGeometryRecord> ().get (pDD);

  for (TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); ++it) 
    {
      DetId detId = ((*it)->geographicalId());

      unsigned int subid = detId.subdetId();
      if ( !((subid==1) || (subid==2)) ) 
	continue; // end subid if

      SiPixelRecHitCollection::range pixelrechitRange = (pixelhits.product())->get(detId);
      SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.first;
      SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.second;
      SiPixelRecHitCollection::const_iterator pixeliter = pixelrechitRangeIteratorBegin;

      //----Loop over rechits for this detId
      for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) 
	{
	  GlobalPoint gp = (*it)->toGlobal(pixeliter->localPosition());

	  detId = (*it)->geographicalId();
	  if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel )
	    {
	      PixelHit_X_->Fill( gp.x() );
	      PixelHit_Y_->Fill( gp.y() );
	      PixelHit_Z_->Fill( gp.z() );
	      PixelHit_XY_->Fill( gp.x(), gp.y() );
	    }
	  else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap )
	    {
	      PXFDetId fdetid(detId);
	      int side  = fdetid.side();
// 	      int disk  = fdetid.disk();
// 	      int blade = fdetid.blade();
	      int panel = fdetid.panel();
// 	      int plaq  = fdetid.module(); // also known as plaquette

	      if ( side==1 ) 
		{
		  if ( panel==1 )
		    {
// 		      mePosxZmPanel1_all_hits->Fill( rechitx );
// 		      mePosyZmPanel1_all_hits->Fill( rechity );
		    }
		  else if ( panel==2 )
		    {
// 		      mePosxZmPanel2_all_hits->Fill( rechitx );
// 		      mePosyZmPanel2_all_hits->Fill( rechity );
		    }
		  else std::cout << "..............................................Wrong panel number !" << std::endl; 
		} // if ( side==1 ) 
	      else if ( side==2 )
		{
		  if ( panel==1 )
		    {
// 		      mePosxZpPanel1_all_hits->Fill( rechitx );
// 		      mePosyZpPanel1_all_hits->Fill( rechity );
		    }
		  else if ( panel==2 )
		    {
// 		      mePosxZpPanel2_all_hits->Fill( rechitx );
// 		      mePosyZpPanel2_all_hits->Fill( rechity );
		    }
		  else std::cout << "..............................................Wrong panel number !" << std::endl; 
		} //else if ( side==2 )
	      else std::cout << ".......................................................Wrong side !" << std::endl;
                
	    } // else if ( detId.subdetId()==PixelSubdetector::PixelEndcap )
	  else std::cout << "We are not in the pixel detector" << (int)detId.subdetId() << endl;

	}
    }
  // ------------------------------------------------ all hits ---------------------------------------------------------------



  // PixelJets
  // ---------
  typedef PixelJetCollection::const_iterator PJ_IT;
  const PixelJetCollection pjs = *(pixeljets.product());
#ifdef DEBUG
  std::cout << "Number of pixel-jets: "<< pjs.size() << " pixel-jets" << std::endl;
#endif
  for (PJ_IT pj_it=pjs.begin(); pj_it!=pjs.end(); ++pj_it) {
    PixelJet_Pt_->Fill(pj_it->pt());
    PixelJet_Eta_->Fill(pj_it->eta());
    PixelJet_Phi_->Fill(pj_it->phi());
    PixelJet_NumTk_->Fill(pj_it->NumTk());
    PixelJet_Vertex_Z_->Fill(pj_it->z());
  }


  if( genJets->size() != 0 ) { 
    for( reco::GenJetCollection::const_iterator genjet_it=genJets->begin(); genjet_it!=genJets->end(); ++genjet_it ) {
      GenJet_Pt_->Fill( genjet_it->pt() );
      GenJet_Eta_->Fill( genjet_it->eta() );
      GenJet_Phi_->Fill( genjet_it->phi() );
    }
  }

  const reco::GenJetCollection genjets = *(genJets.product());

  // Associate PixelJets with themselves
//  Associator<PixelJet, reco::Track> associator( 0.5 );
//  std::auto_ptr<std::map<const PixelJet*, const reco::Track*> > AssocMap( associator.Associate( pjs, tracks ) );

  Associator<PixelJet, reco::GenJet> associator( 0.5 );
  std::auto_ptr<std::map<const PixelJet*, const reco::GenJet*> > AssocMap( associator.Associate( pjs, genjets ) );

#ifdef DEBUG
  std::cout << "AssocMap ptr = " << AssocMap.get() << std::endl;
#endif

  std::map<const PixelJet*, const reco::GenJet*>::const_iterator assoc_it = (*AssocMap).begin();
  for( ; assoc_it != (*AssocMap).end(); ++assoc_it ) {
    PJ_PtRes_->Fill( assoc_it->first->pt(), assoc_it->second->pt() );
#ifdef DEBUG
    std::cout << "PJ_pt = " << assoc_it->first->pt() << " , GenJet_pt = " << assoc_it->second->pt() << std::endl;
#endif
  }

  // Create a new collection starting from L1 Calo Central and Tau jets
  // this is done because tau jets have many fakes and are more numerous than calo central.

  vector<SimpleJet> vec_Jet;

  for ( L1JetParticleCollection::const_iterator cj = l1eCenJets->begin(); cj != l1eCenJets->end(); ++cj ) {
    vec_Jet.push_back( SimpleJet( cj->et(), cj->eta(), cj->phi() ) );
  }
  // Tau jets
  for (L1JetParticleCollection::const_iterator tj = l1eTauJets->begin(); tj != l1eTauJets->end(); ++tj) {
    double tjeta = tj->eta();
// Use a cut in eta, for |eta| > 1.6 the endcaps start (same cut used for pixeljets)
    if ( fabs(tjeta) < 1.6 ) {
      vec_Jet.push_back( SimpleJet( tj->et(), tj->eta(), tj->phi() ) );
    }
  }
  Associator<PixelJet, SimpleJet> PJ_L1J_associator( 0.5 );
  std::auto_ptr<std::map<const PixelJet*, const SimpleJet*> > PJ_L1J_AssocMap( PJ_L1J_associator.Associate( pjs, vec_Jet ) );


  std::map<const PixelJet*, const SimpleJet*>::const_iterator PJ_L1J_assoc_it = (*PJ_L1J_AssocMap).begin();
  for( ; PJ_L1J_assoc_it != (*PJ_L1J_AssocMap).end(); ++PJ_L1J_assoc_it ) {
    PJ_L1J_PtRes_->Fill( PJ_L1J_assoc_it->first->pt(), PJ_L1J_assoc_it->second->pt() );
#ifdef DEBUG
    std::cout << "PJ_pt = " << PJ_L1J_assoc_it->first->pt() << " , SimpleJet_pt = " << PJ_L1J_assoc_it->second->pt() << std::endl;
#endif
  }

  


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1PixelAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1PixelAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelAnalyzer);
