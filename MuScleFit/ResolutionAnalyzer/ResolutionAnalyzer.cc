#ifndef RESOLUTIONANALYZER_CC
#define RESOLUTIONANALYZER_CC

#include "ResolutionAnalyzer.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ResolutionAnalyzer::ResolutionAnalyzer(const edm::ParameterSet& iConfig) :
  theMuonLabel_( iConfig.getParameter<edm::InputTag>( "MuonLabel" ) ),
  theMuonType_( iConfig.getParameter<int>( "MuonType" ) ),
  theRootFileName_( iConfig.getUntrackedParameter<string>("OutputFileName") ),
  debug_( iConfig.getUntrackedParameter<bool>( "Debug" ) )
{
  //now do what ever initialization is needed
  MuScleFitUtils::resfind = iConfig.getParameter<vector<int> >("ResFind");

  outputFile_ = new TFile(theRootFileName_.c_str(), "RECREATE");
  outputFile_->cd();
  fillHistoMap();

  eventCounter_ = 0;
}


ResolutionAnalyzer::~ResolutionAnalyzer()
{
  outputFile_->cd();
  writeHistoMap();
  outputFile_->Close();
}


//
// member functions
//

// ------------ method called to for each event  ------------
void ResolutionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  ++eventCounter_;
  if ( eventCounter_%100 == 0 ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
  }

  Handle<HepMCProduct> evtMC;
  try {
    iEvent.getByLabel ("source", evtMC);
  } catch (...) { 
    cout << "HepMCProduct non existent" << endl;
  }

  Handle<GenParticleCollection> genParticles; 
  try {
    iEvent.getByLabel ("genParticles", genParticles);
    if (debug_>0) cout << "Found genParticles" << endl;
  } catch (...) {
    cout << "GenParticles non existent" << endl;
  }

  Handle<SimTrackContainer> simTracks;
  try {
    iEvent.getByLabel ("g4SimHits",simTracks);
  } catch (...) { 
    cout << "SimTracks not existent" << endl;
  }

  // Take the reco-muons, depending on the type selected in the cfg
  // --------------------------------------------------------------

  vector<reco::LeafCandidate> muons;

  if (theMuonType_==1) { // GlobalMuons
    Handle<reco::MuonCollection> glbMuons;
    iEvent.getByLabel (theMuonLabel_, glbMuons);
    muons = fillMuonCollection(*glbMuons);
  }
  else if (theMuonType_==2) { // StandaloneMuons
    Handle<reco::TrackCollection> saMuons;
    iEvent.getByLabel (theMuonLabel_, saMuons);
    muons = fillMuonCollection(*saMuons);
  }
  else if (theMuonType_==3) { // Tracker tracks
    Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel (theMuonLabel_, tracks);
    muons = fillMuonCollection(*tracks);
  }

  // Find the best reconstructed resonance
  // -------------------------------------

  reco::Particle::LorentzVector recMu1 = reco::Particle::LorentzVector(0,0,0,0);
  reco::Particle::LorentzVector recMu2 = reco::Particle::LorentzVector(0,0,0,0);
  pair<lorentzVector,lorentzVector> recoRes = MuScleFitUtils::findBestRecoRes(muons);
  if (MuScleFitUtils::ResFound) {
    if (debug_>0) {
      cout <<setprecision(9)<< "Pt after findbestrecores: " << (recoRes.first).Pt() << " " 
           << (recoRes.second).Pt() << endl;
      cout << "recMu1 = " << recMu1 << endl;
      cout << "recMu2 = " << recMu2 << endl;
    }
    recMu1 = recoRes.first;
    recMu2 = recoRes.second;
    if (debug_>0) {
      cout << "after recMu1 = " << recMu1 << endl;
      cout << "after recMu2 = " << recMu2 << endl;
      cout << "mu1.pt = " << recMu1.Pt() << endl;
      cout << "mu2.pt = " << recMu2.Pt() << endl;
    }
  }

  // Match the reco muons with the gen and sim tracks
  // ------------------------------------------------

  pair <reco::Particle::LorentzVector, reco::Particle::LorentzVector> genMu = MuScleFitUtils::findGenMuFromRes(evtMC);
  //first is always mu-, second is always mu+
  if(checkDeltaR(genMu.first,recMu1)){
    mapHisto_["hResolPtGenVSMu"]->Fill(genMu.first,(-genMu.first.Pt()+recMu1.Pt())/genMu.first.Pt(),-1);
    mapHisto_["hResolThetaGenVSMu"]->Fill(genMu.first,(-genMu.first.Theta()+recMu1.Theta()),-1);
    mapHisto_["hResolCotgThetaGenVSMu"]->Fill(genMu.first,(-cos(genMu.first.Theta())/sin(genMu.first.Theta())
                                                          +cos(recMu1.Theta())/sin(recMu1.Theta())),-1);
    mapHisto_["hResolEtaGenVSMu"]->Fill(genMu.first,(-genMu.first.Eta()+recMu1.Eta()),-1);
    mapHisto_["hResolPhiGenVSMu"]->Fill(genMu.first,(-genMu.first.Phi()+recMu1.Phi()),-1);
  }
  if(checkDeltaR(genMu.second,recMu2)){
    mapHisto_["hResolPtGenVSMu"]->Fill(genMu.second,(-genMu.second.Pt()+recMu2.Pt())/genMu.second.Pt(),+1);
    mapHisto_["hResolThetaGenVSMu"]->Fill(genMu.second,(-genMu.second.Theta()+recMu2.Theta()),+1);
    mapHisto_["hResolCotgThetaGenVSMu"]->Fill(genMu.second,(-cos(genMu.second.Theta())/sin(genMu.second.Theta())
                                                           +cos(recMu2.Theta())/sin(recMu2.Theta())),+1);
    mapHisto_["hResolEtaGenVSMu"]->Fill(genMu.second,(-genMu.second.Eta()+recMu2.Eta()),+1);
    mapHisto_["hResolPhiGenVSMu"]->Fill(genMu.second,(-genMu.second.Phi()+recMu2.Phi()),+1);
  }
  pair <reco::Particle::LorentzVector, reco::Particle::LorentzVector> simMu = MuScleFitUtils::findSimMuFromRes(evtMC,simTracks);
  //first is always mu-, second is always mu+
  if(checkDeltaR(simMu.first,recMu1)){
    mapHisto_["hResolPtSimVSMu"]->Fill(simMu.first,(-simMu.first.Pt()+recMu1.Pt())/simMu.first.Pt(),-1);
    mapHisto_["hResolThetaSimVSMu"]->Fill(simMu.first,(-simMu.first.Theta()+recMu1.Theta()),-1);
    mapHisto_["hResolCotgThetaSimVSMu"]->Fill(simMu.first,(-cos(simMu.first.Theta())/sin(simMu.first.Theta())
                                                          +cos(recMu1.Theta())/sin(recMu1.Theta())),-1);
    mapHisto_["hResolEtaSimVSMu"]->Fill(simMu.first,(-simMu.first.Eta()+recMu1.Eta()),-1);
    mapHisto_["hResolPhiSimVSMu"]->Fill(simMu.first,(-simMu.first.Phi()+recMu1.Phi()),-1);
  }
  if(checkDeltaR(simMu.second,recMu2)){
    mapHisto_["hResolPtSimVSMu"]->Fill(simMu.second,(-simMu.second.Pt()+recMu2.Pt())/simMu.first.Pt(),+1);
    mapHisto_["hResolThetaSimVSMu"]->Fill(simMu.second,(-simMu.second.Theta()+recMu2.Theta()),+1);
    mapHisto_["hResolCotgThetaSimVSMu"]->Fill(simMu.second,(-cos(simMu.second.Theta())/sin(simMu.second.Theta())
                                                           +cos(recMu2.Theta())/sin(recMu2.Theta())),+1);
    mapHisto_["hResolEtaSimVSMu"]->Fill(simMu.second,(-simMu.second.Eta()+recMu2.Eta()),+1);
    mapHisto_["hResolPhiSimVSMu"]->Fill(simMu.second,(-simMu.second.Phi()+recMu2.Phi()),+1);
  }

  // Fill the mass resolution histograms
  // -----------------------------------
  // check if the recoMuons match the genMuons
  if( MuScleFitUtils::ResFound && checkDeltaR(simMu.first,recMu1) && checkDeltaR(simMu.second,recMu2) ) {
    double resMass = (recMu1+recMu2).mass();
    double genMass = (genMu.first + genMu.second).mass();
    // first is always mu-, second is always mu+
    mapHisto_["hResolMassVSMu"]->Fill(recMu1, -1, recMu2, +1, resMass, genMass);
  }
}

void ResolutionAnalyzer::fillHistoMap() {
  //Reconstructed muon kinematics
  //-----------------------------
  mapHisto_["hRecBestMu"]             = new HParticle         ("hRecBestMu");
  mapHisto_["hRecBestMu_Acc"]         = new HParticle         ("hRecBestMu_Acc"); 
  mapHisto_["hDeltaRecBestMu"]        = new HDelta            ("hDeltaRecBestMu");

  mapHisto_["hRecBestRes"]            = new HParticle         ("hRecBestRes");
  mapHisto_["hRecBestRes_Acc"]        = new HParticle         ("hRecBestRes_Acc"); 
  mapHisto_["hRecBestResVSMu"]        = new HMassVSPart       ("hRecBestResVSMu");

  //Resolution VS muon kinematic
  //----------------------------
  mapHisto_["hResolMassVSMu"]         = new HMassResolutionVSPart ("hResolMassVSMu");
  mapHisto_["hResolPtGenVSMu"]        = new HMassResolutionVSPart ("hResolPtGenVSMu");
  mapHisto_["hResolPtSimVSMu"]        = new HMassResolutionVSPart ("hResolPtSimVSMu");
  mapHisto_["hResolEtaGenVSMu"]       = new HMassResolutionVSPart ("hResolEtaGenVSMu");
  mapHisto_["hResolEtaSimVSMu"]       = new HMassResolutionVSPart ("hResolEtaSimVSMu");
  mapHisto_["hResolThetaGenVSMu"]     = new HMassResolutionVSPart ("hResolThetaGenVSMu");
  mapHisto_["hResolThetaSimVSMu"]     = new HMassResolutionVSPart ("hResolThetaSimVSMu");
  mapHisto_["hResolCotgThetaGenVSMu"] = new HMassResolutionVSPart ("hResolCotgThetaGenVSMu");
  mapHisto_["hResolCotgThetaSimVSMu"] = new HMassResolutionVSPart ("hResolCotgThetaSimVSMu");
  mapHisto_["hResolPhiGenVSMu"]       = new HMassResolutionVSPart ("hResolPhiGenVSMu");
  mapHisto_["hResolPhiSimVSMu"]       = new HMassResolutionVSPart ("hResolPhiSimVSMu");

  //  mapHisto_["hResolRecoMassVSGenMassVSPt"] = new HResolutionVSPart
}

void ResolutionAnalyzer::writeHistoMap() {
  for (map<string, Histograms*>::const_iterator histo=mapHisto_.begin(); 
       histo!=mapHisto_.end(); histo++) {
    (*histo).second->Write();
  }
}

bool ResolutionAnalyzer::checkDeltaR(reco::Particle::LorentzVector& genMu, reco::Particle::LorentzVector& recMu){
  //first is always mu-, second is always mu+
  double deltaR = sqrt(MuScleFitUtils::deltaPhi(recMu.Phi(),genMu.Phi()) * MuScleFitUtils::deltaPhi(recMu.Phi(),genMu.Phi()) +
			 ((recMu.Eta()-genMu.Eta()) * (recMu.Eta()-genMu.Eta())));
  if(deltaR<0.01)
    return true;
  else if (debug_>0)
    cout<<"Reco muon "<<recMu<<" with eta "<<recMu.Eta()<<" and phi "<<recMu.Phi()<<endl
	<<" DOES NOT MATCH with generated muon from resonance: "<<endl
	<<genMu<<" with eta "<<genMu.Eta()<<" and phi "<<genMu.Phi()<<endl;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ResolutionAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ResolutionAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ResolutionAnalyzer);

#endif // RESOLUTIONANALYZER_CC
