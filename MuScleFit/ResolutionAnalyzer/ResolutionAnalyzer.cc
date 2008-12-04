#ifndef RESOLUTIONANALYZER_CC
#define RESOLUTIONANALYZER_CC

#include "ResolutionAnalyzer.h"
#include "Functions.h"

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

  // Initial parameters values
  // -------------------------
  int resolFitType = iConfig.getParameter<int>("ResolFitType");
  MuScleFitUtils::ResolFitType = resolFitType;
  MuScleFitUtils::resolutionFunction = resolutionFunctionArray[resolFitType];
  MuScleFitUtils::resolutionFunctionForVec = resolutionFunctionArrayForVec[resolFitType];

  MuScleFitUtils::parResol = iConfig.getParameter<vector<double> >("parResol");

  MuScleFitUtils::resfind = iConfig.getParameter<vector<int> >("ResFind");

  outputFile_ = new TFile(theRootFileName_.c_str(), "RECREATE");
  outputFile_->cd();
  fillHistoMap();

  eventCounter_ = 0;
  resonance_ = iConfig.getUntrackedParameter<bool>( "Resonance" );
}


ResolutionAnalyzer::~ResolutionAnalyzer()
{
  outputFile_->cd();
  writeHistoMap();
  outputFile_->Close();
  cout << "Total analyzed events = " << eventCounter_ << endl;
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
    cout << "SimTracks not existent, not using them" << endl;
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

  if ( resonance_ ) {

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

    // Histograms with genParticles characteristics
    // --------------------------------------------

    //first is always mu-, second is always mu+
    pair <reco::Particle::LorentzVector, reco::Particle::LorentzVector> genMu = MuScleFitUtils::findGenMuFromRes(evtMC);

    reco::Particle::LorentzVector genMother( genMu.first + genMu.second );

    mapHisto_["GenMother"]->Fill( genMother );
    mapHisto_["DeltaGenMotherMuons"]->Fill( genMu.first, genMu.second );
    mapHisto_["GenMotherMuons"]->Fill( genMu.first );
    mapHisto_["GenMotherMuons"]->Fill( genMu.second );

    // Match the reco muons with the gen and sim tracks
    // ------------------------------------------------
    if(checkDeltaR(genMu.first,recMu1)){
      mapHisto_["PtResolutionGenVSMu"]->Fill(recMu1,(-genMu.first.Pt()+recMu1.Pt())/genMu.first.Pt(),-1);
      mapHisto_["ThetaResolutionGenVSMu"]->Fill(recMu1,(-genMu.first.Theta()+recMu1.Theta()),-1);
      mapHisto_["CotgThetaResolutionGenVSMu"]->Fill(recMu1,(-cos(genMu.first.Theta())/sin(genMu.first.Theta())
                                                            +cos(recMu1.Theta())/sin(recMu1.Theta())),-1);
      mapHisto_["EtaResolutionGenVSMu"]->Fill(recMu1,(-genMu.first.Eta()+recMu1.Eta()),-1);
      mapHisto_["PhiResolutionGenVSMu"]->Fill(recMu1,(-genMu.first.Phi()+recMu1.Phi()),-1);
    }
    if(checkDeltaR(genMu.second,recMu2)){
      mapHisto_["PtResolutionGenVSMu"]->Fill(recMu2,(-genMu.second.Pt()+recMu2.Pt())/genMu.second.Pt(),+1);
      mapHisto_["ThetaResolutionGenVSMu"]->Fill(recMu2,(-genMu.second.Theta()+recMu2.Theta()),+1);
      mapHisto_["CotgThetaResolutionGenVSMu"]->Fill(recMu2,(-cos(genMu.second.Theta())/sin(genMu.second.Theta())
                                                            +cos(recMu2.Theta())/sin(recMu2.Theta())),+1);
      mapHisto_["EtaResolutionGenVSMu"]->Fill(recMu2,(-genMu.second.Eta()+recMu2.Eta()),+1);
      mapHisto_["PhiResolutionGenVSMu"]->Fill(recMu2,(-genMu.second.Phi()+recMu2.Phi()),+1);
    }

    if( simTracks.isValid() ) {
      pair <reco::Particle::LorentzVector, reco::Particle::LorentzVector> simMu = MuScleFitUtils::findSimMuFromRes(evtMC,simTracks);
      reco::Particle::LorentzVector simResonance( simMu.first+simMu.second );
      mapHisto_["SimResonance"]->Fill( simResonance );
      mapHisto_["DeltaSimResonanceMuons"]->Fill( simMu.first, simMu.second );
      mapHisto_["SimResonanceMuons"]->Fill( simMu.first );
      mapHisto_["SimResonanceMuons"]->Fill( simMu.second );

      //first is always mu-, second is always mu+
      if(checkDeltaR(simMu.first,recMu1)){
        mapHisto_["PtResolutionSimVSMu"]->Fill(recMu1,(-simMu.first.Pt()+recMu1.Pt())/simMu.first.Pt(),-1);
        mapHisto_["ThetaResolutionSimVSMu"]->Fill(recMu1,(-simMu.first.Theta()+recMu1.Theta()),-1);
        mapHisto_["CotgThetaResolutionSimVSMu"]->Fill(recMu1,(-cos(simMu.first.Theta())/sin(simMu.first.Theta())
                                                                   +cos(recMu1.Theta())/sin(recMu1.Theta())),-1);
        mapHisto_["EtaResolutionSimVSMu"]->Fill(recMu1,(-simMu.first.Eta()+recMu1.Eta()),-1);
        mapHisto_["PhiResolutionSimVSMu"]->Fill(recMu1,(-simMu.first.Phi()+recMu1.Phi()),-1);
      }
      if(checkDeltaR(simMu.second,recMu2)){
        mapHisto_["PtResolutionSimVSMu"]->Fill(recMu2,(-simMu.second.Pt()+recMu2.Pt())/simMu.first.Pt(),+1);
        mapHisto_["ThetaResolutionSimVSMu"]->Fill(recMu2,(-simMu.second.Theta()+recMu2.Theta()),+1);
        mapHisto_["CotgThetaResolutionSimVSMu"]->Fill(recMu2,(-cos(simMu.second.Theta())/sin(simMu.second.Theta())
                                                                    +cos(recMu2.Theta())/sin(recMu2.Theta())),+1);
        mapHisto_["EtaResolutionSimVSMu"]->Fill(recMu2,(-simMu.second.Eta()+recMu2.Eta()),+1);
        mapHisto_["PhiResolutionSimVSMu"]->Fill(recMu2,(-simMu.second.Phi()+recMu2.Phi()),+1);
      }
    }

    // Fill the mass resolution histograms
    // -----------------------------------
    // check if the recoMuons match the genMuons
    // if( MuScleFitUtils::ResFound && checkDeltaR(simMu.first,recMu1) && checkDeltaR(simMu.second,recMu2) ) {
    if( MuScleFitUtils::ResFound && checkDeltaR(genMu.first,recMu1) && checkDeltaR(genMu.second,recMu2) ) {
      double resMass = (recMu1+recMu2).mass();
      double genMass = (genMu.first + genMu.second).mass();
      // first is always mu-, second is always mu+
      mapHisto_["MassResolution"]->Fill(recMu1, -1, genMu.first, recMu2, +1, genMu.second, resMass, genMass);

      // Fill the reconstructed resonance
      reco::Particle::LorentzVector recoResonance( recMu1+recMu2 );
      mapHisto_["RecoResonance"]->Fill( recoResonance );
      mapHisto_["DeltaRecoResonanceMuons"]->Fill( recMu1, recMu2 );
      mapHisto_["RecoResonanceMuons"]->Fill( recMu1 );
      mapHisto_["RecoResonanceMuons"]->Fill( recMu2 );

      // Fill the covariances histograms
      mapHisto_["Covariances"]->Fill(recMu1, genMu.first, recMu2, genMu.second);

      // Fill the mass resolution (computed from MC), we use the covariance class to compute the variance
      if( genMass != 0 ) {
        double diffMass = (resMass - genMass)/genMass;
        // Fill if for both muons
        double pt1 = recMu1.pt();
        double eta1 = recMu1.eta();
        double pt2 = recMu2.pt();
        double eta2 = recMu2.eta();
        massResolutionVsPtEta_->Fill(pt1, eta1, diffMass, diffMass);
        massResolutionVsPtEta_->Fill(pt2, eta2, diffMass, diffMass);
        // Fill with mass resolution from resolution function
        double massRes = MuScleFitUtils::massResolution(recMu1, recMu2, MuScleFitUtils::parResol);
        // The value given by massRes is already divided by the mass, since the derivative functions have mass at the denominator.
        // We alos take the squared value, since var = sigma^2.
        massResolutionFromFunction_->Fill(pt1, eta1, pow(massRes,2));
        massResolutionFromFunction_->Fill(pt2, eta2, pow(massRes,2));
      }

      // Fill resolution functions for the muons (fill the squared value to make it comparable with the variance)
      mapHisto_["hFunctionResolPt"]->Fill( recMu1, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaPt(recMu1.Pt(), recMu1.Eta(), MuScleFitUtils::parResol ),2), -1 );
      mapHisto_["hFunctionResolCotgTheta"]->Fill( recMu1, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaCotgTh(recMu1.Pt(), recMu1.Eta(), MuScleFitUtils::parResol ),2), -1 );
      mapHisto_["hFunctionResolPhi"]->Fill( recMu1, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaPhi(recMu1.Pt(), recMu1.Eta(), MuScleFitUtils::parResol ),2), -1 );
      mapHisto_["hFunctionResolPt"]->Fill( recMu2, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaPt(recMu2.Pt(), recMu2.Eta(), MuScleFitUtils::parResol ),2), +1 );
      mapHisto_["hFunctionResolCotgTheta"]->Fill( recMu2, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaCotgTh(recMu2.Pt(), recMu2.Eta(), MuScleFitUtils::parResol ),2), +1 );
      mapHisto_["hFunctionResolPhi"]->Fill( recMu2, pow(MuScleFitUtils::resolutionFunctionForVec->sigmaPhi(recMu2.Pt(), recMu2.Eta(), MuScleFitUtils::parResol ),2), +1 );
    }
  } // end if resonance
  else {

    // Loop on the recMuons
    vector<reco::LeafCandidate>::const_iterator recMuon = muons.begin();
    for ( ; recMuon!=muons.end(); ++recMuon ) {  
      int charge = recMuon->charge();

      lorentzVector recMu(recMuon->p4());

      // Find the matching MC muon
      const HepMC::GenEvent* Evt = evtMC->GetEvent();
      //Loop on generated particles
      map<double, lorentzVector> genAssocMap;
      HepMC::GenEvent::particle_const_iterator part = Evt->particles_begin();
      for( ; part!=Evt->particles_end(); ++part ) {
        if (fabs((*part)->pdg_id())==13 && (*part)->status()==1) {
          lorentzVector genMu = (lorentzVector((*part)->momentum().px(),(*part)->momentum().py(),
                                 (*part)->momentum().pz(),(*part)->momentum().e()));

          double deltaR = sqrt(MuScleFitUtils::deltaPhi(recMu.Phi(),genMu.Phi()) * MuScleFitUtils::deltaPhi(recMu.Phi(),genMu.Phi()) +
                               ((recMu.Eta()-genMu.Eta()) * (recMu.Eta()-genMu.Eta())));

          // 13 for the muon (-1) and -13 for the antimuon (+1), thus pdg*charge = -13.
          // Only in this case we consider it matching.
          if( ((*part)->pdg_id())*charge == -13 ) genAssocMap.insert(make_pair(deltaR, genMu));
        }
      }
      // Take the closest in deltaR
      lorentzVector genMu(genAssocMap.begin()->second);

      // Histograms with genParticles characteristics
      // --------------------------------------------

      if(checkDeltaR(genMu,recMu)){
        mapHisto_["PtResolutionGenVSMu"]->Fill(genMu,(-genMu.Pt()+recMu.Pt())/genMu.Pt(),charge);
        mapHisto_["ThetaResolutionGenVSMu"]->Fill(genMu,(-genMu.Theta()+recMu.Theta()),charge);
        mapHisto_["CotgThetaResolutionGenVSMu"]->Fill(genMu,(-cos(genMu.Theta())/sin(genMu.Theta())
                                                             +cos(recMu.Theta())/sin(recMu.Theta())),charge);
        mapHisto_["EtaResolutionGenVSMu"]->Fill(genMu,(-genMu.Eta()+recMu.Eta()),charge);
        mapHisto_["PhiResolutionGenVSMu"]->Fill(genMu,(-genMu.Phi()+recMu.Phi()),charge);
      }

      // Find the matching simMu
      if( simTracks.isValid() ) {
        map<double, lorentzVector> simAssocMap;
        for ( vector<SimTrack>::const_iterator simMuon=simTracks->begin(); simMuon!=simTracks->end(); ++simMuon ) {
          lorentzVector simMu = lorentzVector(simMuon->momentum().px(),simMuon->momentum().py(),
                                              simMuon->momentum().pz(),simMuon->momentum().e());

          double deltaR = sqrt(MuScleFitUtils::deltaPhi(recMu.Phi(),simMu.Phi()) * MuScleFitUtils::deltaPhi(recMu.Phi(),simMu.Phi()) +
                               ((recMu.Eta()-simMu.Eta()) * (recMu.Eta()-simMu.Eta())));

          if( simMuon->charge()*charge == 1 ) simAssocMap.insert(make_pair(deltaR, simMu));
        }
        lorentzVector simMu(genAssocMap.begin()->second);

        //first is always mu-, second is always mu+
        if(checkDeltaR(simMu,recMu)) {
          mapHisto_["PtResolutionSimVSMu"]->Fill(simMu,(-simMu.Pt()+recMu.Pt())/simMu.Pt(),charge);
          mapHisto_["ThetaResolutionSimVSMu"]->Fill(simMu,(-simMu.Theta()+recMu.Theta()),charge);
          mapHisto_["CotgThetaResolutionSimVSMu"]->Fill(simMu,(-cos(simMu.Theta())/sin(simMu.Theta())
                                                               +cos(recMu.Theta())/sin(recMu.Theta())),charge);
          mapHisto_["EtaResolutionSimVSMu"]->Fill(simMu,(-simMu.Eta()+recMu.Eta()),charge);
          mapHisto_["PhiResolutionSimVSMu"]->Fill(simMu,(-simMu.Phi()+recMu.Phi()),charge);
        }
      }
    }
  }

}

void ResolutionAnalyzer::fillHistoMap() {

  outputFile_->cd();
  // Resonances
  // If no Z is required, use a smaller mass range.
  double minMass = 0.;
  double maxMass = 200.;
  if( MuScleFitUtils::resfind[0] != 1 ) {
    maxMass = 30.;
  }
  mapHisto_["GenMother"]               = new HParticle(outputFile_, "GenMother", minMass, maxMass);
  mapHisto_["SimResonance"]            = new HParticle(outputFile_, "SimResonance", minMass, maxMass);
  mapHisto_["RecoResonance"]           = new HParticle(outputFile_, "RecoResonance", minMass, maxMass);

  // Resonance muons
  mapHisto_["GenMotherMuons"]          = new HParticle(outputFile_, "GenMotherMuons", minMass, 1.);
  mapHisto_["SimResonanceMuons"]       = new HParticle(outputFile_, "SimResonanceMuons", minMass, 1.);
  mapHisto_["RecoResonanceMuons"]      = new HParticle(outputFile_, "RecoResonanceMuons", minMass, 1.);

  // Deltas between resonance muons
  mapHisto_["DeltaGenMotherMuons"]     = new HDelta (outputFile_, "DeltaGenMotherMuons");
  mapHisto_["DeltaSimResonanceMuons"]  = new HDelta (outputFile_, "DeltaSimResonanceMuons");
  mapHisto_["DeltaRecoResonanceMuons"] = new HDelta (outputFile_, "DeltaRecoResonanceMuons");

  //   //Reconstructed muon kinematics
  //   //-----------------------------
  //   mapHisto_["hRecBestMu"]             = new HParticle         ("hRecBestMu");
  //   mapHisto_["hRecBestMu_Acc"]         = new HParticle         ("hRecBestMu_Acc"); 
  //   mapHisto_["hDeltaRecBestMu"]        = new HDelta            ("hDeltaRecBestMu");

  //   mapHisto_["hRecBestRes"]            = new HParticle         ("hRecBestRes");
  //   mapHisto_["hRecBestRes_Acc"]        = new HParticle         ("hRecBestRes_Acc"); 
  //   mapHisto_["hRecBestResVSMu"]        = new HMassVSPart       ("hRecBestResVSMu");

  //Resolution VS muon kinematic
  //----------------------------
  mapHisto_["PtResolutionGenVSMu"]        = new HResolutionVSPart (outputFile_, "PtResolutionGenVSMu");
  mapHisto_["PtResolutionSimVSMu"]        = new HResolutionVSPart (outputFile_, "PtResolutionSimVSMu");
  mapHisto_["EtaResolutionGenVSMu"]       = new HResolutionVSPart (outputFile_, "EtaResolutionGenVSMu");
  mapHisto_["EtaResolutionSimVSMu"]       = new HResolutionVSPart (outputFile_, "EtaResolutionSimVSMu");
  mapHisto_["ThetaResolutionGenVSMu"]     = new HResolutionVSPart (outputFile_, "ThetaResolutionGenVSMu");
  mapHisto_["ThetaResolutionSimVSMu"]     = new HResolutionVSPart (outputFile_, "ThetaResolutionSimVSMu");
  mapHisto_["CotgThetaResolutionGenVSMu"] = new HResolutionVSPart (outputFile_, "CotgThetaResolutionGenVSMu");
  mapHisto_["CotgThetaResolutionSimVSMu"] = new HResolutionVSPart (outputFile_, "CotgThetaResolutionSimVSMu");
  mapHisto_["PhiResolutionGenVSMu"]       = new HResolutionVSPart (outputFile_, "PhiResolutionGenVSMu", -0.05, 0.05);
  mapHisto_["PhiResolutionSimVSMu"]       = new HResolutionVSPart (outputFile_, "PhiResolutionSimVSMu");

  // Covariances between muons kinematic quantities
  // ----------------------------------------------
  double ptMax = 200.;
  if (TString(outputFile_->GetName()).Contains("Y")) ptMax = 40.;
  mapHisto_["Covariances"] = new HCovarianceVSParts ( outputFile_, "Covariance", ptMax );

  // Mass resolution
  // ---------------
  mapHisto_["MassResolution"] = new HMassResolutionVSPart (outputFile_,"MassResolution");

  //  mapHisto_["hResolRecoMassVSGenMassVSPt"] = new HResolutionVSPart

  // Mass resolution vs (pt, eta) of the muons from MC
  massResolutionVsPtEta_ = new HCovarianceVSxy ( "Mass", "Mass", 20, 0., ptMax, 20, -3, 3 );
  // Mass resolution vs (pt, eta) of the muons from function
  massResolutionFromFunction_ = new TH2D("massResolutionFromFunction", "mass resolution from function", 40, 0, ptMax, 40, -3, 3);

  // Muons resolutions from resolution functions
  // -------------------------------------------
  mapHisto_["hFunctionResolPt"]        = new HFunctionResolution (outputFile_, "hFunctionResolPt");
  mapHisto_["hFunctionResolCotgTheta"] = new HFunctionResolution (outputFile_, "hFunctionResolCotgTheta");
  mapHisto_["hFunctionResolPhi"]       = new HFunctionResolution (outputFile_, "hFunctionResolPhi");
}

void ResolutionAnalyzer::writeHistoMap() {
  for (map<string, Histograms*>::const_iterator histo=mapHisto_.begin(); 
       histo!=mapHisto_.end(); histo++) {
    (*histo).second->Write();
  }
  outputFile_->cd();
  massResolutionVsPtEta_->Write();
  massResolutionFromFunction_->Write();
}

bool ResolutionAnalyzer::checkDeltaR(const reco::Particle::LorentzVector & genMu, const reco::Particle::LorentzVector & recMu){
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
