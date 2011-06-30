// -*- C++ -*-
//
// Package:    TrackingEfficiencyFromCosmics
// Class:      TrackingEfficiencyFromCosmics
// 
/**\class TrackingEfficiencyFromCosmics TrackingEfficiencyFromCosmics.cc Analysis/TrackingEfficiencyFromCosmics/src/TrackingEfficiencyFromCosmics.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,40 3-B32,+41227671551,
//         Created:  Wed May 25 16:44:02 CEST 2011
// $Id: TrackingEfficiencyFromCosmics.cc,v 1.2 2011/05/30 14:47:33 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// #include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "Analysis/TrackingEfficiencyFromCosmics/interface/AssociatorByDeltaR.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlPlots.h"

#include <boost/foreach.hpp>

//
// class declaration
//

class TrackingEfficiencyFromCosmics : public edm::EDAnalyzer {
public:
  explicit TrackingEfficiencyFromCosmics(const edm::ParameterSet&);
  ~TrackingEfficiencyFromCosmics();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  const reco::GenParticle * takeStableMuon(const reco::GenParticleCollection & genParticles);

  void dumpGenParticleInfo(const reco::GenParticle & genParticle);
  void dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber);

  // ----------member data ---------------------------
  // double maxDeltaR_;
  TH1F * hMinDeltaR_, * hMinTrackToGenDeltaR_, * hMinStaMuonToGenDeltaR_;
  reco::Track::TrackQuality quality_;
  bool useMCtruth_;

  std::auto_ptr<AssociatorByDeltaR> associatorByDeltaR_;
  std::auto_ptr<ControlPlots> controlPlotsGeneralTracks_;
  std::auto_ptr<ControlPlots> controlPlotsStandAloneMuons_;
};

TrackingEfficiencyFromCosmics::TrackingEfficiencyFromCosmics(const edm::ParameterSet& iConfig) :
  useMCtruth_(iConfig.getParameter<bool>("UseMCtruth"))
{
  associatorByDeltaR_.reset(new AssociatorByDeltaR(iConfig.getParameter<double>("MaxDeltaR")));
  // maxDeltaR_ = iConfig.getParameter<double>("MaxDeltaR");
}

TrackingEfficiencyFromCosmics::~TrackingEfficiencyFromCosmics() {}

void TrackingEfficiencyFromCosmics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);
  edm::Handle<reco::TrackCollection> staMuons;
  iEvent.getByLabel("standAloneMuons", staMuons);

  // Gen Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    const reco::GenParticle * stableMuon = takeStableMuon(*genParticles);
    BOOST_FOREACH( const reco::Track & track, *tracks ) {
      hMinTrackToGenDeltaR_->Fill(reco::deltaR(*stableMuon, track));
    }
    BOOST_FOREACH( const reco::Track & staMuon, *staMuons ) {
      hMinStaMuonToGenDeltaR_->Fill(reco::deltaR(*stableMuon, staMuon));
    }
  }

  // Fill all control plots
  controlPlotsGeneralTracks_->fillControlPlots(*tracks);
  controlPlotsStandAloneMuons_->fillControlPlots(*staMuons);

  // Association map of StandAloneMuons and TrackerTracks
  std::map<const reco::Track *, const reco::Track *> matchesMap;

  if( staMuons->size() > 0 ) {
    associatorByDeltaR_->fillAssociationMap(tracks, staMuons, matchesMap, hMinDeltaR_);

    std::map<const reco::Track *, const reco::Track *>::const_iterator it = matchesMap.begin();
    for( ; it != matchesMap.end(); ++it ) {
      if( it->second == 0 ) {
        std::cout << "no match found for standAlone with = " << it->first->pt() << std::endl;
      }
      else {
        std::cout << "standAlone with pt = " << it->first->pt() << ", matches with track of pt = " << it->second->pt() << std::endl;
      }
    }
  }



#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

const reco::GenParticle * TrackingEfficiencyFromCosmics::takeStableMuon(const reco::GenParticleCollection & genParticles)
{
  for( reco::GenParticleCollection::const_iterator it = genParticles.begin();
       it != genParticles.end(); ++it ) {
    // dumpGenParticleInfo(*it);
    if( it->status() == 1 ) return( &*it );
  }
  return 0;
}

// ------------ method called once each job just before starting event loop  ------------
void TrackingEfficiencyFromCosmics::beginJob()
{
  edm::Service<TFileService> fileService;
  hMinDeltaR_ = fileService->make<TH1F>("minDeltaR","#Delta R between standalone muon and closest track",500,0,5);
  hMinTrackToGenDeltaR_ = fileService->make<TH1F>("minTrackGenDeltaR","#Delta R between gen muon and closest track",500,0,5);
  hMinStaMuonToGenDeltaR_ = fileService->make<TH1F>("minStaMunoGenDeltaR","#Delta R between gen muon and closest standAloneMuon",500,0,5);
  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "generalTracks"));
  controlPlotsStandAloneMuons_.reset(new ControlPlots(fileService, "standAloneMuons"));
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackingEfficiencyFromCosmics::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void TrackingEfficiencyFromCosmics::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TrackingEfficiencyFromCosmics::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TrackingEfficiencyFromCosmics::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TrackingEfficiencyFromCosmics::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackingEfficiencyFromCosmics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TrackingEfficiencyFromCosmics::dumpGenParticleInfo(const reco::GenParticle & genParticle)
{
  std::cout << "pdgid = " << genParticle.pdgId() << std::endl;
  std::cout << "status = " << genParticle.status() << std::endl;
  std::cout << "pt = " << genParticle.pt() << ", eta = " << genParticle.eta() << ", phi = " << genParticle.phi() << std::endl;
}

void TrackingEfficiencyFromCosmics::dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber)
{
  track.quality(quality_);
  std::stringstream ss;
  ss << "track["<<trackNumber<<"]" << std::endl <<
        "algoName = " << track.algoName() << std::endl <<
        "pt = " << track.pt() << std::endl <<
        "eta = " << track.eta() << std::endl <<
        "phi = " << track.phi() << std::endl <<
        "number of hits = " << track.recHitsSize() << std::endl <<
        "number of valid hits = " << track.found() << std::endl <<
        "number of invalid hits = " << track.lost() << std::endl <<
        "quality = " << track.qualityName(quality_) << std::endl <<
        "d0 = " << track.d0() << std::endl <<
        "d0/d0Error = " << track.d0Error() << std::endl <<
        "dz = " << track.dz() << std::endl <<
        "dz/dzError = " << track.dzError() << std::endl <<
        "chi^2/ndof = " << track.normalizedChi2() << std::endl;
  edm::LogInfo("Demo") << ss.str();
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingEfficiencyFromCosmics);
