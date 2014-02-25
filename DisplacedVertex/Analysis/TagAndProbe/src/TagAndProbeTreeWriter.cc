// -*- C++ -*-
//
// Package:    TagAndProbeTreeWriter
// Class:      TagAndProbeTreeWriter
//
/**\class TagAndProbeTreeWriter TagAndProbeTreeWriter.cc Analysis/TagAndProbeTreeWriter/src/TagAndProbeTreeWriter.cc

 Description: Can be used to save the information in a tree. The tracks can then be accessed with the macros under test/Macros.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,40 1-A11,
//         Created:  Thu Nov 10 19:02:00 CEST 2013
// $Id: TagAndProbeTreeWriter.cc,v 1.4 2013/11/10 19:02:00 demattia Exp $
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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include <TVector3.h>

#include "Analysis/RootTreeProducers/interface/RootTreeHandler.h"
#include "Analysis/RootTreeProducers/interface/Track.h"
#include "Analysis/RootTreeProducers/interface/GenParticle.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <boost/foreach.hpp>
#include <string>
#include <map>

//
// class declaration
//

class TagAndProbeTreeWriter : public edm::EDAnalyzer {
public:
  explicit TagAndProbeTreeWriter(const edm::ParameterSet&);
  ~TagAndProbeTreeWriter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // Helpers to fill the Track and GenParticle objects.
  void fillTrackToTreeTrack( Track & treeTrack, const reco::Track & track );
  // For the tree
  template <class T>
  void fillTreeTracks( const T & collection, std::vector<Track> & tracks, const double & minPt );
  template <class T>
  GenParticle fillGenParticle(const T & it);
  template <class T>
  Track fillTreeTriggerTrack( const T & triggerObject );
  template <class T>
  double computeIsolation(const T & track);

  // ----------member data ---------------------------
  bool useMCtruth_;

  TString outputName_;
  edm::InputTag muonCollection_;
  edm::InputTag standAloneMuonCollection_;
  edm::InputTag trackCollectionName_;
  const reco::TrackCollection * trackCollection_;
  boost::shared_ptr<RootTreeHandler> treeHandler_;
  std::vector<Track> tracks_;
  std::vector<Track> muons_;
  std::vector<Track> standAloneMuons_;
  std::vector<GenParticle> genParticles_;
  std::vector<std::string> selectedTriggerNames_;
  std::vector<std::string> triggerNamesPassed_;
  std::vector<std::string> selectedFilterNames_;
  std::map<std::string, std::vector<Track> > triggerFilterObjectsMap_;
  double minTrackPt_;
  double minMuonPt_;
  double minStandAloneMuonPt_;
};

TagAndProbeTreeWriter::TagAndProbeTreeWriter(const edm::ParameterSet& iConfig) :
  useMCtruth_(iConfig.getParameter<bool>("UseMCtruth")),
  outputName_(iConfig.getParameter<std::string>("OutputName")),
  muonCollection_(iConfig.getParameter<edm::InputTag>("MuonCollection")),
  standAloneMuonCollection_(iConfig.getParameter<edm::InputTag>("StandAloneMuonCollection")),
  // triggerMuonCollection_(iConfig.getParameter<edm::InputTag>("TriggerMuonCollection")),
  trackCollectionName_(iConfig.getParameter<edm::InputTag>("TrackCollection")),
  selectedTriggerNames_(iConfig.getParameter<std::vector<std::string> >("TriggerNames")),
  selectedFilterNames_(iConfig.getParameter<std::vector<std::string> >("FilterNames")),
  minTrackPt_(iConfig.getParameter<double>("MinTrackPt")),
  minMuonPt_(iConfig.getParameter<double>("MinMuonPt")),
  minStandAloneMuonPt_(iConfig.getParameter<double>("MinStandAloneMuonPt"))
{
  treeHandler_.reset(new RootTreeHandler(outputName_, "MuonsTree"));
  treeHandler_->addBranch("tracks", "std::vector<Track>", &tracks_);
  treeHandler_->addBranch("muons", "std::vector<Track>", &muons_);
  treeHandler_->addBranch("refittedStandAloneMuons", "std::vector<Track>", &standAloneMuons_);
  treeHandler_->addBranch("genParticles", "std::vector<GenParticle>", &genParticles_);
  treeHandler_->addBranch("triggerNames", "std::vector<std::string> >", &triggerNamesPassed_);

  // Create a separate collection for each triggerFilter. We were doing it with a map<string, vector<Track> >,
  // but had to change to a simpler way due to unpredictable ROOT behavior.
  for( auto filterName = selectedFilterNames_.begin(); filterName != selectedFilterNames_.end(); ++filterName ) {
    triggerFilterObjectsMap_.insert(std::make_pair(*filterName, std::vector<Track>()));
    treeHandler_->addBranch(*filterName+"Objects", "std::vector<Track>", &(triggerFilterObjectsMap_.at(*filterName)));
  }
}

TagAndProbeTreeWriter::~TagAndProbeTreeWriter()
{
  std::cout << "Saving trees" << std::endl;
  treeHandler_->writeTree();
}

void TagAndProbeTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Trigger information
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigResultsTag("TriggerResults","","HLT");
  iEvent.getByLabel(trigResultsTag,trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

  // Find the trigger names
  for( auto it = trigNames.triggerNames().begin(); it != trigNames.triggerNames().end(); ++it ) {
    for( auto name = selectedTriggerNames_.begin(); name != selectedTriggerNames_.end(); ++name ) {
      if( it->find(*name) != std::string::npos ) {
        if( trigResults->accept(trigNames.triggerIndex(*it)) ) {
          triggerNamesPassed_.push_back(*name);
        }
      }
    }
  }
  if( triggerNamesPassed_.size() == 0 ) return;

  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent;
  iEvent.getByLabel(trigEventTag, trigEvent);

  for( auto filterName = selectedFilterNames_.begin(); filterName != selectedFilterNames_.end(); ++filterName ) {
    trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(*filterName,"",trigEventTag.process()));
    if(filterIndex<trigEvent->sizeFilters()){
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex);
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop on the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
        const trigger::TriggerObject& obj = trigObjColl[*keyIt];
        triggerFilterObjectsMap_.at(*filterName).push_back(fillTreeTriggerTrack(obj));
      }
    }
  }

  // Fill the tracks
  try {
    edm::Handle<reco::TrackCollection> allTracks;
    iEvent.getByLabel(trackCollectionName_, allTracks);
    trackCollection_ = allTracks.product();
    fillTreeTracks(*(allTracks.product()), tracks_, minTrackPt_);
  }
  catch (cms::Exception & ex) {
    std::cerr << ex;
  }

  // Fill the muons
  try {
    edm::Handle<reco::TrackCollection> muons;
    iEvent.getByLabel(muonCollection_, muons);
    fillTreeTracks(*(muons.product()), muons_, minMuonPt_);
  }
  catch (cms::Exception & ex) {
    std::cerr << ex;
  }

  // Fill the standAloneMuons
  try {
    edm::Handle<reco::TrackCollection> standAloneMuons;
    iEvent.getByLabel(standAloneMuonCollection_, standAloneMuons);
    fillTreeTracks(*(standAloneMuons.product()), standAloneMuons_, minStandAloneMuonPt_);
  }
  catch (cms::Exception & ex) {
    std::cerr << ex;
  }


  // Gen Particles. Note: selecting stable muons only
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
      if( it->status() == 1 && fabs(it->pdgId()) == 13 ) {
        genParticles_.push_back(fillGenParticle(it));
      }
    }
  }

  // Fill the trees
  treeHandler_->saveToTree(iEvent.id().event(), iEvent.run());
  tracks_.clear();
  muons_.clear();
  standAloneMuons_.clear();
  genParticles_.clear();
  triggerNamesPassed_.clear();
  for( auto filterName = selectedFilterNames_.begin(); filterName != selectedFilterNames_.end(); ++filterName ) {
    triggerFilterObjectsMap_.at(*filterName).clear();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void TagAndProbeTreeWriter::beginJob()
{
  // edm::Service<TFileService> fileService;
}

// ------------ method called once each job just after ending the event loop  ------------
void TagAndProbeTreeWriter::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void TagAndProbeTreeWriter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TagAndProbeTreeWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TagAndProbeTreeWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TagAndProbeTreeWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TagAndProbeTreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template <class T>
Track TagAndProbeTreeWriter::fillTreeTriggerTrack( const T & triggerObject )
{
  Track treeTrack;
  treeTrack.pt = triggerObject.pt();
  treeTrack.eta = triggerObject.eta();
  treeTrack.phi = triggerObject.phi();
  return treeTrack;
}

template <class T>
void TagAndProbeTreeWriter::fillTreeTracks( const T & collection, std::vector<Track> & tracks, const double & minPt )
{
  typename T::const_iterator it = collection.begin();
  for( ; it != collection.end(); ++it ) {
    if( it->pt() < minPt ) continue;
    Track treeTrack;
    fillTrackToTreeTrack(treeTrack, *it);
    tracks.push_back(treeTrack);
  }
}

template <class T>
GenParticle TagAndProbeTreeWriter::fillGenParticle(const T & it)
{
  GenParticle p;
  p.pt = it->pt();
  p.eta = it->eta();
  p.phi = it->phi();
  p.charge = it->charge();
  p.dxy = 0.;
  p.dz = 0.;
  p.referencePointRadius = 0.;
  p.referencePointZ = 0.;
  p.vx = it->vertex().x();
  p.vy = it->vertex().y();
  p.vz = it->vertex().z();
  p.pid = it->pdgId();
  return p;
}

template <class T>
double TagAndProbeTreeWriter::computeIsolation(const T & track)
{
  double isolation = 0.;
  for( auto tk = trackCollection_->begin(); tk != trackCollection_->end(); ++tk ) {
    if( tk->pt() <= 1. ) continue;
    double deltaR = reco::deltaR(track, *tk);
    if( deltaR < 0.3 && deltaR > 0.03 ) {
      isolation += tk->pt();
    }
  }
  return isolation;
}

void TagAndProbeTreeWriter::fillTrackToTreeTrack( Track & treeTrack, const reco::Track & track )
{
  treeTrack.pt = track.pt();
  treeTrack.ptError = track.ptError();
  treeTrack.eta = track.eta();
  treeTrack.etaError = track.etaError();
  treeTrack.phi = track.phi();
  treeTrack.phiError = track.phiError();
  treeTrack.charge = track.charge();
  treeTrack.chi2 = track.chi2();
  treeTrack.normalizedChi2 = track.normalizedChi2();
  treeTrack.isolation = computeIsolation(track);
  treeTrack.trackerLayersWithMeasurement = track.hitPattern().trackerLayersWithMeasurement();
  treeTrack.pixelLayersWithMeasurement = track.hitPattern().pixelLayersWithMeasurement();
  treeTrack.numberOfValidStripLayersWithMonoAndStereo = track.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
  treeTrack.referencePointRadius = std::sqrt(std::pow(track.referencePoint().x(),2)+std::pow(track.referencePoint().y(),2));
  treeTrack.referencePointZ = track.referencePoint().z();
  // treeTrack.nHits = track.recHitsSize();
  treeTrack.nValidHits = track.found();
  treeTrack.nValidPlusInvalidHits = track.found()+track.lost();
  treeTrack.trackAlgorithm = track.algo();
  reco::TrackBase::TrackQuality trackQualityHighPurity = reco::TrackBase::qualityByName("highPurity");
  treeTrack.trackQuality = track.quality(trackQualityHighPurity);
//  try {
//    treeTrack.innermostHitRadius = track.innerPosition().r();
//    treeTrack.innermostHitZ = track.innerPosition().z();
//    treeTrack.muonStationsWithAnyHits = track.hitPattern().muonStationsWithAnyHits();
//    treeTrack.dtStationsWithValidHits = track.hitPattern().dtStationsWithValidHits();
//    treeTrack.cscStationsWithValidHits = track.hitPattern().cscStationsWithValidHits();
//  }
//  catch (cms::Exception & ex) {
//    // Not present in AOD
//    // std::cerr << ex;
//  }
  treeTrack.dxy = track.dxy();
  treeTrack.dxyError = track.dxyError();
  treeTrack.dz = track.dz();
  treeTrack.dzError = track.dzError();

  // Removed
  treeTrack.nHits = 0.;
  treeTrack.innermostHitRadius = 0.;
  treeTrack.innermostHitZ = 0.;
  treeTrack.muonStationsWithAnyHits = 0.;
  treeTrack.dtStationsWithValidHits = 0.;
  treeTrack.cscStationsWithValidHits = 0.;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagAndProbeTreeWriter);
