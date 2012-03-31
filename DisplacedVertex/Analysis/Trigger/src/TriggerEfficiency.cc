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
// $Id: TriggerEfficiency.cc,v 1.4 2012/03/07 17:09:13 demattia Exp $
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
#include <TCanvas.h>
#include <TGraph.h>

// #include "Analysis/TrackingEfficiencyFromCosmics/interface/AssociatorByDeltaR.h"
// #include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlPlots.h"
// #include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlDeltaPlots.h"
// #include "Analysis/TrackingEfficiencyFromCosmics/interface/Efficiency.h"
// #include "Analysis/TrackingEfficiencyFromCosmics/interface/EfficiencyTree.h"
#include "Analysis/SmartPropagatorWithIP/interface/SmartPropagatorWithIP.h"
#include "Analysis/Records/interface/SmartPropagatorWithIPComponentsRecord.h"
#include "Analysis/RootTreeProducers/interface/Track.h"
#include "Analysis/RootTreeProducers/interface/GenParticle.h"
#include "Analysis/RootTreeProducers/interface/RootTreeHandler.h"

#include <boost/foreach.hpp>

//
// class declaration
//

class TriggerEfficiency : public edm::EDAnalyzer {
public:
  explicit TriggerEfficiency(const edm::ParameterSet&);
  ~TriggerEfficiency();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  template <class T>
  Track fillTrack(const T & itTrk);
  template <class T>
  GenParticle fillGenParticle(const T & it, const SmartPropagatorWithIP::IP & ip);

  // ----------member data ---------------------------
  edm::ESHandle<TransientTrackBuilder> theB_;
  const SmartPropagatorWithIP * smartPropIP_;
  edm::InputTag trackCollection_;
  bool useMCtruth_;
  boost::shared_ptr<RootTreeHandler> treeHandler_;
  std::vector<Track> tracks_;
  std::vector<GenParticle> genParticles_;
};

TriggerEfficiency::TriggerEfficiency(const edm::ParameterSet& iConfig) :
  smartPropIP_(0),
  trackCollection_(iConfig.getParameter<edm::InputTag>("TrackCollection")),
  useMCtruth_(iConfig.getParameter<bool>("UseMCtruth"))
{
  treeHandler_.reset(new RootTreeHandler(trackCollection_.label()+".root", "L2Muons"));
  treeHandler_->addBranch("tracks", "std::vector<Track>", &tracks_);
  treeHandler_->addBranch("genParticles", "std::vector<GenParticle>", &genParticles_);
}

TriggerEfficiency::~TriggerEfficiency()
{
  treeHandler_->writeTree();
}

template <class T>
Track TriggerEfficiency::fillTrack(const T & itTrk)
{
  Track t;
  t.pt = itTrk->pt();
  t.ptError = itTrk->ptError();
  t.eta = itTrk->eta();
  t.etaError = itTrk->etaError();
  t.phi = itTrk->phi();
  t.phiError = itTrk->phiError();
  t.charge = itTrk->charge();
  t.dxy = itTrk->dxy();
  t.dxyError = itTrk->dxyError();
  t.dz = itTrk->dz();
  t.dzError = itTrk->dzError();
  t.vx = itTrk->vx();
  t.vy = itTrk->vy();
  t.vz = itTrk->vz();
  t.chi2 = itTrk->chi2();
  t.normalizedChi2 = itTrk->normalizedChi2();
  t.referencePointRadius = sqrt(pow(itTrk->referencePoint().x(),2) + pow(itTrk->referencePoint().y(),2));
  t.referencePointZ = itTrk->referencePoint().z();
  t.nHits = itTrk->recHitsSize();
  t.nValidHits = itTrk->found();
  t.nValidPlusInvalidHits = itTrk->found()+itTrk->lost();
  t.innermostHitRadius = itTrk->innerPosition().r();
  t.innermostHitZ = itTrk->innerPosition().z();
  t.muonStationsWithAnyHits = itTrk->hitPattern().muonStationsWithAnyHits();
  t.dtStationsWithAnyHits = itTrk->hitPattern().dtStationsWithAnyHits();
  t.cscStationsWithAnyHits = itTrk->hitPattern().cscStationsWithAnyHits();
  return t;
}

template <class T>
GenParticle TriggerEfficiency::fillGenParticle(const T & it, const SmartPropagatorWithIP::IP & ip)
{
  GenParticle p;
  p.pt = it->pt();
  p.eta = it->eta();
  p.phi = it->phi();
  p.charge = it->charge();
  p.dxy = ip.dxyValue;
  p.dz = ip.dzValue;
  p.referencePointRadius = 0.;
  p.referencePointZ = 0.;
  p.vx = it->vertex().x();
  p.vy = it->vertex().y();
  p.vz = it->vertex().z();
  p.pid = it->pdgId();
  return p;
}

void TriggerEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // eventNum_ = iEvent.id().event();
  // std::cout << "event = " << eventNum_ << std::endl;

  // Load the transient track builder
  // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB_);

  // Load the propagator and IP calculator
  edm::ESHandle<Propagator> smartPropIPHandle;
  iSetup.get<SmartPropagatorWithIPComponentsRecord>().get("SmartPropagatorWithIP", smartPropIPHandle);
  smartPropIP_ = dynamic_cast<const SmartPropagatorWithIP*>(&*smartPropIPHandle);

  try {
    edm::Handle<reco::TrackCollection> allTracks;
    iEvent.getByLabel(trackCollection_, allTracks);
    reco::TrackCollection::const_iterator itTrk = allTracks->begin();
    for( ; itTrk != allTracks->end(); ++itTrk ) {
      tracks_.push_back(fillTrack(itTrk));
    }
  }
  catch (cms::Exception & ex) {
    std::cerr << ex;
  }

  // Gen Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
      if( it->status() == 1 && fabs(it->pdgId()) == 13 ) {
        // Compute impact parameters for generator particle
        SmartPropagatorWithIP::IP ip;
        if( sqrt(pow(it->vertex().x(),2) + pow(it->vertex().y(),2)) < 110. && fabs(it->vertex().z()) < 280. ) {
          ip = smartPropIP_->computeGenImpactParametersInsideTkVol( *it, it->vertex(), it->charge(), GlobalPoint(0,0,0) );
        }
        else {
          ip = smartPropIP_->computeGenImpactParametersOutsideTkVol( *it, it->vertex(), it->charge(), GlobalPoint(0,0,0) );
        }

        genParticles_.push_back(fillGenParticle(it, ip));
      }
    }
  }

  treeHandler_->saveToTree(iEvent.id().event(), iEvent.run());
  tracks_.clear();
  genParticles_.clear();
}

// ------------ method called once each job just before starting event loop  ------------
void TriggerEfficiency::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerEfficiency::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void TriggerEfficiency::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TriggerEfficiency::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TriggerEfficiency::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TriggerEfficiency::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEfficiency);
