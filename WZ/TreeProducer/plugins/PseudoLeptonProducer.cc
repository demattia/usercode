// -*- C++ -*-
//
// Package:    PseudoLeptonProducer
// Class:      PseudoLeptonProducer
// 
/**\class PseudoLeptonProducer PseudoLeptonProducer.cc HarderAnalysis/DisplacedDileptons/plugins/PseudoLeptonProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kristian Harder
//         Created:  Tue Jul 02 21:16:02 GMT 2011
// $Id: PseudoLeptonProducer.cc,v 1.1 2012/05/24 11:53:57 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

#include "TreeProducer/TreeProducer/interface/PseudoLepton.h"


//
// class declaration
//

class PseudoLeptonProducer : public edm::EDProducer {
public:
  explicit PseudoLeptonProducer(const edm::ParameterSet&);
  ~PseudoLeptonProducer();
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  // ----------member data ---------------------------

  edm::InputTag trackTag_;
  edm::InputTag muonTag_;
  edm::InputTag patMuonTag_;

  double minPt_;
  bool useStandAloneMuons_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PseudoLeptonProducer::PseudoLeptonProducer(const edm::ParameterSet& iConfig) :
  trackTag_(iConfig.getParameter<edm::InputTag>("trackSrc")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonSrc")),
  patMuonTag_(iConfig.getParameter<edm::InputTag>("patMuonSrc")),
  minPt_(iConfig.getParameter<double>("minPt")),
  useStandAloneMuons_(iConfig.getParameter<bool>("useStandAloneMuons"))
{
  produces<std::vector<PseudoLepton> >();
}


PseudoLeptonProducer::~PseudoLeptonProducer()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
PseudoLeptonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // access track collection
   Handle<RefVector<reco::TrackCollection> > tracks;
   iEvent.getByLabel(trackTag_,tracks);

   // copy tracks into PseudoLepton collection
   std::vector<PseudoLepton>* leptons = new std::vector<PseudoLepton>();
   for (unsigned i=0; i<tracks->size(); i++) {
     const Ref<reco::TrackCollection> trk = (*tracks)[i];
     PseudoLepton pl(*trk);
     pl.isStandAloneMuon(false);
     pl.isGlobalMuon(false);
     pl.isTrackerMuon(false);
     leptons->push_back(pl);
   }

   if( useStandAloneMuons_ ) {
     // add high p_t standalone muons to list of high_pt tracks
     Handle<reco::TrackCollection> muons;
     iEvent.getByLabel(muonTag_,muons);
     for (reco::TrackCollection::const_iterator trk=muons->begin(); trk!=muons->end(); trk++) {
       if (trk->pt()>minPt_) {
	 PseudoLepton pl(*trk);
	 pl.isStandAloneMuon(true);
	 pl.isGlobalMuon(false);
	 pl.isTrackerMuon(false);
	 leptons->push_back(pl);
       }
     }
   }

   // add high p_t globalMuons to list of high_pt tracks
   Handle<pat::MuonCollection> patMuons;
   iEvent.getByLabel(patMuonTag_,patMuons);
   for (pat::MuonCollection::const_iterator trk=patMuons->begin(); trk!=patMuons->end(); trk++) {
     if (trk->pt()>minPt_ && (trk->isGlobalMuon() || trk->isTrackerMuon())) {
       PseudoLepton pl(*(trk->innerTrack()));
       pl.isStandAloneMuon(false);
       if( trk->isTrackerMuon() )  pl.isTrackerMuon(true);
       pl.setTrackIso(trk->trackIso());
       pl.setEcalIso(trk->ecalIso());
       pl.setHcalIso(trk->hcalIso());

       if( trk->isGlobalMuon() )  {
	 pl.isGlobalMuon(true);
	 pl.setNormChi2( trk->globalTrack()->normalizedChi2() );
	 pl.setNumberOfValidTrackerHits( trk->innerTrack()->hitPattern().numberOfValidHits() );
	 pl.setNumberOfValidPixelHits( trk->innerTrack()->hitPattern().numberOfValidPixelHits() );
	 pl.setNumberOfValidMuonHits( trk->globalTrack()->hitPattern().numberOfValidMuonHits() );
	 pl.setNumberOfMatchedStations( trk->numberOfMatchedStations() );
       }

       leptons->push_back(pl);
     }
   }

   std::auto_ptr<std::vector<PseudoLepton> > ptr(leptons);
   iEvent.put(ptr);

}

//define this as a plug-in
DEFINE_FWK_MODULE(PseudoLeptonProducer);
