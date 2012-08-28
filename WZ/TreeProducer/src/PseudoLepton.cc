#include "TreeProducer/TreeProducer/interface/PseudoLepton.h"

/// default constructor
PseudoLepton::PseudoLepton() :
  reco::Track(),
  genPartSize_(0),
  isStandAloneMuon_(false),
  isGlobalMuon_(false),
  isTrackerMuon_(false),
  trackIso_(0.),
  ecalIso_(0.),
  hcalIso_(0.),
  normChi2_(0.),
  numberOfValidTrackerHits_(0.),
  numberOfValidPixelHits_(0.),
  numberOfValidMuonHits_(0.),
  numberOfMatchedStations_(0.)
{
}


/// constructor from reco::Track
PseudoLepton::PseudoLepton(const reco::Track & aTrack) :
  reco::Track(aTrack),
  genPartSize_(0),
  isStandAloneMuon_(false),
  isGlobalMuon_(false),
  isTrackerMuon_(false),
  trackIso_(0.),
  ecalIso_(0.),
  hcalIso_(0.),
  normChi2_(0.),
  numberOfValidTrackerHits_(0.),
  numberOfValidPixelHits_(0.),
  numberOfValidMuonHits_(0.),
  numberOfMatchedStations_(0.)
{
  p4_.SetCoordinates(aTrack.momentum().x(),aTrack.momentum().y(),
		     aTrack.momentum().z(),sqrt(aTrack.momentum().mag2()+0.14*0.14));
}


/// destructor
PseudoLepton::~PseudoLepton() {
}


