#include "TreeProducer/TreeProducer/interface/PseudoLepton.h"

/// default constructor
PseudoLepton::PseudoLepton() :
  reco::Track(),
  genPartSize_(0),
  isStandAloneMuon_(false),
  isGlobalMuon_(false),
  isTrackerMuon_(false)
{
}


/// constructor from reco::Track
PseudoLepton::PseudoLepton(const reco::Track & aTrack) :
  reco::Track(aTrack),
  genPartSize_(0),
  isStandAloneMuon_(false),
  isGlobalMuon_(false),
  isTrackerMuon_(false)
{
  p4_.SetCoordinates(aTrack.momentum().x(),aTrack.momentum().y(),
		     aTrack.momentum().z(),sqrt(aTrack.momentum().mag2()+0.14*0.14));
}


/// destructor
PseudoLepton::~PseudoLepton() {
}


