#ifndef CheckHitPattern_H
#define CheckHitPattern_H

/*
 * Determine if a track has hits in front of its assumed production point.
 * Also determine if it misses hits between its assumed production point and its innermost hit.
 *
 * this class was written by Ian Tomalin, RAL
 * small modification by Kristian Harder: force geometry initialisation in constructor
 */

#define DEBUG_CHECKHITPATTERN

// standard EDAnalyser include files
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include <utility>
#include <map>

class DetId;

class CheckHitPattern {

public:

  CheckHitPattern(const edm::EventSetup& iSetup);
  
  ~CheckHitPattern() {}

  // Check if hit pattern of this track is consistent with it being produced
  // at given vertex. Pair.first gives number of hits on track in front of vertex.
  // Pair.second gives number of missing hits between vertex and innermost hit
  // on track.
  std::pair<unsigned int, unsigned int>
    analyze(const reco::Track& track, const TransientVertex& vert);

  // Print hit pattern on track
  void print(const reco::Track& track) const;

private:
  // Return a pair<uint32, uint32> consisting of the numbers used by HitPattern to 
  // identify subdetector and layer number respectively.
  typedef std::pair<uint32_t, uint32_t> DetInfo;
  static DetInfo interpretDetId(DetId detId);

  // Return a bool indicating if a given subdetector is in the barrel.
  static bool barrel(uint32_t subDet);

  void print(const reco::HitPattern& hp) const;

private:
  // For a given subdetector & layer number, this stores the minimum and maximum
  // r (or z) values if it is barrel (or endcap) respectively.
  typedef std::map< DetInfo, std::pair< double, double> > RZrangeMap;
  static RZrangeMap rangeRorZ_;
};

#endif
