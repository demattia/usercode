#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <map>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include "Analysis/TrackingEfficiencyFromCosmics/interface/TreeTrack.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"
//#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/**
  * Namespace containing some utility functions. They are declared static to limit their scope to this file and avoid linking errors
  * when including it in multiple plugins.
  */
namespace utils
{
  TH1F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name,
                              const TString & nameForTitle, const TString axisTitle, const TString & unit,
		       const int bins, const double & min, const double & max);

  TH2F * bookHistogram(edm::Service<TFileService> & fileService, const TString & name, const TString & nameForTitle,
                              const TString xAxisTitle, const TString & xUnit,
                              const int xBins, const double & xMin, const double & xMax,
                              const TString yAxisTitle, const TString & yUnit,
                              const int yBins, const double & yMin, const double & yMax);

  /// Polynomial function fitted to the dxy_error vs Pt distribution
  double polynomial( const double & pt );

  /// Returns the maximum dxy error value allowed for the cuts
  double dxyErrMax( const double & pt );

  /// Fills the TreeTrack
  void fillTrackToTreeTrack(TreeTrack & treeTrack, const reco::Track & track);
  void fillTrackToTreeTrack(TreeTrack & treeTrack, const reco::Track & track, const SmartPropagatorWithIP::IP & ip);
  void fillGenToTreeTrack(TreeTrack & treeTrack, const reco::GenParticle & track, const SmartPropagatorWithIP::IP & genIp);
}

#endif
