#ifndef TRACK_H
#define TRACK_H

#include <TObject.h>
#include <TLorentzVector.h>

/**
 * Simple class used to save a track in a root tree. <br>
 */

/// Class to be saved in the tree
class Track : public TObject
{
public:
  Track() {}

  double pt, ptError, eta, etaError, phi, phiError;
  int charge;
  double dxy, dxyError, dz, dzError;
  double vx, vy, vz;
  double chi2, normalizedChi2;
  double referencePointRadius, referencePointZ;
  int nHits, nValidHits, nValidPlusInvalidHits;
  double innermostHitRadius, innermostHitZ;
  int muonStationsWithAnyHits;
  int dtStationsWithAnyHits;
  int cscStationsWithAnyHits;
  int trackAlgorithm;
  bool trackQuality;

  ClassDef(Track, 1)
};
ClassImp(Track)

#endif // TRACK_H
