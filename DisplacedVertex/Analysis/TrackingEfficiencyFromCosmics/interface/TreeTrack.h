#ifndef TREETRACK_H
#define TREETRACK_H

#include <TObject.h>
#include <TLorentzVector.h>

/**
 * Simple class used to save a track in a root tree. <br>
 * Includes the information on the genTrack.
 */

/// Class to be saved in the tree
class TreeTrack : public TObject
{
public:
  TreeTrack() {}

  double pt, ptError, eta, etaError, phi, phiError;
  int charge;
  double dxy, dxyError, dz, dzError;
  double vx, vy, vz;
  double chi2, normalizedChi2;
  double referencePointRadius, referencePointZ;
  int nHits, nValidHits, nValidPlusInvalidHits;
  double innermostHitRadius, innermostHitZ;

  double genPt, genEta, genPhi;
  int genCharge;
  double genDxy, genDxyError, genDz, genDzError;
  double genVx, genVy, genVz;

  ClassDef(TreeTrack, 1)
};
ClassImp(TreeTrack)

#endif // TREETRACK_H
