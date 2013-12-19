#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include <TObject.h>
#include <TLorentzVector.h>

/**
 * Simple class used to save a gen particle in a root tree. <br>
 * Contains a vector of indexes to daughters and an index to the mother.
 */

/// Class to be saved in the tree
class GenParticle : public TObject
{
public:
  GenParticle() {}

  double pt, eta, phi;
  int charge;
  double dxy, dxyError, dz, dzError; // They have an error because they are the result of a propagation.
  double referencePointRadius, referencePointZ;
  double vx, vy, vz;
  int pid, motherid;

  ClassDef(GenParticle, 1)
};
ClassImp(GenParticle)

#endif // GENPARTICLE_H
