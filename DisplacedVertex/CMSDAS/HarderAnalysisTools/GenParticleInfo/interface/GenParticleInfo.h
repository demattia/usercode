#ifndef GenParticleInfo_h
#define GenParticleInfo_h 1

#include "DataFormats/Candidate/interface/Candidate.h"


class GenParticleInfo {
  
public:
 
  GenParticleInfo(const reco::Candidate* cand);
  ~GenParticleInfo() { };
  inline double flightLength3D() { return flightLength3D_; };
  inline double flightLength2D() { return flightLength2D_; };
  inline double nDaughters() { return nDaughters_; };
  inline double fractionalDiffE() { return fractionalDiffE_; };

private:

  double flightLength3D_;
  double flightLength2D_;
  unsigned nDaughters_;
  double fractionalDiffE_;

} ;

#endif
