#include "HarderAnalysisTools/GenParticleInfo/interface/GenParticleInfo.h"


GenParticleInfo::GenParticleInfo(const reco::Candidate* cand) {

  nDaughters_=0;
  fractionalDiffE_=0;
  double dx=100000,dy=100000,dz=100000;
  for (unsigned i=0; i<cand->numberOfDaughters(); i++) {
    // I observed that a status 3 boson "decays" into dilepton plus a status 2 boson copy.
    // The status 2 boson copy does not have any daughters.
    if (cand->daughter(i)->pdgId()==cand->pdgId()) {
      fractionalDiffE_=(cand->energy()-cand->daughter(i)->energy())/cand->energy();
    } else {
      ++nDaughters_;
      // the status 3 daughters do not get proper vertex info assigned, it seems
      // (tested with Pythia in CMSSW_3_3_2). Need to follow through to the
      // status 2 daughters.
      const reco::Candidate* daughter = cand->daughter(i);
      if (daughter->status()==3 && daughter->numberOfDaughters()>0) {
	dx=daughter->daughter(0)->vertex().X();
	dy=daughter->daughter(0)->vertex().Y();
	dz=daughter->daughter(0)->vertex().Z();
      } else {
	dx=daughter->vertex().X();
	dy=daughter->vertex().Y();
	dz=daughter->vertex().Z();
      }
      dx-=cand->vertex().X();
      dy-=cand->vertex().Y();
      dz-=cand->vertex().Z();
    }
    flightLength3D_=std::sqrt(dx*dx+dy*dy+dz*dz);
    flightLength2D_=std::sqrt(dx*dx+dy*dy);
  }
}


