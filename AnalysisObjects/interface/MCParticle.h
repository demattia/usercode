#ifndef MCPARTICLE_H
#define MCPARTICLE_H

#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

/**
 *
 * Stores a MCParticle. Inherits from BaseParticle and
 * adds the methods for mass, pid and mother's pid.
 *
 * Author M. De Mattia - 10/11/2007
 *
 */

namespace anaobj {

  class MCParticle : public BaseParticle {
  public:
    /// Default constructor, only needed for classes.h
    MCParticle() : BaseParticle() {
      mass_ = 0.;
      pid_ = 0;
      mPid_ = 0;
    }
    MCParticle( const double & PT, const double & ETA, const double & PHI, const double & MASS, const int PID, const int MPID ) : BaseParticle( PT, ETA, PHI ) {
      mass_ = MASS;
      pid_ = PID;
      mPid_ = MPID;
    }
    double mass() const { return mass_; }
    int pid() const { return pid_; }
    int mPid() const { return mPid_; }
    void setMass( const double & MASS ) { mass_ = MASS; }
    void setPid( const int PID ) { pid_ = PID; }
    void setMpid( const int MPID ) { mPid_ = MPID; }
  protected:
    double mass_;
    int pid_;
    int mPid_;
  };

  typedef std::vector<MCParticle> MCParticleCollection;

}

#endif //MCPARTICLE_H
