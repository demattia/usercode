#ifndef SIMPLETRACK_H
#define SIMPLETRACK_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Used for tracks with vertex.
   * Inherits from BaseParticle.
   *
   * Author M. De Mattia - 10/4/2008
   *
   * Method zError()
   * Modified by Roberto Casagrande - 15/04/2008
   *
   */
  
  class SimpleTrack : public BaseParticle {
  public:
    SimpleTrack( const double & PT, const double & ETA, const double & PHI, const double & Z, const double & ZERROR ) : BaseParticle( PT, ETA, PHI ) {
      z_ = Z;
      zerror_ = ZERROR;
    }
      // Default constructor, only needed for classes.h
      SimpleTrack() : BaseParticle( 0., 0., 0. ) {
	z_ = 0.;
	zerror_ = 0.;
      }
	double z() const { return z_; }
	double zError() const { return zerror_; }
	void setZ( const double & Z ) { z_ = Z; }
	void setZError( const double & ZERROR ) { zerror_ = ZERROR; }
	
  protected:
	double z_;
	double zerror_;
  };
  
  typedef std::vector<SimpleTrack> SimpleTrackCollection;
  
}

#endif // SIMPLETRACK_H
