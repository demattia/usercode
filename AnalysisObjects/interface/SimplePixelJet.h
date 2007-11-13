#ifndef SIMPLEPIXELJET_H
#define SIMPLEPIXELJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

#include <cmath>
#include <vector>

/**
 *
 * Used for pixel jets. Same as a PixelJet but without the
 * vector of references to pixel tracks.
 * Inherits from BaseParticle, since it has pt, not et.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class SimplePixelJet : public BaseParticle {
  public:
    SimplePixelJet( const double & PT, const double & ETA, const double & PHI, const double & Z, const int TKNUM ) : BaseParticle( PT, ETA, PHI ) {
      z_ = Z;
      tkNum_ = TKNUM;
    }
    // Default constructor, only needed for classes.h
    SimplePixelJet() : BaseParticle( 0., 0., 0. ) {
      z_ = 0.;
      tkNum_ = 0;
    }
    double z() const { return z_; }
    int tkNum() const { return tkNum_; }
    void setZ( const double & Z ) { z_ = Z; }
    void setTkNum( const int TKNUM ) { tkNum_ = TKNUM; }
  protected:
    int tkNum_;
    double z_;
  };

  typedef std::vector<SimplePixelJet> SimplePixelJetCollection;

}

#endif // SIMPLEPIXELJET_H
