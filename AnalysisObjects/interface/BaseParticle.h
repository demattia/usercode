#ifndef BASEPARTICLE_H
#define BASEPARTICLE_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseAll.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Simple BaseParticle class used as base class for
   * MCParticles.
   * It is exactly like BaseJet, but et -> pt.
   *
   * Author M. De Mattia - 8/11/2007
   *
   */

  class BaseParticle : public BaseAll {
  public:
    BaseParticle( const double & PT, const double & ETA, const double & PHI ) : BaseAll( ETA, PHI ) {
      pt_ = PT;
    }
    /// Default constructor, only needed for classes.h
    BaseParticle() : BaseAll( 0., 0. ) {
      pt_ = 0.;
    }
    double pt() const {
      return pt_;
    }
    float px() const {
      return pt_*std::cos(phi_);
    }
    float py() const {
      return pt_*std::sin(phi_);
    }
    float pz() const {
      return pt_*(1-std::exp(-2*eta_))/(2*std::exp(-eta_));
    }
    float p() const {
      return std::sqrt(std::pow(px(),2)+std::pow(py(),2)+std::pow(pz(),2));
    }
    /// To sort the BaseParticles
    bool operator< ( const BaseParticle& b ) const {
      return pt() < b.pt();
    }
  protected:
    double pt_;
  };

}

#endif // BASEPARTICLE_H
