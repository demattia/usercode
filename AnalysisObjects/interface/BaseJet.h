#ifndef BASEJET_H
#define BASEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseAll.h"

#include <cmath>
#include <vector>

/**
 *
 * Simple BaseJet class used for L1Jets and as base class
 * for IC5 jets.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class BaseJet : public BaseAll {
  public:
    BaseJet( const double & ET, const double & ETA, const double & PHI ) : BaseAll( ETA, PHI ) {
      et_ = ET;
    }
    /// Default constructor, only needed for classes.h
    BaseJet() : BaseAll( 0., 0. ) {
      et_ = 0.;
    }
    double et() const {
      return et_;
    }
    float ex() const {
      return et_*std::cos(phi_);
    }
    float ey() const {
      return et_*std::sin(phi_);
    }
    float ez() const {
      return et_*(1-std::exp(-2*eta_))/(2*std::exp(-eta_));
    }
    float e() const {
      return std::sqrt(std::pow(ex(),2)+std::pow(ey(),2)+std::pow(ez(),2));
    }
    /// To sort the BaseJets
    bool operator< ( const BaseJet& b ) const {
      return et() < b.et();
    }
  protected:
    double et_;
  };

  typedef std::vector<BaseJet> BaseJetCollection;

}

#endif // BASEJET_H
