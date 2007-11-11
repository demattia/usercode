#ifndef SIMPLEELECTRON_H
#define SIMPLEELECTRON_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

#include <cmath>
#include <vector>

/**
 *
 * Stores an offline electron. Inlcudes Et, had/em fraction and
 * isolation.
 *
 * Inherits from BaseParticle.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class SimpleElectron : public BaseParticle {
  public:
    SimpleElectron( const double & PT, const double & ETA, const double & PHI,
                    const double & ET, const double & HADOVEREM, const double & ISOVAL ) : BaseParticle( PT, ETA, PHI ) {
      et_ = ET;
      hadOverEm_ = HADOVEREM;
      isoVal_ = ISOVAL;
    }
    /// Default constructor, only needed for classes.h
    SimpleElectron() : BaseParticle( 0., 0., 0. ) {
      et_ = 0.;
      hadOverEm_ = 0.;
      isoVal_ = 0.;
    }
    double et() const { return et_; }
    double hadOverEm() const { return hadOverEm_; }
    double isoVal() const { return isoVal_; }
    void setEt( const double & ET ) { et_ = ET; }
    void setHadOverEm( const double & HADOVEREM ) { hadOverEm_ = HADOVEREM; }
    void setIsoVal( const double & ISOVAL ) { isoVal_ = ISOVAL; }
  protected:
    double et_;
    double hadOverEm_;
    double isoVal_;
  };

  typedef std::vector<SimpleElectron> SimpleElectronCollection;

}

#endif // GLOBALMUON_H
