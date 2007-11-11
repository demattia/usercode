#ifndef SIMPLETAU_H
#define SIMPLETAU_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"

#include <cmath>
#include <vector>

/**
 *
 * Used for taus.
 * Inherits from BaseJet.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class SimpleTau : public BaseJet {
  public:
    SimpleTau( const double & ET, const double & ETA, const double & PHI, const double & TKMASS, const int TKNUM ) : BaseJet( ET, ETA, PHI ) {
      tkMass_ = TKMASS;
      tkNum_ = TKNUM;
    }
    // Default constructor, only needed for classes.h
    SimpleTau() : BaseJet( 0., 0., 0. ) {
      tkMass_ = 0.;
      tkNum_ = 0;
    }
    double tkMass() const { return tkMass_; }
    int tkNum() const { return tkNum_; }
    void setTkMass( const double & TKMASS ) { tkMass_ = TKMASS; }
    void setTkNum( const int TKNUM ) { tkNum_ = TKNUM; }
  protected:
    double tkMass_;
    int tkNum_;
  };

  typedef std::vector<SimpleTau> SimpleTauCollection;

}

#endif // SIMPLETAU_H
