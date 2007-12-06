#ifndef SIMPLETAU_H
#define SIMPLETAU_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Used for taus.
   * Inherits from BaseJet.
   *
   * Author M. De Mattia - 8/11/2007
   *
   */

  class SimpleTau : public BaseJet {
  public:
    SimpleTau( const double & ET, const double & ETA, const double & PHI, const double & TKMASS, const int TKNUM,
               const double & ISOTKMASS, const int ISOTKNUM ) : BaseJet( ET, ETA, PHI ) {
      tkMass_ = TKMASS;
      tkNum_ = TKNUM;
      isoTkMass_ = ISOTKMASS;
      isoTkNum_ = ISOTKNUM;
    }
    // Default constructor, only needed for classes.h
    SimpleTau() : BaseJet( 0., 0., 0. ) {
      tkMass_ = 0.;
      tkNum_ = 0;
      isoTkMass_ = 0.;
      isoTkNum_ = 0;
    }
    double tkMass() const { return tkMass_; }
    int tkNum() const { return tkNum_; }
    double isoTkMass() const { return isoTkMass_; }
    int isoTkNum() const { return isoTkNum_; }
    void setTkMass( const double & TKMASS ) { tkMass_ = TKMASS; }
    void setTkNum( const int TKNUM ) { tkNum_ = TKNUM; }
    void setIsoTkMass( const double & ISOTKMASS ) { isoTkMass_ = ISOTKMASS; }
    void setIsoTkNum( const int ISOTKNUM ) { isoTkNum_ = ISOTKNUM; }
  protected:
    double tkMass_;
    int tkNum_;
    double isoTkMass_;
    int isoTkNum_;
  };

  typedef std::vector<SimpleTau> SimpleTauCollection;

}

#endif // SIMPLETAU_H
