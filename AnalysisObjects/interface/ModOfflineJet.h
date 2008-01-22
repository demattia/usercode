#ifndef MODOFFLINEJET_H
#define MODOFFLINEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"

namespace anaobj {

  /**
   *
   * Used for multiplication. It is the same as the template Jet, but stores
   * also the original et of the jet (origEt) and the et of the associated jet, GenJet for example, (refEt).
   * Inherits from the template parameter class (OfflineJet or BaseJet for example).
   *
   * The template parameter must be a class with a et() method.
   *
   * Author M. De Mattia - 6/1/2008
   *
   */

  class ModOfflineJet : public OfflineJet {
  public:
    /// Constructor from OfflineJet
    ModOfflineJet( const OfflineJet & JET ) : OfflineJet( JET ) {
      origEt_ = JET.et();
      refEt_ = 0.;
      refEta_ = 0.;
    }
    double origEt() const {
      return origEt_;
    }
    double refEt() const {
      return refEt_;
    }
    double refEta() const {
      return refEta_;
    }
    void setOrigEt( const double & ORIGET ) {
      origEt_ = ORIGET;
    }
    void setRefEt( const double & REFET ) {
      refEt_ = REFET;
    }
    void setRefEta( const double & REFETA ) {
      refEta_ = REFETA;
    }
    void setEt( const double & ET ) {
      et_ = ET;
    }

  protected:
    double origEt_;
    double refEt_;
    double refEta_;
  };
}

#endif // MODOFFLINEJET_H
