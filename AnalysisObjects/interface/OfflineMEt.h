#ifndef OFFLINEMET_H
#define OFFLINEMET_H

#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"

/**
 *
 * MEt class used to store offline MEt, MEt phi,
 * sumEt, MEtSignificance and minimum DPhi between
 * MEt and the closest offline Jet.
 * Inherits from BaseMEt.
 *
 * Author M. De Mattia - 9/11/2007
 *
 */

namespace anaobj {

  class OfflineMEt : public BaseMEt
  {
  public:
    /// Default empty constructor, needed to make it become a product
    OfflineMEt() {}
    OfflineMEt( const double & MET, const double & PHI, const double & SUMET, const double & METSIGNIFICANCE, const double & DPHIMIN ) : BaseMEt( MET, PHI, SUMET ) {
      mEtSignificance_ = METSIGNIFICANCE;
      dPhiMin_ = DPHIMIN;
    }
    double mEtSignificance() const { return mEtSignificance_; }
    double dPhiMin() const { return dPhiMin_; }
    void setMEtSignificance( const double & METSIGNIFICANCE ) { mEtSignificance_ = METSIGNIFICANCE; }
    void setDPhiMin( const double & DPHIMIN ) { dPhiMin_ = DPHIMIN; }

  protected:
    double mEtSignificance_;
    double dPhiMin_;
  };

}

#endif // OFFLINEMET_H
