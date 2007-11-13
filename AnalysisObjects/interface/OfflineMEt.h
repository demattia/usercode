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
    OfflineMEt( const double & MET, const double & PHI, const double & SUMET, const double & METSIG, const double & DPHIMIN ) : BaseMEt( MET, PHI, SUMET ) {
      mEtSig_ = METSIG;
      dPhiMin_ = DPHIMIN;
    }
    double mEtSig() const { return mEtSig_; }
    double dPhiMin() const { return dPhiMin_; }
    void setMEtSig( const double & METSIG ) { mEtSig_ = METSIG; }
    void setDPhiMin( const double & DPHIMIN ) { dPhiMin_ = DPHIMIN; }

  protected:
    double mEtSig_;
    double dPhiMin_;
  };

}

#endif // OFFLINEMET_H
