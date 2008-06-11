#ifndef OFFLINEMET_H
#define OFFLINEMET_H

#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"

namespace anaobj {

  /**
   *
   * MEt class used to store offline MEt, MEt phi,
   * sumEt, MEtSignificance and minimum DPhi between
   * MEt and the closest offline Jet.
   * Inherits from BaseMEt.
   *
   * Author M. De Mattia - 9/11/2007
   *
   * Added corrL2et and corrL3et for the MEt corrections coming from L2 and L3 corrected jets - 11/6/2008
   *
   */

  class OfflineMEt : public BaseMEt
  {
  public:
    /// Default empty constructor, needed to make it become a product
    OfflineMEt() {}
//    OfflineMEt( const double & MET, const double & PHI, const double & SUMET, const double & METSIG, const double & DPHIMIN ) : BaseMEt( MET, PHI, SUMET ) {
    OfflineMEt( const double & MET, const double & CORRL2ET, const double & CORRL3ET, const double & PHI, const double & SUMET, const double & METSIG ) : BaseMEt( MET, PHI, SUMET ) {
      corrL2et_ = CORRL2ET;
      corrL3et_ = CORRL3ET;
      mEtSig_ = METSIG;
//      dPhiMin_ = DPHIMIN;
    }
    double mEtSig() const { return mEtSig_; }
    double corrL2et() const { return corrL2et_; }
    double corrL3et() const { return corrL3et_; }
//     double dPhiMin() const { return dPhiMin_; }
    void setMEtSig( const double & METSIG ) { mEtSig_ = METSIG; }
    void setCorrL2et( const double & CORRL2ET ) { corrL2et_ = CORRL2ET; }
    void setCorrL3et( const double & CORRL3ET ) { corrL3et_ = CORRL3ET; }
//     void setDPhiMin( const double & DPHIMIN ) { dPhiMin_ = DPHIMIN; }

  protected:
    double mEtSig_;
    double corrL2et_;
    double corrL3et_;
//     double dPhiMin_;
  };

}

#endif // OFFLINEMET_H
