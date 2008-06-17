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
    OfflineMEt( const double & MET, const double & CORRL2ET, const double & CORRL3ET, const double & PHI, const double & CORRL2PHI, const double & CORRL3PHI, const double & SUMET, const double & CORRL2SUMET, const double & CORRL3SUMET, const double & METSIG, const double & CORRL2METSIG, const double & CORRL3METSIG ) : BaseMEt( MET, PHI, SUMET ) {
      corrL2et_ = CORRL2ET;
      corrL3et_ = CORRL3ET;
      corrL2phi_ = CORRL2PHI;
      corrL3phi_ = CORRL3PHI;
      corrL2sumEt_ = CORRL2SUMET;
      corrL3sumEt_ = CORRL3SUMET;
      mEtSig_ = METSIG;
      corrL2mEtSig_ = CORRL2METSIG;
      corrL3mEtSig_ = CORRL3METSIG;
//      dPhiMin_ = DPHIMIN;
    }
    double corrL2et() const { return corrL2et_; }
    double corrL3et() const { return corrL3et_; }
    double corrL2phi() const { return corrL2phi_; }
    double corrL3phi() const { return corrL3phi_; }
    double corrL2sumEt() const { return corrL2sumEt_; }
    double corrL3sumEt() const { return corrL3sumEt_; }
    double mEtSig() const { return mEtSig_; }
    double corrL2mEtSig() const { return corrL2mEtSig_; }
    double corrL3mEtSig() const { return corrL3mEtSig_; }
//     double dPhiMin() const { return dPhiMin_; }
    void setCorrL2et( const double & CORRL2ET ) { corrL2et_ = CORRL2ET; }
    void setCorrL3et( const double & CORRL3ET ) { corrL3et_ = CORRL3ET; }
    void setCorrL2phi( const double & CORRL2PHI ) { corrL2phi_ = CORRL2PHI; }
    void setCorrL3phi( const double & CORRL3PHI ) { corrL3phi_ = CORRL3PHI; }
    void setCorrL2sumEt( const double & CORRL2SUMET ) { corrL2sumEt_ = CORRL2SUMET; }
    void setCorrL3sumEt( const double & CORRL3SUMET ) { corrL3sumEt_ = CORRL3SUMET; }
    void setMEtSig( const double & METSIG ) { mEtSig_ = METSIG; }
    void setCorrL2mEtSig( const double & CORRL2METSIG ) { corrL2mEtSig_ = CORRL2METSIG; }
    void setCorrL3mEtSig( const double & CORRL3METSIG ) { corrL3mEtSig_ = CORRL3METSIG; }
//     void setDPhiMin( const double & DPHIMIN ) { dPhiMin_ = DPHIMIN; }

  protected:
    double corrL2et_;
    double corrL3et_;
    double corrL2phi_;
    double corrL3phi_;
    double corrL2sumEt_;
    double corrL3sumEt_;
    double mEtSig_;
    double corrL2mEtSig_;
    double corrL3mEtSig_;
//     double dPhiMin_;
  };

}

#endif // OFFLINEMET_H
