#ifndef BASEMET_H
#define BASEMET_H

/**
 *
 * Simple MEt class, used to store L1 MEt and SumEt and
 * as a base class for the offline MEt.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class BaseMEt {
  public:
    /// Default empty constructor, needed to make it become a product
    BaseMEt() {};
    /// Receives MEt, MEt phi and SumEt
    BaseMEt( const double & MET, const double & PHI, const double & SUMET ) {
      mEt_ = MET;
      phi_ = PHI;
      sumEt_ = SUMET;
    }
    double MEt() const { return mEt_; }
    double phi() const { return phi_; }
    double sumEt() const { return sumEt_; }
    void setMEt( const double & MET ) { mEt_ = MET; }
    void setPhi( const double & PHI ) { phi_ = PHI; }
    void setSumEt( const double & SUMET ) { sumEt_ = SUMET; }
  protected:
    double mEt_;
    double phi_;
    double sumEt_;
  };

}

#endif //BASEMET_H
