#ifndef BASEALL_H
#define BASEALL_H

/**
 *
 * Base class for all jets and particles to be stored.
 * It has only eta and phi.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class BaseAll {
  public:
    BaseAll( const double & ETA, const double & PHI ) {
      eta_ = ETA;
      phi_ = PHI;
    }
    double eta() const { return eta_; }
    double phi() const { return phi_; }
    void setEta( const double & ETA ) { eta_ = ETA; }
    void setPhi( const double & PHI ) { phi_ = PHI; }
  protected:
    double eta_;
    double phi_;
  };

}

#endif //BASEALL_H
