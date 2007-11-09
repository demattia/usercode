#ifndef SIMPLEJET_H
#define SIMPLEJET_H

#include <cmath>

// Add a simple Jet class to allow merging of CenJets and TauJets to be used with the Associator.

class SimpleJet {
 public:
  SimpleJet( double PT, double ETA, double PHI ) {
    pt_ = PT;
    eta_ = ETA;
    phi_ = PHI;
  }
  double pt() const {
    return pt_;
  }
  double eta() const {
    return eta_;
  }
  double phi() const {
    return phi_;
  }
  float px() const {
    return pt_*std::cos(phi_);
  }
  float py() const {
    return pt_*std::sin(phi_);
  }
  float pz() const {
    return pt_*(1-std::exp(-2*eta_))/(2*std::exp(-eta_));
  }
  float E() const {
    return std::sqrt(std::pow(px(),2)+std::pow(py(),2)+std::pow(pz(),2));
  }
  // To sort the SimpleJets
  bool operator< ( const SimpleJet& b ) const {
    return pt() < b.pt();
  }
 private:
  double pt_;
  double eta_;
  double phi_;
};


#endif // SIMPLEJET_H
