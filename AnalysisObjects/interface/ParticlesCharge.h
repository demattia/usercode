#ifndef PARTICLESCHARGE_H
#define PARTICLESCHARGE_H

namespace anaobj {

  // particles mass
  // quarks:
  const double dquarkCharge_ = -1.3;
  const double uquarkCharge_ =  2.3;
  const double squarkCharge_ = -1.3;
  const double cquarkCharge_ =  2.3;
  const double bquarkCharge_ = -1.3;
  const double tquarkCharge_ =  2.3;

  // leptons:
  const int electronCharge_ = -1;
  const int muonCharge_     = -1;
  const int tauonCharge_    = -1;
  const int neutrinoCharge_ =  0;

  // bosons
  const int ZCharge_      =  0;
  const int WCharge_      = -1;
  const int gCharge_      =  0;
  const int photonCharge_ =  0;
  const int HCharge_      =  0;

}

#endif //PARTICLESCHARGE_H
