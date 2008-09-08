#ifndef DELTAPHI_H
#define DELTAPHI_H

#include "TMath.h"
#include <cmath>

namespace {
  double pi = TMath::Pi();
}

template<class T1>
inline double DeltaPhi ( T1* particle1, T1* particle2 ) {

  double particle1Phi = particle1->phi();
  double particle2Phi = particle2->phi();
  return ( pi - std::fabs( std::fabs( particle1Phi - particle2Phi ) - pi ) );

}

inline double DeltaPhi ( double Phi1, double Phi2 ) {

  return ( pi - std::fabs( std::fabs( Phi1 - Phi2 ) - pi ) );

}

inline float DeltaPhi ( float Phi1, float Phi2 ) {

  return ( pi - std::fabs( std::fabs( Phi1 - Phi2 ) - pi ) );

}

#endif // DELTAPHI_H
