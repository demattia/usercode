#ifndef DELTAPHI_H
#define DELTAPHI_H

#include <cmath>

inline double DeltaPhi ( double Phi1, double Phi2 ) {

  return ( 3.141593 - std::fabs( std::fabs( Phi1 - Phi2 ) - 3.141593 ) );

}

inline float DeltaPhi ( float Phi1, float Phi2 ) {

  return ( 3.141593 - std::fabs( std::fabs( Phi1 - Phi2 ) - 3.141593 ) );

}

#endif // DELTAPHI_H
