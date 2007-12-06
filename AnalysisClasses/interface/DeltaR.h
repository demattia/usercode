#ifndef DELTAR_H
#define DELTAR_H

#include <cmath>
#include "AnalysisExamples/AnalysisClasses/interface/DeltaPhi.h"

inline double DeltaR ( double Eta1,
		       double Phi1, 
		       double Eta2, 
		       double Phi2) {

  //  double PI_ = 3.141593;

  double delPhi;
  double delEta;

  //  delPhi = PI_ - std::fabs( std::fabs( Phi1 - Phi2 ) - PI_ ) );
  delPhi = DeltaPhi( Phi1, Phi2 );
  delEta = Eta1 - Eta2;

  return ( std::sqrt(delPhi*delPhi + delEta*delEta ) );


}

#endif // DELTAR_H
