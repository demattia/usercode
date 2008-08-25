#ifndef SUMJETS_H
#define SUMJETS_H

/**
 * Functor to evaluate the sum of two jets.
 * Takes the Px, Py, Pz ed E of the jets.
 * Evaluates the quantities of the resulting jet as:
 * Px = Px1+Px2, Py = Py1+Py2, Pz = Pz1+Pz2
 * Phi = atan(Py/Px)
 * Pt = sqrt(Px^2+Py^2)
 * P = sqrt(Pt^2+Pz^2)
 * Eta = 0.5*log((P+Pz)/(P-Pz))
 * assumes a constructor of the form: Jet( const double & ET, const double & ETA, const double & PHI )
 * ATTENTION:
 * Since the classes inheriting from BaseJet evaluate the E from the sqrt(Ex^2+Ey^2+Ez^2), it coincides
 * with P and cannot be used to evaluate the mass from the difference M = sqrt(E^2-P^2).
 * This means that E != E1+E2 and it cannot be assigned with the current form of the class.
 * It should be changed in the next production, but for now all the masses evaluated in this way would be = 0.
 * To evaluate the mass simply take the E from the resulting jet and the E1 and E2 from the added jets and
 * perform the computation by hand.
 *
 * It returns an object of type BaseJet by value.
 *
 * Author: M. De Mattia
 * Date: 1/7/2008
 */

#include <cmath>
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"

template <class T>
class SumJets {
public:
  BaseJet operator()(const T * jet1, const T * jet2) const;
};

template <class T>
BaseJet SumJets<T>::operator()(const T * jet1, const T * jet2) const {
  double px = jet1->ex()+jet2->ex();
  double py = jet1->ey()+jet2->ey();
  double pz = jet1->ez()+jet2->ez();
  double phi = atan2(py,px);
  double pt = sqrt(pow(px,2)+pow(py,2));
  double p = sqrt(pow(pt,2)+pow(pz,2));
  double eta = 0.5*log((p+pz)/(p-pz));

  return( BaseJet(pt, eta, phi) );
}



#endif // SUMJETS_H
