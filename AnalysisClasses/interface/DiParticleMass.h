#ifndef DIPARTICLEMASS_H
#define DIPARTICLEMASS_H

#include <cmath>

template<class T1>
inline double DiParticleMass ( const T1 * particle1, const T1 * particle2 ) {

  double particle1E  = particle1->energy();
  double particle1Px = particle1->px();
  double particle1Py = particle1->py();
  double particle1Pz = particle1->pz();
  double particle2E  = particle2->energy();
  double particle2Px = particle2->px();
  double particle2Py = particle2->py();
  double particle2Pz = particle2->pz();

  double mass = pow(particle1E+particle2E,2) - 
    (pow(particle1Px+particle2Px,2)+
     pow(particle1Py+particle2Py,2)+
     pow(particle1Pz+particle2Pz,2));

  if (mass > 0. ) mass = TMath::Sqrt(mass);
  else mass = 0.;
  return mass;
//  return dimass(particle1E,particle2E,
//		particle1Px,particle2Px,
//		particle1Py,particle2Py,
//		particle1Pz,particle2Pz);
}

//template<class T1>
////inline double DiParticleMass ( const T1 & particle1, const T1 & particle2 ) {
//inline double DiParticleMass ( T1 & particle1, T1 & particle2 ) {
//
//  double particle1E  = particle1.energy();
//  double particle1Px = particle1.px();
//  double particle1Py = particle1.py();
//  double particle1Pz = particle1.pz();
//  double particle2E  = particle2.energy();
//  double particle2Px = particle2.px();
//  double particle2Py = particle2.py();
//  double particle2Pz = particle2.pz();
//
//  double mass = pow(particle1E+particle2E,2) - 
//    (pow(particle1Px+particle2Px,2)+
//     pow(particle1Py+particle2Py,2)+
//     pow(particle1Pz+particle2Pz,2));
//  
//  if (mass > 0. ) mass = TMath::Sqrt(mass);
//  else mass = 0.;
//  return mass;
////  return dimass(particle1E,particle2E,
////		particle1Px,particle2Px,
////		particle1Py,particle2Py,
////		particle1Pz,particle2Pz);
//}

#endif // DIPARTICLEMASS_H
