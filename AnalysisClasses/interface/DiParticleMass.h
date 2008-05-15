#ifndef DIPARTICLEMASS_H
#define DIPARTICLEMASS_H

#include <cmath>

inline double DiMass(double E1,  double E2,
		     double px1, double px2,
		     double py1, double py2,
		     double pz1, double pz2){
  double mass = pow(E1+E2,2) - 
    (pow(px1+px2,2)+
     pow(py1+py2,2)+
     pow(pz1+pz2,2));
  
  if (mass > 0. ) mass = TMath::Sqrt(mass);
  else mass = 0.;
  return mass;
}


inline double DiMass ( double pT1, double phi1, double eta1,
		       double pT2, double phi2, double eta2) {

  double theta1 = 2*TMath::ATan(TMath::Exp(-eta1));
  double theta2 = 2*TMath::ATan(TMath::Exp(-eta2));
  double px1 = pT1*TMath::Cos(phi1);
  double px2 = pT1*TMath::Cos(phi2);
  double py1 = pT1*TMath::Sin(phi1);
  double py2 = pT1*TMath::Sin(phi2);
  double pz1 = pT1/TMath::Tan(theta1);
  double pz2 = pT1/TMath::Tan(theta2);

  double energy1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  double energy2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);

  return DiMass(energy1,energy2,
		px1,px2,
		py1,py2,
		pz1,pz2);
}

template<class T1>
inline double DiParticleMass ( const T1 * particle1, 
			       const T1 * particle2 ) {

  double particle1E  = particle1->energy();
  double particle1Px = particle1->px();
  double particle1Py = particle1->py();
  double particle1Pz = particle1->pz();
  double particle2E  = particle2->energy();
  double particle2Px = particle2->px();
  double particle2Py = particle2->py();
  double particle2Pz = particle2->pz();

  return DiMass(particle1E,particle2E,
		particle1Px,particle2Px,
		particle1Py,particle2Py,
		particle1Pz,particle2Pz);
}

#endif // DIPARTICLEMASS_H
