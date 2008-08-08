#ifndef PARTICLE_H
#define PARTICLE_H

/**
 * Particle class used to store decayed particles.
 * It contains a vector of pointers to the jets from which
 * the particle has been reconstructed.
 * Inherits from BaseParticle.
 * It can be built from any jet objects that have the method:
 * et(), eta(), phi().
 *
 * Author: M. De Mattia
 * Date: 5/7/2008
 */

#include <vector>
#include <cmath>
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

using namespace std;

namespace anaobj {

  template <class T>
  class Particle : public BaseParticle {
  public:
    /// Default constructor, initializes to zero
    Particle();
    /// Constructor receives a pointer to the first jet
    Particle(const T * jet);
    /** Method to add another jet. The quadri-momentum is summed
     * to that of the preceding ones (if any).
     */
    void add(const T * jet);
    /** Method to add another particle. The quadri-momentum is summed
     * to that of the preceding ones (if any).
     */
    void add(const Particle<T> * particle);
    double e() const { return e_; }
    double mass() const { return mass_; }
    vector<const T *> components() const { return decayJets_; }
    vector<const Particle<T> *> componentParticles() const { return decayParticles_; }
    /// Scalar product: corresponds to the scalar product of the spacelike part of the quadrimomentum.
    double operator*( const BaseParticle * second ) const;
  protected:
    double px_;
    double py_;
    double pz_;
    double e_;
    double mass_;
    vector<const T *> decayJets_;
    vector<const Particle<T> *> decayParticles_;
  };

  template <class T>
  Particle<T>::Particle() : BaseParticle( 0., 0., 0. ) {
    px_ = 0.;
    py_ = 0.;
    pz_ = 0.;
    e_ = 0.;
    mass_ = 0.;
  }

  template <class T>
  Particle<T>::Particle( const T * jet ) : BaseParticle( jet->et(), jet->eta(), jet->phi() ) {
    decayJets_.push_back(jet);
    px_ = jet->ex();
    py_ = jet->ey();
    pz_ = jet->ez();
    e_ = jet->e();
    mass_ = sqrt( pow(e_,2) - pow(p(),2) );
  }

  template <class T>
  void Particle<T>::add( const T * jet ) {
    decayJets_.push_back(jet);
    px_ += jet->ex();
    py_ += jet->ey();
    pz_ += jet->ez();
    pt_ = sqrt( pow(px_,2) + pow(py_,2) );
    e_ += jet->e();
    // Evaluate the p here since the p() method evaluates px, py and pz from pt, eta, phi
    double p = sqrt(pow(pt_,2) + pow(pz_,2));
    mass_ = sqrt( pow(e_,2) - pow(p,2) );
    phi_ = atan2(py_,px_);
    eta_ = 0.5*log((p+pz_)/(p-pz_));
  }

  template <class T>
  void Particle<T>::add( const Particle<T> * particle ) {
    // Add all the jet of the particle to this particle's list
    const vector<const T *> & particleComponents = particle->components();
    typename vector<const T *>::const_iterator jetIter = particleComponents.begin();
    for ( ; jetIter != particleComponents.end(); ++jetIter ) {
      decayJets_.push_back(*jetIter);
    }
    // Add this particle to the particles list
    decayParticles_.push_back(particle);
    px_ += particle->px();
    py_ += particle->py();
    pz_ += particle->pz();
    pt_ = sqrt( pow(px_,2) + pow(py_,2) );
    e_ += particle->e();
    // Evaluate the p here since the p() method evaluates px, py and pz from pt, eta, phi
    double p = sqrt(pow(pt_,2) + pow(pz_,2));
    mass_ = sqrt( pow(e_,2) - pow(p,2) );
    phi_ = atan2(py_,px_);
    eta_ = 0.5*log((p+pz_)/(p-pz_));
  }

  template <class T>
  double Particle<T>::operator*( const BaseParticle * second ) const {
    return(px_*second->px() + py_*second->py() + pz_*second->pz());
  }
}

#endif // PARTICLE_H
