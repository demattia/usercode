#ifndef FlatRandomD0GunProducer_H
#define FlatRandomD0GunProducer_H

/** \class FlatRandomD0GunProducer
 *
 * Generates single particle gun in HepMC format
 * Julia Yarba 12/2005 
 ***************************************/

#include "HarderAnalysisTools/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm
{
  
  class FlatRandomD0GunProducer : public BaseFlatGunProducer
  {
  
  public:
    FlatRandomD0GunProducer(const ParameterSet & pset);
    virtual ~FlatRandomD0GunProducer();
   
    virtual void produce(Event & e, const EventSetup& es);

  private:
    
    // data members
    
    double            fMinD0   ;
    double            fMaxD0   ;
    double            fMinPt   ;
    double            fMaxPt   ;

  };
} 

#endif
