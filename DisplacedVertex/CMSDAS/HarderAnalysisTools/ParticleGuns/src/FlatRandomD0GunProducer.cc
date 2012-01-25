/*
 *  $Date: 2010/03/02 13:46:35 $
 *  $Revision: 1.1 $
 *  \author Julia Yarba
 */

#include <ostream>

#include "HarderAnalysisTools/ParticleGuns/interface/FlatRandomD0GunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


using namespace edm;
using namespace std;

FlatRandomD0GunProducer::FlatRandomD0GunProducer(const ParameterSet& pset) : 
   BaseFlatGunProducer(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
     pset.getParameter<ParameterSet>("PGunParameters") ;
   
   fMinD0 = pgun_params.getParameter<double>("MinD0");
   fMaxD0 = pgun_params.getParameter<double>("MaxD0");
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");
  
   produces<HepMCProduct>();
   produces<GenEventInfoProduct>();
}

FlatRandomD0GunProducer::~FlatRandomD0GunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

void FlatRandomD0GunProducer::produce(Event &e, const EventSetup& es) 
{

   if ( fVerbosity > 0 )
   {
      cout << " FlatRandomD0GunProducer : Begin New Event Generation" << endl ; 
   }
   // event loop (well, another step in it...)
          
   // no need to clean up GenEvent memory - done in HepMCProduct
   // 
   
   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
   
   // now actualy, cook up the event from PDGTable and gun parameters
   //
   // 1st, primary vertex
   //
   double d0     = fRandomGenerator->fire(fMinD0, fMaxD0) ;
   double phi    = fRandomGenerator->fire(fMinPhi, fMaxPhi) ;
   HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(-d0*sin(phi),d0*cos(phi),0.));

   // loop over particles
   //
   int barcode = 1 ;
   for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
   {

       double pt     = fRandomGenerator->fire(fMinPt, fMaxPt) ;
       double eta    = fRandomGenerator->fire(fMinEta, fMaxEta) ;
       int PartID = fPartIDs[ip] ;
       const HepPDT::ParticleData* 
          PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
       double mass   = PData->mass().value() ;
       double theta  = 2.*atan(exp(-eta)) ;
       double mom    = pt/sin(theta) ;
       double px     = pt*cos(phi) ;
       double py     = pt*sin(phi) ;
       double pz     = mom*cos(theta) ;
       double energy2= mom*mom + mass*mass ;
       double energy = sqrt(energy2) ; 
       HepMC::FourVector p(px,py,pz,energy) ;
       HepMC::GenParticle* Part = 
           new HepMC::GenParticle(p,PartID,1);
       Part->suggest_barcode( barcode ) ;
       barcode++ ;
       Vtx->add_particle_out(Part);

       if ( fAddAntiParticle )
       {
          HepMC::FourVector ap(-px,-py,-pz,energy) ;
	  int APartID = -PartID ;
	  if ( PartID == 22 || PartID == 23 )
	  {
	     APartID = PartID ;
	  }	  
	  HepMC::GenParticle* APart =
	     new HepMC::GenParticle(ap,APartID,1);
	  APart->suggest_barcode( barcode ) ;
	  barcode++ ;
	  Vtx->add_particle_out(APart) ;
       }

   }

   fEvt->add_vertex(Vtx) ;
   fEvt->set_event_number(e.id().event()) ;
   fEvt->set_signal_process_id(20) ; 
        
   if ( fVerbosity > 0 )
   {
      fEvt->print() ;  
   }

   auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
   BProduct->addHepMCData( fEvt );
   e.put(BProduct);

   auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
   e.put(genEventInfo);
    
   if ( fVerbosity > 0 )
   {
      // for testing purpose only
      // fEvt->print() ; // prints empty info after it's made into edm::Event
      cout << " FlatRandomD0GunProducer : Event Generation Done " << endl;
   }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FlatRandomD0GunProducer);
