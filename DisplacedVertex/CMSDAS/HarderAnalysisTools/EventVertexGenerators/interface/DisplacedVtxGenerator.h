#ifndef HarderAnalysisTools_DisplacedVtxGenerator_H
#define HarderAnalysisTools_DisplacedVtxGenerator_H

/**
 * Generate event vertices according to a Flat distribution in cylindrical coordinates. 
 * Attention: All values are assumed to be cm!
 *
 * $Id: DisplacedVtxGenerator.h,v 1.1 2009/12/11 17:05:23 harder Exp $
 */

#include "HarderAnalysisTools/EventVertexGenerators/interface/BaseEvtVtxGenerator.h"

namespace CLHEP {
   class RandFlat;
}

class DisplacedVtxGenerator : public BaseEvtVtxGenerator 
{
public:
  DisplacedVtxGenerator(const edm::ParameterSet & p);
  virtual ~DisplacedVtxGenerator();

  /// return a new event vertex
  virtual HepMC::FourVector* newVertex() ;

  virtual TMatrixD* GetInvLorentzBoost() {
	  return 0;
  }

    
  /// set min in R in cm
  void minR(double m=0.0);
  /// set min in Phi in radians
  void minPhi(double m=0.0);
  /// set min in Z in cm
  void minZ(double m=0.0);

  /// set max in R in cm
  void maxR(double m=0);
  /// set max in Phi in radians
  void maxPhi(double m=0);
  /// set max in Z in cm
  void maxZ(double m=0);
  
private:
  /** Copy constructor */
  DisplacedVtxGenerator(const DisplacedVtxGenerator &p);
  /** Copy assignment operator */
  DisplacedVtxGenerator&  operator = (const DisplacedVtxGenerator & rhs );
private:
  double fMinR, fMinPhi, fMinZ;
  double fMaxR, fMaxPhi, fMaxZ;
  CLHEP::RandFlat*  fRandom ;
  double fTimeOffset;
};

#endif
