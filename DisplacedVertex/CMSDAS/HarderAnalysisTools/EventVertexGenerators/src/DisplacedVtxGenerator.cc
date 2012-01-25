
// $Id: DisplacedVtxGenerator.cc,v 1.1 2009/12/11 17:05:24 harder Exp $

#include "HarderAnalysisTools/EventVertexGenerators/interface/DisplacedVtxGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
//#include "CLHEP/Vector/ThreeVector.h"
#include "HepMC/SimpleVector.h"

DisplacedVtxGenerator::DisplacedVtxGenerator(const edm::ParameterSet& p )
: BaseEvtVtxGenerator(p)
{ 
  
  fRandom = new CLHEP::RandFlat(getEngine()) ;
  
  fMinR = p.getParameter<double>("MinR")*cm;
  fMinPhi = p.getParameter<double>("MinPhi");
  fMinZ = p.getParameter<double>("MinZ")*cm;
  fMaxR = p.getParameter<double>("MaxR")*cm;
  fMaxPhi = p.getParameter<double>("MaxPhi");
  fMaxZ = p.getParameter<double>("MaxZ")*cm;     
  fTimeOffset = p.getParameter<double>("TimeOffset")*ns*c_light;
  
  if (fMinR > fMaxR) {
    throw cms::Exception("Configuration")
      << "Error in DisplacedVtxGenerator: "
      << "MinR is greater than MaxR";
  }
  if (fMinPhi > fMaxPhi) {
    throw cms::Exception("Configuration")
      << "Error in DisplacedVtxGenerator: "
      << "MinPhi is greater than MaxPhi";
  }
  if (fMinZ > fMaxZ) {
    throw cms::Exception("Configuration")
      << "Error in DisplacedVtxGenerator: "
      << "MinZ is greater than MaxZ";
  }
}

DisplacedVtxGenerator::~DisplacedVtxGenerator() 
{
  delete fRandom; 
}


//Hep3Vector * DisplacedVtxGenerator::newVertex() {
HepMC::FourVector* DisplacedVtxGenerator::newVertex() {
  double aR,aPhi,aZ;
  aR = fRandom->fire(fMinR,fMaxR) ;
  aPhi = fRandom->fire(fMinPhi,fMaxPhi) ;
  aZ = fRandom->fire(fMinZ,fMaxZ) ;

  double aX,aY;
  aX = aR*std::cos(aPhi);
  aY = aR*std::sin(aPhi);

  //if (fVertex == 0) fVertex = new CLHEP::Hep3Vector;
  //fVertex->set(aX,aY,aZ);
  if ( fVertex == 0 ) fVertex = new HepMC::FourVector() ;
  fVertex->set(aX,aY,aZ,fTimeOffset);

  return fVertex;
}

void DisplacedVtxGenerator::minR(double min) 
{
  fMinR = min;
}

void DisplacedVtxGenerator::minPhi(double min) 
{
  fMinPhi = min;
}

void DisplacedVtxGenerator::minZ(double min) 
{
  fMinZ = min;
}

void DisplacedVtxGenerator::maxR(double max) 
{
  fMaxR = max;
}

void DisplacedVtxGenerator::maxPhi(double max) 
{
  fMaxPhi = max;
}

void DisplacedVtxGenerator::maxZ(double max) 
{
  fMaxZ = max;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DisplacedVtxGenerator);
