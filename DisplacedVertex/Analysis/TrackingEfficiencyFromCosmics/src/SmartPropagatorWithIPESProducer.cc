/** \class SmartPropagatorWithIPESProducer
 *  ES producer needed to put the SmartPropagatorWithIP inside the EventSetup
 *
 *  $Date: 2011/10/27 17:30:00 $
 *  $Revision: 1.0 $
 *  \author M. De Mattia - <m.de.mattia@cern.ch>
 */

#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIPESProducer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

using namespace edm;
using namespace std;

SmartPropagatorWithIPESProducer::SmartPropagatorWithIPESProducer(const ParameterSet& parameterSet)
{
  string myname = parameterSet.getParameter<string>("ComponentName");

  string propDir = parameterSet.getParameter<string>("PropagationDirection");

  if (propDir == "oppositeToMomentum") thePropagationDirection = oppositeToMomentum;
  else if (propDir == "alongMomentum") thePropagationDirection = alongMomentum;
  else if (propDir == "anyDirection") thePropagationDirection = anyDirection;
  else
    throw cms::Exception("SmartPropagatorWithIPESProducer")
      << "Wrong fit direction chosen in SmartPropagatorWithIPESProducer";


  theEpsilon = parameterSet.getParameter<double>("Epsilon");

  theTrackerPropagatorName = parameterSet.getParameter<string>("TrackerPropagator");
  theMuonPropagatorName = parameterSet.getParameter<string>("MuonPropagator");

  setWhatProduced(this,myname);
}

SmartPropagatorWithIPESProducer::~SmartPropagatorWithIPESProducer() {}

boost::shared_ptr<Propagator>
SmartPropagatorWithIPESProducer::produce(const SmartPropagatorWithIPComponentsRecord& iRecord){

  std::cout << "Building the SmartPropagatorWithIP" << std::endl;

  ESHandle<MagneticField> magField;
  iRecord.getRecord<TrackingComponentsRecord>().getRecord<IdealMagneticFieldRecord>().get(magField);

  ESHandle<Propagator> trackerPropagator;
  iRecord.getRecord<TrackingComponentsRecord>().get(theTrackerPropagatorName,trackerPropagator);

  ESHandle<Propagator> muonPropagator;
  iRecord.getRecord<TrackingComponentsRecord>().get(theMuonPropagatorName,muonPropagator);

  ESHandle<TransientTrackBuilder> theBuilder;
  iRecord.getRecord<TransientTrackRecord>().get("TransientTrackBuilder",theBuilder);

  thePropagator  = boost::shared_ptr<Propagator>(new SmartPropagatorWithIP(*trackerPropagator, *muonPropagator, &*magField,
                                                                           &*theBuilder,
                                                                           thePropagationDirection, theEpsilon));
  return thePropagator;
}
