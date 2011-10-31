#ifndef Analysis_TrackingEfficiencyFromCosmics_SmartPropagatorWithIPESProducer_H
#define Analysis_TrackingEfficiencyFromCosmics_SmartPropagatorWithIPESProducer_H

/** \class SmartPropagatorWithIPESProducer
 *  ES producer needed to put the SmartPropagatorWithIP inside the EventSetup
 *
 *  $Date: 2011/10/27 17:30:00 $
 *  $Revision: 1.0 $
 *  \author M. De Mattia - <m.de.mattia@cern.ch>
 */

#include "FWCore/Framework/interface/ESProducer.h"

#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include <boost/shared_ptr.hpp>

#include "Analysis/Records/interface/SmartPropagatorWithIPComponentsRecord.h"

namespace edm {class ParameterSet;}

class TrackingComponentsRecord;

class  SmartPropagatorWithIPESProducer: public edm::ESProducer{

 public:

  /// Constructor
  SmartPropagatorWithIPESProducer(const edm::ParameterSet &);

  /// Destructor
  virtual ~SmartPropagatorWithIPESProducer();

  // Operations
  boost::shared_ptr<Propagator> produce(const SmartPropagatorWithIPComponentsRecord &);
  // boost::shared_ptr<Propagator> produce(const TrackingComponentsRecord &);

 private:
  boost::shared_ptr<Propagator> thePropagator;
  PropagationDirection thePropagationDirection;
  std::string theTrackerPropagatorName;
  std::string theMuonPropagatorName;
  double theEpsilon;
};

#endif
