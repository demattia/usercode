#ifndef SMARTPROPAGATORWITHIPCOMPONENTSRECORD_H
#define SMARTPROPAGATORWITHIPCOMPONENTSRECORD_H

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
// #include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
// #include "Geometry/Records/interface/IdealGeometryRecord.h"
// #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "boost/mpl/vector.hpp"

class SmartPropagatorWithIPComponentsRecord: public edm::eventsetup::DependentRecordImplementation<SmartPropagatorWithIPComponentsRecord,
					     // boost::mpl::vector<IdealMagneticFieldRecord, GlobalTrackingGeometryRecord, TransientTrackRecord> > {};
  boost::mpl::vector<TrackingComponentsRecord, TransientTrackRecord> > {};

#endif // SMARTPROPAGATORWITHIPCOMPONENTSRECORD_H
