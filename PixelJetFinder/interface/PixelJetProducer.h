//
// Package:     AnalysisExample/PixelJetFinder
// Class:       PixelJetProducer
//
//
// Description: Produces PixelJets from PixelTracks
//
// Original     Authors: M. De Mattia
// Created:     6/8/2007
//

#ifndef AnalysisExamples_PixelJetFinder_h
#define AnalysisExamples_PixelJetFinder_h

//#include <map>

// Needed to use the DEFINE_FW_MODULE directive at the end of the cc file
// and register this class as a framework module
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
// #include "Geometry/Vector/interface/GlobalVector.h"
// #include "Geometry/Vector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
/* #include "DataFormats/TrackReco/interface/TrackExtra.h" */
/* #include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h" */
/* #include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h" */
/* #include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h" */
/* #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h" */
/* #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" */
/* #include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h" */
/* #include "TrackingTools/KalmanUpdators/interface/KFUpdator.h" */
/* #include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h" */
/* #include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h" */
/* #include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h" */
/* #include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"  */
/* #include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h" */
/* #include "FWCore/MessageLogger/interface/MessageLogger.h" */
/* #include "DataFormats/DetId/interface/DetId.h" */
/* #include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h" */
/* #include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h" */
/* #include "AnalysisDataFormats/SiStripClusterInfo/interface/SiStripClusterInfo.h" */

#include "AnalysisExamples/PixelJet/interface/PixelJet.h"

// Added for the trackassociator
/* #include "SimTracker/TrackAssociation/test/testTrackAssociator.h" */
/* #include "SimTracker/Records/interface/TrackAssociatorRecord.h" */
/* #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h" */
/* #include "DataFormats/TrackReco/interface/TrackFwd.h" */
/* #include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" */

//#include <TMath.h>


class PixelJetProducer : public edm::EDProducer
{
 public:
  
  explicit PixelJetProducer(const edm::ParameterSet& conf);
  
  virtual ~PixelJetProducer();
  
  virtual void beginJob(const edm::EventSetup& c);
  
  virtual void endJob(); 
  
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
 private:

  edm::ParameterSet conf_;
  std::string filename_;
  std::string pixeljet_;
  double eta_cut_;
  double ConeR_cut_;
  int NumTk_cut_;
  float Eta_cut_Eff_;
  unsigned int eventcounter_;
};

#endif // AnalysisExamples_PixelJetFinder_h
