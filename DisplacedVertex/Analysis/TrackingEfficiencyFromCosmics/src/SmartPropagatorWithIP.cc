/* This Class Header */
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"

/* Collaborating Class Header */
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"

#include <DataFormats/GeometrySurface/interface/BoundCylinder.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

/* Base Class Headers */

/* C++ Headers */

/* ====================================================================== */

/* Constructor */
SmartPropagatorWithIP::SmartPropagatorWithIP(Propagator* aTkProp, Propagator* aGenProp, const MagneticField* field,
                                             const TransientTrackBuilder * theBuilder,
                                             PropagationDirection dir, float epsilon) :
  Propagator(dir), theTkProp_(aTkProp->clone()), theGenProp_(aGenProp->clone()), theField_(field), epsilon_(epsilon),
  // transverseExtrapolator_(new TransverseImpactPointExtrapolator(theField_)),
  // analyticalExtrapolator_(new AnalyticalImpactPointExtrapolator(theField_)),
  theBuilder_(theBuilder)
{
  transverseExtrapolator_ = new TransverseImpactPointExtrapolator(theField_);
  analyticalExtrapolator_ = new AnalyticalImpactPointExtrapolator(theField_);
  // Initialize the nullCovariance matrix
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<5; ++j ) {
      nullCovariance_(i,j) = 0;
    }
  }
}

SmartPropagatorWithIP::SmartPropagatorWithIP(const Propagator& aTkProp, const Propagator& aGenProp,const MagneticField* field,
                                             const TransientTrackBuilder * theBuilder,
                                             PropagationDirection dir, float epsilon) :
  Propagator(dir), theTkProp_(aTkProp.clone()), theGenProp_(aGenProp.clone()), theField_(field), epsilon_(epsilon),
  theBuilder_(theBuilder)
{
  transverseExtrapolator_ = new TransverseImpactPointExtrapolator(theField_);
  analyticalExtrapolator_ = new AnalyticalImpactPointExtrapolator(theField_);
  // Initialize the nullCovariance matrix
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<5; ++j ) {
      nullCovariance_(i,j) = 0;
    }
  }
}

SmartPropagatorWithIP::SmartPropagatorWithIP(const SmartPropagatorWithIP& aProp) :
  Propagator(aProp.propagationDirection()), theTkProp_(0), theGenProp_(0)
{
  if (aProp.theTkProp_)
    theTkProp_=aProp.getTkPropagator()->clone();
  if (aProp.theGenProp_)
    theTkProp_=aProp.getGenPropagator()->clone();

  epsilon_ = aProp.epsilon();
  transverseExtrapolator_ = aProp.transverseExtrapolator();
  analyticalExtrapolator_ = aProp.analyticalExtrapolator();
  theBuilder_ = aProp.theBuilder();
  // Initialize the nullCovariance matrix
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<5; ++j ) {
      nullCovariance_(i,j) = 0;
    }
  }
}

/* Destructor */
SmartPropagatorWithIP::~SmartPropagatorWithIP()
{
  delete theTkProp_;
  delete theGenProp_;
}

///* Operations */
TrajectoryStateOnSurface SmartPropagatorWithIP::propagate(const FreeTrajectoryState& fts,
                                                          const Surface& surface) const
{
  return Propagator::propagate( fts, surface);
}

TrajectoryStateOnSurface SmartPropagatorWithIP::propagate(const FreeTrajectoryState& fts,
                                                          const Plane& plane) const
{
  if (insideTkVol(fts) && insideTkVol(plane)) {
    return getTkPropagator()->propagate(fts, plane);
  } else {
    return getGenPropagator()->propagate(fts, plane);
  }
}

TrajectoryStateOnSurface SmartPropagatorWithIP::propagate(const FreeTrajectoryState& fts,
                                                          const Cylinder& cylinder) const
{
  if (insideTkVol(fts) && insideTkVol(cylinder)) {
    return getTkPropagator()->propagate(fts, cylinder);
  } else {
    return getGenPropagator()->propagate(fts, cylinder);
  }
}

std::pair<TrajectoryStateOnSurface,double>
SmartPropagatorWithIP::propagateWithPath(const FreeTrajectoryState& fts,
                                         const Plane& plane) const
{
  if (insideTkVol(fts) && insideTkVol(plane)) {
    return getTkPropagator()->propagateWithPath(fts, plane);
  } else {
    return getGenPropagator()->propagateWithPath(fts, plane);
  }
}

std::pair<TrajectoryStateOnSurface,double>
SmartPropagatorWithIP::propagateWithPath(const FreeTrajectoryState& fts,
                                         const Cylinder& cylinder) const
{
  if (insideTkVol(fts) && insideTkVol(cylinder)) {
    return getTkPropagator()->propagateWithPath(fts, cylinder);
  } else {
    return getGenPropagator()->propagateWithPath(fts, cylinder);
  }
}

bool SmartPropagatorWithIP::insideTkVol(const FreeTrajectoryState& fts) const
{
  return ( insideTkVol(fts.position()) );
}

bool SmartPropagatorWithIP::insideTkVol(const Surface& surface) const
{
  return ( insideTkVol(surface.position()) );
}

bool SmartPropagatorWithIP::insideTkVol( const BoundCylinder& cylin) const
{
  return ( insideTkVol(GlobalPoint(cylin.radius(),0.,(cylin.bounds().length())/2.)) );
}

bool SmartPropagatorWithIP::insideTkVol( const Plane& plane) const
{
  return ( insideTkVol(plane.position()) );
}

Propagator* SmartPropagatorWithIP::getTkPropagator() const
{
  return theTkProp_;
}

Propagator* SmartPropagatorWithIP::getGenPropagator() const
{
  return theGenProp_;
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{
  const reco::TransientTrack transientTrack = theBuilder_->build(&track);
  TrajectoryStateClosestToPoint traj = transientTrack.trajectoryStateClosestToPoint(vertex);
  if( traj.isValid() ) {
    return IP(traj.theState().momentum(),
              traj.perigeeParameters().transverseImpactParameter(),
              traj.perigeeError().transverseImpactParameterError(),
              traj.perigeeParameters().longitudinalImpactParameter(),
              traj.perigeeError().longitudinalImpactParameterError());
  }
  else {
    std::cout << "Invalid trajectoryStateClosestToPoint" << std::endl;
    return IP();
  }
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersOutsideTkVol(const FreeTrajectoryState & fts, const GlobalPoint & vertex) const
{
  SteppingHelixPropagator * steppingHelixProp = dynamic_cast<SteppingHelixPropagator*>(theGenProp_);
  if( steppingHelixProp != 0 ) {
    FreeTrajectoryState ftsPCA(steppingHelixProp->propagate(fts, vertex));
    GlobalPoint pca(ftsPCA.position());
    double dxy = sqrt(pow(pca.x()-vertex.x(), 2) + pow(pca.y()-vertex.y(), 2));
    double dz = pca.z() - vertex.z();
    // Ignore the errors for now
    return IP(ftsPCA.momentum(), dxy, 0., dz, 0.);
    // return IP(dxy, 0., dz, 0.);
  }
  return IP();
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersOutsideInTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{

  return computeImpactParametersOutsideTkVol(FreeTrajectoryState(GlobalPoint(track.innerPosition().x(),track.innerPosition().y(),track.innerPosition().z()),
								 GlobalVector(track.innerMomentum().x(),track.innerMomentum().y(),track.innerMomentum().z()),
								 TrackCharge(track.charge()), theField_),
                                             vertex);
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideOutTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{
  return computeImpactParametersOutsideTkVol(FreeTrajectoryState(GlobalPoint(track.outerPosition().x(),track.outerPosition().y(),track.outerPosition().z()),
								 GlobalVector(track.outerMomentum().x(),track.outerMomentum().y(),track.outerMomentum().z()),
								 TrackCharge(track.charge()), theField_),
					     vertex);
}
