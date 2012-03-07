/* This Class Header */
#include "Analysis/SmartPropagatorWithIP/interface/SmartPropagatorWithIP.h"

/* Collaborating Class Header */
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"

#include <DataFormats/GeometrySurface/interface/BoundCylinder.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

/* Base Class Headers */

/* C++ Headers */

/* ====================================================================== */

ReferenceCountingPointer<Cylinder> & SmartPropagatorWithIP::theTkVolume()
{
  static ReferenceCountingPointer<Cylinder> local=0;
  return local;
}


/* Constructor */
SmartPropagatorWithIP::SmartPropagatorWithIP(Propagator* aTkProp, Propagator* aGenProp, const MagneticField* field,
                                             const TransientTrackBuilder * theBuilder,
                                             PropagationDirection dir, float epsilon) :
  Propagator(dir), theTkProp_(aTkProp->clone()), theGenProp_(aGenProp->clone()), theField_(field), epsilon_(epsilon),
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
  if (theTkVolume()==0) initTkVolume(epsilon);
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
  if (theTkVolume()==0) initTkVolume(epsilon);
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
void SmartPropagatorWithIP::initTkVolume(float epsilon)
{
  // fill tracker dimensions
  Surface::PositionType pos(0,0,0); // centered at the global origin
  Surface::RotationType rot; // unit matrix - barrel cylinder orientation
  theTkVolume() = Cylinder::build(pos, rot, TrackerBounds::radius()+epsilon);
}

TrajectoryStateOnSurface SmartPropagatorWithIP::propagate(const FreeTrajectoryState& fts,
                                                          const Surface& surface) const
{
  return Propagator::propagate(fts, surface);
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

//template <class T>
//SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideTkVol( const T & track, const GlobalPoint & vertex ) const
SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{
  // std::cout << "propagating inside tkVolume" << std::endl;
  const reco::TransientTrack transientTrack = theBuilder_->build(track);
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
  // std::cout << "propagating outside tkVolume" << std::endl;
  SteppingHelixPropagator * steppingHelixProp = dynamic_cast<SteppingHelixPropagator*>(theGenProp_);
  if( steppingHelixProp != 0 ) {
    FreeTrajectoryState ftsPCA(steppingHelixProp->propagate(fts, vertex));
    GlobalPoint pca(ftsPCA.position());
    // Taking only the diagonal terms. This is a semplification, should be checked.
    double dxy = sqrt(pow(pca.x()-vertex.x(), 2) + pow(pca.y()-vertex.y(), 2));
    double dz = pca.z() - vertex.z();
    // Ignore the errors for now
    return IP(ftsPCA.momentum(), dxy, 0., dz, 0.);
    // return IP(dxy, 0., dz, 0.);
  }
  return IP();
}

//SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersOutsideInTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
//{
//  return computeImpactParametersOutsideTkVol(FreeTrajectoryState(GlobalPoint(track.innerPosition().x(),track.innerPosition().y(),track.innerPosition().z()),
//                                                                 GlobalVector(track.innerMomentum().x(),track.innerMomentum().y(),track.innerMomentum().z()),
//                                                                 TrackCharge(track.charge()), theField_),
//                                             vertex);
//}

//SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideOutTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
//{
//  return computeImpactParametersOutsideTkVol(FreeTrajectoryState(GlobalPoint(track.outerPosition().x(),track.outerPosition().y(),track.outerPosition().z()),
//                                                                 GlobalVector(track.outerMomentum().x(),track.outerMomentum().y(),track.outerMomentum().z()),
//                                                                 TrackCharge(track.charge()), theField_),
//                                             vertex);
//}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersOutsideInTkVol( const FreeTrajectoryState & fts, const GlobalPoint & vertex ) const
{
  TrajectoryStateOnSurface tsos(propagate(fts, *(theTkVolume().get())));
  if( tsos.isValid() ) {
    // This does not work because the errors are missing. For now use the outsideTkVol.
    // return computeImpactParametersInsideTkVol(*(tsos.freeTrajectoryState()), vertex);
    return computeImpactParametersOutsideTkVol(fts, vertex);
  }
  // std::cout << "Mixed propagation failed, trying again with stepping helix only" << std::endl;
  return computeImpactParametersOutsideTkVol(fts, vertex);
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersOutsideInTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{
  return computeImpactParametersOutsideInTkVol(FreeTrajectoryState(GlobalPoint(track.innerPosition().x(),track.innerPosition().y(),track.innerPosition().z()),
                                                                   GlobalVector(track.innerMomentum().x(),track.innerMomentum().y(),track.innerMomentum().z()),
                                                                   TrackCharge(track.charge()), theField_), vertex);
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideOutTkVol( const FreeTrajectoryState & fts, const GlobalPoint & vertex ) const
{
  TrajectoryStateOnSurface tsos(propagate(fts, *(theTkVolume().get())));
  if( tsos.isValid() ) {
    return computeImpactParametersOutsideTkVol(*(tsos.freeTrajectoryState()), vertex);
  }
  // std::cout << "Mixed propagation failed, trying again with stepping helix only" << std::endl;
  return computeImpactParametersOutsideTkVol(fts, vertex);
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParametersInsideOutTkVol( const reco::Track & track, const GlobalPoint & vertex ) const
{
  return computeImpactParametersInsideOutTkVol(FreeTrajectoryState(GlobalPoint(track.outerPosition().x(),track.outerPosition().y(),track.outerPosition().z()),
                                                                   GlobalVector(track.outerMomentum().x(),track.outerMomentum().y(),track.outerMomentum().z()),
                                                                   TrackCharge(track.charge()), theField_), vertex);
}

SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeImpactParameters( const reco::Track & track, const GlobalPoint & vertex ) const
{
  // When the vertex is inside the tkVolume
  if( insideTkVol(vertex) ) {
    // std::cout << "Vertex inside tkVolume" << std::endl;
    FreeTrajectoryState innerTsos(GlobalPoint(track.innerPosition().x(),track.innerPosition().y(),track.innerPosition().z()),
                                  GlobalVector(track.innerMomentum().x(),track.innerMomentum().y(),track.innerMomentum().z()),
                                  TrackCharge(track.charge()), theField_);
    if( insideTkVol(innerTsos) ) return computeImpactParametersInsideTkVol(track, vertex);
    // return computeImpactParametersInsideTkVol(track, vertex);
    return computeImpactParametersOutsideInTkVol(innerTsos, vertex);
  }
  else {
    // std::cout << "Vertex outside tkVolume" << std::endl;
    FreeTrajectoryState outerTsos(GlobalPoint(track.outerPosition().x(),track.outerPosition().y(),track.outerPosition().z()),
                                  GlobalVector(track.outerMomentum().x(),track.outerMomentum().y(),track.outerMomentum().z()),
                                  TrackCharge(track.charge()), theField_);
    if( insideTkVol(outerTsos) ) return computeImpactParametersInsideOutTkVol(outerTsos, vertex);
    return computeImpactParametersOutsideTkVol(outerTsos, vertex);
  }
}
