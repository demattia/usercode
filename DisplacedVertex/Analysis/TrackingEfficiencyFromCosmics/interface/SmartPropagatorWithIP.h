#ifndef SmartPropagatorWithIP_H
#define SmartPropagatorWithIP_H

/** \class SmartPropagatorWithIP
 *
 * Starting from the SmartPropagatorWithIP in TrackingTools/GeomPropagators
 * this class adds the methods to compute the point of closest approach and the IP.
 *
 * \author  M. De Mattia - 26/10/2011
 *
 */

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/GeomPropagators/interface/TrackerBounds.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include <TVector3.h>

class BoundCylinder;
class BoundPlane;

/* Class SmartPropagatorWithIP Interface */

class SmartPropagatorWithIP : public Propagator
{
public:

  /* Constructor */
  ///Defines which propagator is used inside Tk and which outside
  SmartPropagatorWithIP(Propagator* aTkProp, Propagator* aGenProp, const MagneticField* field,
                        const TransientTrackBuilder * theBuilder,
                        PropagationDirection dir = alongMomentum, float epsilon = 5);

  ///Defines which propagator is used inside Tk and which outside
  SmartPropagatorWithIP(const Propagator& aTkProp, const Propagator& aGenProp,const MagneticField* field,
                        const TransientTrackBuilder * theBuilder,
                        PropagationDirection dir = alongMomentum, float epsilon = 5);

  ///Copy constructor
  SmartPropagatorWithIP( const SmartPropagatorWithIP& );

  /** virtual destructor */
  virtual ~SmartPropagatorWithIP();

  ///Virtual constructor (using copy c'tor)
  virtual SmartPropagatorWithIP* clone() const
  {
    return new SmartPropagatorWithIP(getTkPropagator(),getGenPropagator(),magneticField(),theBuilder_,propagationDirection(),epsilon_);
  }

  ///setting the direction fo both components
  void setPropagationDirection (PropagationDirection dir) const
  {
    Propagator::setPropagationDirection (dir);
    getTkPropagator()->setPropagationDirection(dir);
    getGenPropagator()->setPropagationDirection(dir);
  }

  /* Operations as propagator*/
  TrajectoryStateOnSurface propagate(const FreeTrajectoryState& fts,
                                     const Surface& surface) const;

  TrajectoryStateOnSurface propagate(const TrajectoryStateOnSurface& tsos,
                                     const Surface& surface) const
  {
    return Propagator::propagate(tsos,surface);
  }

  TrajectoryStateOnSurface propagate(const FreeTrajectoryState& fts,
                                     const Plane& plane) const;

  TrajectoryStateOnSurface propagate(const TrajectoryStateOnSurface& tsos,
                                     const Plane& plane) const
  {
    return Propagator::propagate(tsos, plane);
  }

  TrajectoryStateOnSurface propagate(const FreeTrajectoryState& fts,
                                     const Cylinder& cylinder) const;

  TrajectoryStateOnSurface propagate(const TrajectoryStateOnSurface& tsos,
                                     const Cylinder& cylinder) const
  {
    return Propagator::propagate(tsos, cylinder);
  }

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const FreeTrajectoryState& fts,
                                                               const Surface& surface) const
  {
    return Propagator::propagateWithPath(fts,surface);
  }

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const TrajectoryStateOnSurface& tsos,
                                                               const Surface& surface) const
  {
    return Propagator::propagateWithPath(tsos,surface);
  }

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const FreeTrajectoryState& fts,
                                                               const Plane& plane) const;

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const TrajectoryStateOnSurface& tsos,
                                                               const Plane& plane) const
  {
    return Propagator::propagateWithPath(tsos, plane);
  }

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const FreeTrajectoryState& fts,
                                                               const Cylinder& cylinder) const;

  std::pair<TrajectoryStateOnSurface,double> propagateWithPath(const TrajectoryStateOnSurface& tsos,
                                                               const Cylinder& cylinder) const
  {
    return Propagator::propagateWithPath(tsos, cylinder);
  }

  ///true if a fts is inside tracker volume
  bool insideTkVol(const FreeTrajectoryState& fts) const;
  ///true if a surface is inside tracker volume
  bool insideTkVol(const Surface& surface) const;
  ///true if a cylinder is inside tracker volume
  bool insideTkVol(const BoundCylinder& cylin) const;
  ///true if a plane is inside tracker volume
  bool insideTkVol(const Plane& plane) const;
  ///true if a point is inside tracker volume
  inline bool insideTkVol(const GlobalPoint & gp) const
  {
    return (gp.perp()<= TrackerBounds::radius()+epsilon_) && (fabs(gp.z())<= TrackerBounds::halfLength()+epsilon_);
  }

  ///return the propagator used inside tracker
  Propagator* getTkPropagator() const;
  ///return the propagator used outside tracker
  Propagator* getGenPropagator() const;
  ///return the magneticField
  virtual const MagneticField* magneticField() const {return theField_;}

  struct IP
  {
    IP(): dxyValue(65535.), dxyError(65535.), dzValue(65535.), dzError(65535.), pt(0.), eta(0.), phi(0.), ptError(0.), etaError(0.), phiError(0.)
    {}
    IP(const double & inputPt, const double & inputPtError,
       const double & inputEta, const double & inputEtaError,
       const double & inputPhi, const double & inputPhiError,
       const double & inputDxyValue, const double & inputDxyError,
       const double & inputDzValue, const double & inputDzError) :
      dxyValue(inputDxyValue), dxyError(inputDxyError), dzValue(inputDzValue), dzError(inputDzError),
      pt(inputPt), eta(inputEta), phi(inputPhi), ptError(inputPtError), etaError(inputEtaError), phiError(inputPhiError)
    {}
    template <class T>
    IP(const T & momentum, const double & inputDxyValue, const double & inputDxyError, const double & inputDzValue, const double & inputDzError) :
      dxyValue(inputDxyValue), dxyError(inputDxyError), dzValue(inputDzValue), dzError(inputDzError)
    {
      pt = momentum.perp();
      eta = momentum.eta();
      phi = momentum.phi();
      // Do not have the error available in the momentum, need to check how to get it from the propagated track.
      ptError = 0.;
      etaError = 0.;
      phiError = 0.;
    }
    double dxyValue, dxyError;
    double dzValue, dzError;
    double pt, eta, phi;
    double ptError, etaError, phiError;
  };

  // IP methods
  // ----------
  template <class T>
  SmartPropagatorWithIP::IP computeGenImpactParametersInsideTkVol( const T & track, const math::XYZPoint & genVertex,
                                                                   const int genCharge, const GlobalPoint & vertex ) const;
  template <class T>
  SmartPropagatorWithIP::IP computeGenImpactParametersOutsideTkVol( const T & track, const math::XYZPoint & genVertex,
                                                                    const int genCharge, const GlobalPoint & vertex ) const;
  // template <class T>
  // SmartPropagatorWithIP::IP computeImpactParametersInsideTkVol( const T & track, const GlobalPoint & vertex ) const;
  SmartPropagatorWithIP::IP computeImpactParametersInsideTkVol( const reco::Track & track, const GlobalPoint & vertex ) const;
  /// This is more for internal use, but leaving it public because it is generic enough that it might be useful
  SmartPropagatorWithIP::IP computeImpactParametersOutsideTkVol(const FreeTrajectoryState & fts, const GlobalPoint & vertex) const;
  SmartPropagatorWithIP::IP computeImpactParametersOutsideInTkVol( const FreeTrajectoryState & tsos, const GlobalPoint & vertex ) const;
  SmartPropagatorWithIP::IP computeImpactParametersOutsideInTkVol( const reco::Track & track, const GlobalPoint & vertex ) const;
  SmartPropagatorWithIP::IP computeImpactParametersInsideOutTkVol( const FreeTrajectoryState & tsos, const GlobalPoint & vertex ) const;
  SmartPropagatorWithIP::IP computeImpactParametersInsideOutTkVol( const reco::Track & track, const GlobalPoint & vertex ) const;
  SmartPropagatorWithIP::IP computeImpactParameters( const reco::Track & track, const GlobalPoint & vertex ) const;

  static void initTkVolume(float epsilon);

  inline double epsilon() const {return epsilon_;}
  inline TransverseImpactPointExtrapolator * transverseExtrapolator() const {return transverseExtrapolator_;}
  inline AnalyticalImpactPointExtrapolator * analyticalExtrapolator() const {return analyticalExtrapolator_;}
  inline const TransientTrackBuilder * theBuilder() const {return theBuilder_;}

private:

  template <class T>
  FreeTrajectoryState ftsAtProduction(const T & track, const math::XYZPoint & genVertex, const int genCharge) const;

  mutable Propagator * theTkProp_;
  mutable Propagator * theGenProp_;
  const MagneticField * theField_;
  mutable float epsilon_;
  mutable TransverseImpactPointExtrapolator * transverseExtrapolator_;
  mutable AnalyticalImpactPointExtrapolator * analyticalExtrapolator_;
  const TransientTrackBuilder * theBuilder_;
  mutable AlgebraicSymMatrix55 nullCovariance_;
  static ReferenceCountingPointer<Cylinder> & theTkVolume();

protected:
};

template <class T>
FreeTrajectoryState SmartPropagatorWithIP::ftsAtProduction(const T & track, const math::XYZPoint & genVertex, const int genCharge) const
{
  TVector3 genMomentum(0,0,0);
  genMomentum.SetPtEtaPhi(track.pt(),track.eta(),track.phi());
  return FreeTrajectoryState(GlobalPoint(genVertex.x(),genVertex.y(),genVertex.z()),
                             GlobalVector(genMomentum.x(),genMomentum.y(),genMomentum.z()),
                             TrackCharge(genCharge), theField_);
}

template <class T>
SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeGenImpactParametersInsideTkVol( const T & track, const math::XYZPoint & genVertex,
                                                                                        const int genCharge, const GlobalPoint & vertex ) const
{
  VertexDistanceXY distXY;
  // VertexDistance3D dist3D;

  FreeTrajectoryState ftsAtProd = ftsAtProduction(track, genVertex, genCharge);
  // TrajectoryStateOnSurface analyticalTSOS_ = analyticalExtrapolator_->extrapolate(ftsAtProd, vertex);
  TrajectoryStateOnSurface transverseTSOS_ = transverseExtrapolator_->extrapolate(ftsAtProd, vertex);

  if(transverseTSOS_.isValid()) {
    TrajectoryStateOnSurface transverseTSOS(transverseTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_),
                                            transverseTSOS_.surface(), transverseTSOS_.magneticField(), transverseTSOS_.weight());
    // Create a dummy vertex
    reco::Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    reco::Vertex::Point p(vertex.x(), vertex.y(), vertex.z());
    reco::Vertex dummyVertex(p, e, 0, 0, 0);

    std::pair<bool,Measurement1D> dxy = IPTools::absoluteImpactParameter(transverseTSOS, dummyVertex, distXY);
    return IP(transverseTSOS.freeTrajectoryState()->momentum(),
              dxy.second.value(), dxy.second.error(),
              transverseTSOS.globalPosition().z(), transverseTSOS.cartesianError().position().czz());
  }
  return IP();

//  // This is for the 3D impact parameter
//  // activate if needed.

//  if(analyticalTSOS_.isValid()) {
//    // analyticalTSOS_ has no errors defined. Explicitly set the errors to 0 for the genparticle state
//    TrajectoryStateOnSurface analyticalTSOS(analyticalTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_),
//                                            analyticalTSOS_.surface(), analyticalTSOS_.magneticField(), analyticalTSOS_.weight());
//    std::pair<bool,Measurement1D> dxyz = IPTools::absoluteImpactParameter(analyticalTSOS, reco::Vertex(vertex), dist3D);
//    // dxyz_.first = dxyz.second.value();
//    // dxyz_.second = dxyz.second.error();
//  }
//  else {
//    std::cout << "Invalid trajectoryStateClosestToPoint for analyticalExtrapolator for GEN" << std::endl;
//    // dxyz_.first = 65535;
//    // dxyz_.second = 65535;
//  }
}

template <class T>
SmartPropagatorWithIP::IP SmartPropagatorWithIP::computeGenImpactParametersOutsideTkVol( const T & track, const math::XYZPoint & genVertex,
                                                                                         const int genCharge, const GlobalPoint & vertex ) const
{
  return computeImpactParametersOutsideTkVol(ftsAtProduction(track, genVertex, genCharge), vertex);
}
#endif // SmartPropagatorWithIP_H
