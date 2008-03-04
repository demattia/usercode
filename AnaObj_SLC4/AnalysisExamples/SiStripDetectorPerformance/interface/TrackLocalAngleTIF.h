// ---------------------------------------------------------------------------
// This class performs a simple disentaglment of the matched rechits
// providing a vector of rechits and angles.
//
// Class objects need to be initialized with the event setup. It needs
// to access the tracker geometry to convert the local directions.
// The SeparateHits method requires a reco::TrackInfoRef. It takes it
// by reference. It returns a vector<pair<const TrackingRecHit*,float> >,
// which contains all the hits (matched hits are divided in mono and stereo)
// and the corresponding track angles.
// It also stores the phi and theta angles and the local and global direction
// vectors which can be taken by the corresponding get methods.
//
// M.De Mattia 13/2/2007
// ---------------------------------------------------------------------------

#ifndef AnalysisExamples_SiStripDetectorPerformance_TrackLocalAngleTIF_h
#define AnalysisExamples_SiStripDetectorPerformance_TrackLocalAngleTIF_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"

// For the auto pointer
#include <memory>

class TrackLocalAngleTIF 
{
 public:
  typedef std::vector<std::pair<const TrackingRecHit *, float> > HitAngleAssociation;
  typedef std::vector<std::pair<const TrackingRecHit *, LocalVector > > HitLclDirAssociation;
  typedef std::vector<std::pair<const TrackingRecHit *, GlobalVector> > HitGlbDirAssociation;
  typedef std::vector<std::pair<const TrackingRecHit *, LocalPoint > > LclPosAssociation;
  typedef std::vector<std::pair<int, int > > HitValidationAssociation;
  
 private:

  std::auto_ptr<HitAngleAssociation> oXZHitAngle;
  std::auto_ptr<HitAngleAssociation> oYZHitAngle;
  std::auto_ptr<HitAngleAssociation> oTrackLocalPitch;

  std::auto_ptr<HitLclDirAssociation> oLocalDir;
  std::auto_ptr<HitGlbDirAssociation> oGlobalDir;

  std::auto_ptr<LclPosAssociation> oLocalHitPos;
  std::auto_ptr<LclPosAssociation> oLocalTrackPos;

  std::auto_ptr<HitAngleAssociation> Hit3DAngle;
  std::auto_ptr<HitAngleAssociation> HitPhiAngle;

  std::auto_ptr<HitValidationAssociation> HitValidation;
  std::auto_ptr<HitAngleAssociation> HitType;

  int numberProjectedRecHits;
  int numberMatchedRecHits;
  int numberSingleRecHits;

  const TrackerGeometry * _tracker;
  reco::TrackInfo::TrajectoryInfo::const_iterator _tkinfoiter;

 public:
  TrackLocalAngleTIF();
  
  virtual ~TrackLocalAngleTIF();

  void init ( const edm::EventSetup& es  );

  std::vector<std::pair<const TrackingRecHit*,float> > SeparateHits(reco::TrackInfoRef & trackinforef);

  LocalPoint project(const GeomDet *det,const GeomDet* projdet,LocalPoint position,LocalVector trackdirection)const;

/*   inline HitAngleAssociation getXZHitAngle() const throw() {  */
/*     return (oXZHitAngle); } */
/*   inline HitAngleAssociation getYZHitAngle() const throw() {  */
/*     return (oYZHitAngle); } */

/*   inline HitLclDirAssociation getLocalDir() const throw() { */
/*     return (oLocalDir); } */
/*   inline HitGlbDirAssociation getGlobalDir() const throw() { */
/*     return (oGlobalDir); } */

/*   inline HitAngleAssociation getHit3DAngle() const throw() { */
/*     return (Hit3DAngle); } */

  inline std::auto_ptr<HitAngleAssociation> getXZHitAngle() { 
    return oXZHitAngle; }
  inline std::auto_ptr<HitAngleAssociation> getYZHitAngle() { 
    return (oYZHitAngle); }
  inline std::auto_ptr<HitAngleAssociation> getTrackLocalPitch() { 
    return (oTrackLocalPitch); }

  inline std::auto_ptr<HitLclDirAssociation> getLocalDir() {
    return (oLocalDir); }
  inline std::auto_ptr<HitGlbDirAssociation> getGlobalDir() {
    return (oGlobalDir); }

  inline std::auto_ptr<LclPosAssociation> getLocalHitPos() {
    return (oLocalHitPos); }
  inline std::auto_ptr<LclPosAssociation> getLocalTrackPos() {
    return (oLocalTrackPos); }

  inline std::auto_ptr<HitAngleAssociation> getHit3DAngle() {
    return (Hit3DAngle); }
  inline std::auto_ptr<HitAngleAssociation> getHitPhiAngle() {
    return (HitPhiAngle); }

  inline std::auto_ptr<HitValidationAssociation> getHitValidation() {
    return (HitValidation); }
  inline std::auto_ptr<HitAngleAssociation> getHitType() { 
    return HitType; }

  inline int getNrProjectedRecHits() {
    return numberProjectedRecHits; }
  inline int getNrMatchedRecHits() {
    return numberMatchedRecHits; }
  inline int getNrSingleRecHits() {
    return numberSingleRecHits; }
};


#endif // AnalysisExamples_SiStripDetectorPerformance_TrackLocalAngleTIF_h
