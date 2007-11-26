
#include <memory>
#include <string>
#include <iostream>
#include <TMath.h>

#include "AnalysisExamples/SiStripDetectorPerformance/interface/TrackLocalAngleTIF.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
//#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoEnum.h"

using namespace std;

TrackLocalAngleTIF::TrackLocalAngleTIF() {
}

TrackLocalAngleTIF::~TrackLocalAngleTIF() {
}

void TrackLocalAngleTIF::init( const edm::EventSetup& es ) {
  //
  // get geometry
  //
  edm::ESHandle<TrackerGeometry> estracker;
  es.get<TrackerDigiGeometryRecord>().get(estracker);
  _tracker=&(* estracker);
}

std::vector<std::pair<const TrackingRecHit*,float> > TrackLocalAngleTIF::SeparateHits(reco::TrackInfoRef & trackinforef) {
  std::vector<std::pair<const TrackingRecHit*,float> >hitangleassociation;

  // Create the other objects. Use auto_ptr so that ownership is passed when the vector is taken and its memory is correctly freed
  oXZHitAngle      = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );
  oYZHitAngle      = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );
  oTrackLocalPitch = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );
             
  oLocalDir  = std::auto_ptr<HitLclDirAssociation>( new HitLclDirAssociation );
  oGlobalDir = std::auto_ptr<HitGlbDirAssociation>( new HitGlbDirAssociation );

  oLocalHitPos   = std::auto_ptr<LclPosAssociation>( new LclPosAssociation );
  oLocalTrackPos = std::auto_ptr<LclPosAssociation>( new LclPosAssociation );

  Hit3DAngle  = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );
  HitPhiAngle = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );

  HitValidation = std::auto_ptr<HitValidationAssociation>( new HitValidationAssociation );
  HitType       = std::auto_ptr<HitAngleAssociation>( new HitAngleAssociation );

  numberProjectedRecHits = 0;
  numberMatchedRecHits   = 0;
  numberSingleRecHits    = 0;

  for(_tkinfoiter=trackinforef->trajStateMap().begin();_tkinfoiter!=trackinforef->trajStateMap().end();++_tkinfoiter) {

    DetId detectorId = DetId(((*_tkinfoiter).first)->geographicalId());

    // Check if it is a valid hit
    if (((*_tkinfoiter).first)->isValid()) {

      HitValidation->push_back( std::make_pair( detectorId.rawId(), 1));

      const ProjectedSiStripRecHit2D* phit=dynamic_cast<const ProjectedSiStripRecHit2D*>(&(*(_tkinfoiter->first)));
      const SiStripMatchedRecHit2D* matchedhit=dynamic_cast<const SiStripMatchedRecHit2D*>(&(*(_tkinfoiter->first)));
      const SiStripRecHit2D* hit=dynamic_cast<const SiStripRecHit2D*>(&(*(_tkinfoiter->first)));
      //     LocalVector trackdirection=(_tkinfoiter->second.parameters()).momentum();
      //       LocalVector trackdirection=(trackinforef->stateOnDet((*_tkinfoiter).first).parameters()).momentum();
      //       LocalPoint  trackposition =(trackinforef->stateOnDet((*_tkinfoiter).first).parameters()).position();

      LocalVector trackdirection=(trackinforef->stateOnDet(reco::Updated, (*_tkinfoiter).first)->parameters()).momentum();
      LocalPoint  trackposition =(trackinforef->stateOnDet(reco::Updated, (*_tkinfoiter).first)->parameters()).position();
      // Projected Hit
      ////////////////
      if (phit) {
	//phit = POINTER TO THE PROJECTED RECHIT
	hit=&(phit->originalHit());

#ifdef DEBUG
	std::cout << "ProjectedHit found" << std::endl;
#endif

      }

      // Matched Hit
      //////////////
      if(matchedhit){//if matched hit...
	numberMatchedRecHits++;

#ifdef DEBUG
	std::cout<<"MatchedHit found"<<std::endl;
#endif

	GluedGeomDet * gdet=(GluedGeomDet *)_tracker->idToDet(matchedhit->geographicalId());

	GlobalVector gtrkdir=gdet->toGlobal(trackdirection);

#ifdef DEBUG
	std::cout<<"Track direction trasformed in global direction"<<std::endl;
#endif	

	//cluster and trackdirection on mono det
	
	// THIS THE POINTER TO THE MONO HIT OF A MATCHED HIT 
	const SiStripRecHit2D *monohit=matchedhit->monoHit();
 
	HitType->push_back( std::make_pair( monohit, 2 ) ); // 2==matched rechit (mono)
     
	const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > monocluster=monohit->cluster();
	const GeomDetUnit * monodet=gdet->monoDet();
      
	LocalVector monotkdir=monodet->toLocal(gtrkdir);
	//size=(monocluster->amplitudes()).size();
	// THE LOCAL ANGLE (MONO)
	if(monotkdir.z()!=0) {
	  float angle = atan(monotkdir.x()/ monotkdir.z())*180/TMath::Pi();
	  hitangleassociation.push_back(std::make_pair(monohit, angle));
	  oXZHitAngle->push_back( std::make_pair( monohit, atan( monotkdir.x()/ monotkdir.z())));
	  oYZHitAngle->push_back( std::make_pair( monohit, atan( monotkdir.y()/ monotkdir.z())));
	}
	else {

	    if ( monotkdir.x() != 0 ) {
	      float angle = ( monotkdir.x()/fabs(monotkdir.x()) )*180/TMath::Pi();
	      hitangleassociation.push_back(std::make_pair(monohit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( monohit, ( monotkdir.x()/fabs(monotkdir.x()) ) ) );
	    }
	    else { // angle and XZAngle is not defined in this case, fill with some fixed non physical value
	      float angle = -9999;
	      hitangleassociation.push_back(std::make_pair(monohit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( monohit, -9999 ) );
	    }
	    if ( monotkdir.y() != 0 ) {
	      oYZHitAngle->push_back( std::make_pair( monohit, ( monotkdir.y()/fabs(monotkdir.y()) ) ) );
	    }
	    else { // YZAngle is not defined in this case, fill with some fixed non physical value
	      oYZHitAngle->push_back( std::make_pair( monohit, -9999 ) );
	    }
	}
	//std::cout<<"Angle="<<atan(monotkdir.x(), monotkdir.z())*180/TMath::Pi()<<std::endl;

	// 3D angle and Phi angle
	Hit3DAngle->push_back( std::make_pair( monohit, acos( monotkdir.z()/ monotkdir.mag() ) ) );
	HitPhiAngle->push_back( std::make_pair( monohit, monotkdir.phi().value() ) );

	oLocalDir->push_back( std::make_pair( monohit, monotkdir));
	oGlobalDir->push_back( std::make_pair( monohit, gtrkdir));

	LocalPoint projectedPoint;
	projectedPoint = project(gdet,monodet,trackposition,monotkdir);
	oLocalHitPos->push_back( std::make_pair( monohit, monohit->localPosition()));
	oLocalTrackPos->push_back( std::make_pair( monohit, projectedPoint));

	const StripGeomDetUnit* theStripDetmono = dynamic_cast<const StripGeomDetUnit*>( (_tracker->idToDet(monohit->geographicalId())) );
	if (theStripDetmono)
	  {
	    const StripTopology* theStripTopol = dynamic_cast<const StripTopology*>( &(theStripDetmono->specificTopology()) );
	    oTrackLocalPitch->push_back( std::make_pair( monohit, theStripTopol->localPitch(projectedPoint)));
	  }
	else
	  {
	    oTrackLocalPitch->push_back( std::make_pair( monohit, -9999));
	  }
	
	//cluster and trackdirection on stereo det
	
	// THIS THE POINTER TO THE STEREO HIT OF A MATCHED HIT 
	const SiStripRecHit2D *stereohit=matchedhit->stereoHit();
	
	HitType->push_back( std::make_pair( stereohit, 3 ) ); // 3==matched rechit (stereo)

	const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > stereocluster=stereohit->cluster();
	const GeomDetUnit * stereodet=gdet->stereoDet(); 
	LocalVector stereotkdir=stereodet->toLocal(gtrkdir);
	//size=(stereocluster->amplitudes()).size();
	// THE LOCAL ANGLE (STEREO)
	if(stereotkdir.z()!=0) {
	  float angle = atan(stereotkdir.x()/ stereotkdir.z())*180/TMath::Pi();
	  hitangleassociation.push_back(std::make_pair(stereohit, angle));
	  oXZHitAngle->push_back( std::make_pair( stereohit, atan( stereotkdir.x()/ stereotkdir.z())));
	  oYZHitAngle->push_back( std::make_pair( stereohit, atan( stereotkdir.y()/ stereotkdir.z())));
	}
	else {

	    if ( stereotkdir.x() != 0 ) {
	      float angle = ( stereotkdir.x()/fabs(stereotkdir.x()) )*180/TMath::Pi();
	      hitangleassociation.push_back(std::make_pair(stereohit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( stereohit, ( stereotkdir.x()/fabs(stereotkdir.x()) ) ) );
	    }
	    else { // angle and XZAngle is not defined in this case, fill with some fixed non physical value
	      float angle = -9999;
	      hitangleassociation.push_back(std::make_pair(stereohit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( stereohit, -9999 ) );
	    }
	    if ( stereotkdir.y() != 0 ) {
	      oYZHitAngle->push_back( std::make_pair( stereohit, ( stereotkdir.y()/fabs(stereotkdir.y()) ) ) );
	    }
	    else { // YZAngle is not defined in this case, fill with some fixed non physical value
	      oYZHitAngle->push_back( std::make_pair( stereohit, -9999 ) );
	    }
	}
	//std::cout<<"Angle="<<atan(monotkdir.x(), monotkdir.z())*180/TMath::Pi()<<std::endl;

	// 3D angle and Phi angle
	Hit3DAngle->push_back( std::make_pair( stereohit, acos( stereotkdir.z()/ stereotkdir.mag() ) ) );
	HitPhiAngle->push_back( std::make_pair( stereohit, stereotkdir.phi().value() ) );

	oLocalDir->push_back( std::make_pair( stereohit, stereotkdir));
	oGlobalDir->push_back( std::make_pair( stereohit, gtrkdir));

	projectedPoint = project(gdet,stereodet,trackposition,stereotkdir);
	oLocalHitPos->push_back( std::make_pair( stereohit, stereohit->localPosition()));
	oLocalTrackPos->push_back( std::make_pair( stereohit, projectedPoint));

	const StripGeomDetUnit*  theStripDetstereo = dynamic_cast<const StripGeomDetUnit*>( (_tracker->idToDet(stereohit->geographicalId())) );
	if (theStripDetstereo)
	  {
	    const StripTopology* theStripTopol = dynamic_cast<const StripTopology*>( &(theStripDetstereo->specificTopology()) );
	    oTrackLocalPitch->push_back( std::make_pair( stereohit, theStripTopol->localPitch(projectedPoint)));
	  }
	else
	  {
	    oTrackLocalPitch->push_back( std::make_pair( stereohit, -9999));
	  }
      }

      // Single sided detector Hit
      ////////////////////////////
      else if( hit ) {
	//  hit= POINTER TO THE RECHIT
	if (phit) {
	  HitType->push_back( std::make_pair( hit, 1 ) ); // 1==projected rechit
	  numberProjectedRecHits++;
	}
	else {
	  HitType->push_back( std::make_pair( hit, 4 ) ); // 4==single rechit
	  numberSingleRecHits++;
	}

	const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > cluster=hit->cluster();
	GeomDet * gdet=(GeomDet *)_tracker->idToDet(hit->geographicalId());
	//size=(cluster->amplitudes()).size();
	GlobalVector gtrkdir=gdet->toGlobal(trackdirection);
      
	// THE LOCAL ANGLE
	if(trackdirection.z()!=0) {
	  float angle = atan(trackdirection.x()/ trackdirection.z())*180/TMath::Pi();
	  hitangleassociation.push_back(std::make_pair(hit, angle));
	  oXZHitAngle->push_back( std::make_pair( hit, atan( trackdirection.x()/ trackdirection.z())));
	  oYZHitAngle->push_back( std::make_pair( hit, atan( trackdirection.y()/ trackdirection.z())));
	}
	else {

	    if ( trackdirection.x() != 0 ) {
	      float angle = ( trackdirection.x()/fabs(trackdirection.x()) )*180/TMath::Pi();
	      hitangleassociation.push_back(std::make_pair(hit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( hit, ( trackdirection.x()/fabs(trackdirection.x()) ) ) );
	    }
	    else { // angle and XZAngle is not defined in this case, fill with some fixed non physical value
	      float angle = -9999;
	      hitangleassociation.push_back(std::make_pair(hit, angle)); 
	      oXZHitAngle->push_back( std::make_pair( hit, -9999 ) );
	    }
	    if ( trackdirection.y() != 0 ) {
	      oYZHitAngle->push_back( std::make_pair( hit, ( trackdirection.y()/fabs(trackdirection.y()) ) ) );
	    }
	    else { // YZAngle is not defined in this case, fill with some fixed non physical value
	      oYZHitAngle->push_back( std::make_pair( hit, -9999 ) );
	    }
	}

	// 3D angle and Phi angle
	Hit3DAngle->push_back( std::make_pair( hit, acos( trackdirection.z()/ trackdirection.mag() ) ) );
	HitPhiAngle->push_back( std::make_pair( hit, trackdirection.phi().value() ) );

	oLocalDir->push_back( std::make_pair( hit, trackdirection));
	oGlobalDir->push_back( std::make_pair( hit, gtrkdir));

	oLocalHitPos->push_back( std::make_pair( hit, hit->localPosition()));
	oLocalTrackPos->push_back( std::make_pair( hit, trackposition));

	const StripGeomDetUnit*  theStripDet = dynamic_cast<const StripGeomDetUnit*>( (_tracker->idToDet(detectorId)) );
	if (theStripDet)
	  {
	    const StripTopology* theStripTopol = dynamic_cast<const StripTopology*>( &(theStripDet->specificTopology()) );
	    oTrackLocalPitch->push_back( std::make_pair( hit, theStripTopol->localPitch(trackposition)));
	  }
	else
	  {
	    oTrackLocalPitch->push_back( std::make_pair( hit, -9999));
	  }
      }
      else {
#ifdef DEBUG
	std::cout << "not matched, mono or projected rechit" << std::endl;
#endif
      }
    } // end if is valid hit
    else { // hit is not valid
      HitValidation->push_back( std::make_pair( detectorId.rawId(), 0));
    }
  } // end loop on rechits
  return (hitangleassociation);
}

LocalPoint TrackLocalAngleTIF::project(const GeomDet *det,const GeomDet* projdet,LocalPoint position,LocalVector trackdirection)const
{
  
  GlobalPoint globalpoint=(det->surface()).toGlobal(position);
  
  // position of the initial and final point of the strip in glued local coordinates
  LocalPoint projposition=(projdet->surface()).toLocal(globalpoint);
  
  //correct the position with the track direction
  
  float scale=-projposition.z()/trackdirection.z();
  
  projposition+= scale*trackdirection;
  
  return projposition;
}
