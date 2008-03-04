#ifndef AnalysisExamples_SiStripDetectorPerformance_AnalyzedCluster_h
#define AnalysisExamples_SiStripDetectorPerformance_AnalyzedCluster_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/AnalyzedTrack.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <vector>
#include <boost/cstdint.hpp>

namespace anaobj {

  struct AnalyzedTrack;

  struct AnalyzedCluster {

    int run;
    int event;
    int size;
    int rawId;
    int module;
    int string;
    int rod;
    int extint;
    int bwfw;
    int ring;
    int wheel;
    int type;
    int layer;
    int monostereo;

    std::vector<int> clu_id;

    float    clusterpos;
    float    clustereta;
    float    clusterchg;
    float    clusterchgl;
    float    clusterchgr;
    float    clusterchg_tickcorr;
    float    clustertick;
    float    clusternoise;
    float    clustermaxchg;
    float    clustercrosstalk;
    uint32_t geoId;
    uint16_t firstStrip;
    float    clusterbarycenter;
    float    clusterseednoise;
    double   thickness;
    double   pitch;
    std::vector<float> clusterstripnoises;

    // Cluster position
    float LclPos_X;
    float LclPos_Y;
    float LclPos_Z;
    float GlbPos_X;
    float GlbPos_Y;
    float GlbPos_Z;

    // Vectors of digis to the left and right of the cluster
    std::vector<float> rawDigiAmplitudesL;
    std::vector<float> rawDigiAmplitudesR;

    // Maps with key tk_id
    std::map<int, int>   sign;
    std::map<int, float> angle;
    std::map<int, float> angleXZ;
    std::map<int, float> angleYZ;
    std::map<int, float> angle3D;
    std::map<int, float> anglePhi;
    std::map<int, float> tk_phi;
    std::map<int, float> tk_theta;
//    std::map<int, int>   tk_id;
    std::vector<int>     tk_id;
    std::map<int, float> stereocorrection;
    std::map<int, float> localmagfield;
    // Positions
/*     std::map<int, float> LclPos_X; */
/*     std::map<int, float> LclPos_Y; */
/*     std::map<int, float> LclPos_Z; */
/*     std::map<int, float> GlbPos_X; */
/*     std::map<int, float> GlbPos_Y; */
/*     std::map<int, float> GlbPos_Z; */
    // Directions
    std::map<int, float> LclDir_X;
    std::map<int, float> LclDir_Y;
    std::map<int, float> LclDir_Z;
    std::map<int, float> GlbDir_X;
    std::map<int, float> GlbDir_Y;
    std::map<int, float> GlbDir_Z;

    std::map<int, float> HitLclPos_X; // for TIB and TOB in local frame, for TID and TEC in measurement frame
    std::map<int, float> HitLclPos_Y;
    std::map<int, float> HitLclPos_Z;
    std::map<int, float> TrackLclPos_X;
    std::map<int, float> TrackLclPos_Y;
    std::map<int, float> TrackLclPos_Z;
    std::map<int, float> TrackLclPitch;
    std::map<int, int>   HitType; // 1==projected, 2==matched (mono), 3==matched (stereo), 4==single

    //    std::vector< edm::Ref< std::vector<AnalyzedTrack*> > > vecRefTrackBelonged;
    edm::RefVector<std::vector<AnalyzedTrack> > vecRefTrack;

    AnalyzedCluster() {

      run          = -99;
      event        = -99;
      size         = -99;
      rawId        = -99;
      module       = -99;
      string       = -99;
      rod          = -99;
      extint       = -99;
      bwfw         = -99;
      ring         = -99;
      wheel        = -99;
      type         = -99;
      layer        = -99;
      monostereo   = -99;       

//      clu_id = -1;

      clusterpos          = -99;
      clustereta          = -99;
      clusterchg          = -99;
      clusterchgl         = -99;
      clusterchgr         = -99;
      clusterchg_tickcorr = -99;
      clustertick         = -99;
      clusternoise        = -99;
      clustermaxchg       = -99;
      clustercrosstalk    = -99;
      geoId               = 0;
      firstStrip          = 0;
      clusterbarycenter   = -99;

      LclPos_X = -999;
      LclPos_Y = -999;
      LclPos_Z = -999;
      GlbPos_X = -999;
      GlbPos_Y = -999;
      GlbPos_Z = -999;

      //   clusterstripnoises.clear();
      // 
      //   sign.clear();
      //   angle.clear();            
      //   tk_phi.clear();           
      //   tk_theta.clear();         
      //   tk_id.clear();            
      //   stereocorrection.clear(); 
      //   localmagfield.clear();    
      //   dLclX.clear();
      //   dLclY.clear();
      //   dLclZ.clear();
      //   dGlbX.clear();
      //   dGlbY.clear();
      //   dGlbZ.clear();

    }

    void TrackRef( edm::Ref<std::vector<AnalyzedTrack> > anatk_ref ) {
      vecRefTrack.push_back(anatk_ref);
    }

    const edm::RefVector<std::vector<AnalyzedTrack> > GetTrackRefVec() const {
      return vecRefTrack;
    }

  };

  // AnalyzedTrack collection typedef
  typedef std::vector<AnalyzedCluster> AnalyzedClusterCollection;
  typedef edm::Ref<AnalyzedClusterCollection> AnalyzedClusterRef;
  typedef edm::RefProd<AnalyzedClusterCollection> AnalyzedClusterRefProd;
  typedef edm::RefVector<AnalyzedClusterCollection> AnalyzedClusterRefVector;
}

#endif
