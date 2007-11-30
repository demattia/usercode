#ifndef AnalysisExamples_SiStripDetectorPerformance_AnalyzedTrack_h
#define AnalysisExamples_SiStripDetectorPerformance_AnalyzedTrack_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/AnalyzedCluster.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <vector>
#include <boost/cstdint.hpp>

namespace anaobj {

  struct AnalyzedCluster;

  struct AnalyzedTrack{

    //  int   run;
    //  int   event;
    float momentum;
    float pt;
    int   charge;
    float eta;
    float phi;
    int   hitspertrack;
    float normchi2;
    float chi2;
    float ndof;

    int tk_id;
    std::vector<int> clu_id;
  
    float d0;
    float vx;
    float vy;
    float vz;
    float outerPt;

    float innermostPoint_X;
    float innermostPoint_Y;
    float innermostPoint_Z;
    float outermostPoint_X;
    float outermostPoint_Y;
    float outermostPoint_Z;

    int nrProjectedRecHits;
    int nrMatchedRecHits;
    int nrSingleRecHits;

    // valid and lost hits numbers
    short found;
    short lost;
    std::vector<std::pair<int, int> > hitvalidation;

    // Tracking particle quantities
    float simE;
    float simEt;
    float simP;
    float simPt;
    float simEta;
    float simPhi;
    float simAssocChi2;
    float simpdgId;
    float simptrackId;
  
    // Vector of Refs to AnalyzedClusters of the track
    //    std::vector< edm::Ref< std::vector <AnalyzedCluster> > >  vecRefClusterOwned;
    edm::RefVector<std::vector<AnalyzedCluster> > vecRefCluster;

    AnalyzedTrack() {

      // Initializations per track
      // -------------------------
      momentum        = -99.;  
      pt	      = -99.;  
      charge          = -99;  
      eta	      = -99.;  
      phi	      = -99.;  
      hitspertrack    = -99;  
      normchi2        = -99.;  
      chi2            = -99.;  
      ndof            = -99.;  

      tk_id           = -1;

      d0              = -99.;
      vx              = -99.;
      vy              = -99.;
      vz              = -99.;
      outerPt         = -99.;

      innermostPoint_X = -99;
      innermostPoint_Y = -99;
      innermostPoint_Z = -99;
      outermostPoint_X = -99;
      outermostPoint_Y = -99;
      outermostPoint_Z = -99;

      found           = -99;
      found           = -99;

      simP            = -99.;
      simPt           = -99.;
      simEta          = -99.;
      simPhi          = -99.;
      simAssocChi2    = -99.;
      simpdgId        = -99;
      simptrackId     = -99;
    }

    void ClusterRef( edm::Ref<std::vector<AnalyzedCluster> > anaclu_ref ) {
      vecRefCluster.push_back(anaclu_ref);
    }

    const edm::RefVector<std::vector<AnalyzedCluster> > GetClusterRefVec() const {
      return vecRefCluster;
    }

  };

  // AnalyzedTrack collection typedef
  typedef std::vector<AnalyzedTrack> AnalyzedTrackCollection;
  typedef edm::Ref<AnalyzedTrackCollection> AnalyzedTrackRef;
  typedef edm::RefProd<AnalyzedTrackCollection> AnalyzedTrackRefProd;
  typedef edm::RefVector<AnalyzedTrackCollection> AnalyzedTrackRefVector;

}

#endif
