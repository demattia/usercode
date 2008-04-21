#ifndef REJECTIONALG_H
#define REJECTIONALG_H

#include <vector>
#include <algorithm>

#include "AnalysisExamples/L1PixelAnalyzer/interface/WAverager.h"
#include "AnalysisExamples/L1PixelAnalyzer/interface/RejectionAlg.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"

using namespace std;  
using namespace anaobj;

/**
 *
 *Used to reject tracks from a SimpleTrackCollection.
 *Tracks that are far from the value of weighted 
 *average more than i * Weighted Average Error are
 *rejected and the Weighted Average recalculated.
 *
 *Define a crescent sort for vector<SimpleTrack>.
 *
 *Uses function wAverager, class SimpleTrack.
 *
 *
 *Authors M.De Mattia - R. Casagrande - 16/04/2008
 *
 */









  //definizione dell'operatore minore per il sort di una SimpleTrackCollection 
  bool simpleTkCollSort(const SimpleTrack & a, const SimpleTrack & b) {
    return fabs(a.z()) < fabs(b.z());
  }


  class RejectionAlg {
  private:
    SimpleTrackCollection temp_;
    double wAvgRecoTkError_;
    double wAvgRecoTkValue_;
    double left_;
    double right_;

  public:
    RejectionAlg(const SimpleTrackCollection & recoTk) {
      temp_.insert( temp_.end(), recoTk.begin(), recoTk.end() );
      sort( temp_.begin(), temp_.end(), simpleTkCollSort );
      wAvgRecoTkError_ = 0.;
      wAvgRecoTkValue_ = 0.;
      left_ = 0.;
      right_ = 0.;
    }

    SimpleTrackCollection eval(int i){

      if(temp_.size()>2){
      
	pair< int, pair< double, double > > wAvgRecoTk( wAverager(temp_));
	wAvgRecoTkValue_ = wAvgRecoTk.second.first; //coordinata z media delle tracce
	wAvgRecoTkError_ = wAvgRecoTk.second.second; //errore sulla media pesata
      
	left_  = fabs(temp_.front().z() - wAvgRecoTkValue_);
	right_ = fabs(temp_.back().z() - wAvgRecoTkValue_);
      
	if ( left_ > right_ && left_ > i*wAvgRecoTkError_) {
	  // remove left and iterate
	  temp_.erase(temp_.begin());
	  //
	  return eval(i);
	}
	else if ( right_ > i*wAvgRecoTkError_) {
	  // remove right and iterate
	  temp_.pop_back(); 
	  //
	  return eval(i);
	}
      }
      return temp_;
    }
  };


#endif


