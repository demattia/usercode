#ifndef REJECTIONALG_CC
#define REJECTIONALG_CC

#include "AnalysisExamples/AnalysisClasses/interface/RejectionAlg.h"

#include <vector>
#include <algorithm>

#include "AnalysisExamples/AnalysisClasses/interface/WAverager.h"
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



#endif


