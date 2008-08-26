#ifndef WAVERAGER_H
#define WAVERAGER_H

#include <vector>
#include <algorithm>

#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"

using namespace std;

/**
 * Used to make a weighted average for a simpletrackcollection z value
 * 
 *
 *  Author R.Casagrande - 15/04/2008
 *
 */

pair< int, pair< double, double > > wAverager(const anaobj::SimpleTrackCollection & vecTk);  


#endif
