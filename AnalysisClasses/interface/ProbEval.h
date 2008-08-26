#ifndef PROBEVAL_H
#define PROBEVAL_H

#include <vector>
#include <map>
#include <iostream>
#include "AnalysisExamples/AnalysisObjects/interface/SimpleCaloJet.h"
#include "AnalysisExamples/AnalysisClasses/interface/VerticesJetCounter.h"
#include "AnalysisExamples/AnalysisClasses/interface/ProbCal.h"
#include "AnalysisExamples/AnalysisClasses/interface/SortingDef.h"
#include "TMath.h"

/**
 * Evaluates the probability for a cofiguration of jets associated to the
 * verteces in the event. <BR>
 * It uses a recursive method which fills a member vector with all the
 * probabilities.
 *
 * 28/5/2008
 */

using namespace std;
using namespace anaobj;

class ProbEval {
 protected:
  //  vector<double> probVec_;
  int index_;
  SimpleCaloJet* jetArray_;
  int jetNum_;
  int * jetVtxId_;
  VerticesJetCounter vtxCounter_;
  ProbCal probCal_;
  int maxProbId_;
  double maxProb_;
  int jetNumCut_;

  // protected method called by evalProb
  void eval_(int iJet, double tempProb);
 public:
  ProbEval() {
    index_ = 0;
    jetArray_ = 0;
    jetNum_ = 0;
    jetVtxId_ = 0;
    maxProbId_ = -1;
    maxProb_ = -10000.;
    jetNumCut_ = 0;
  }
  //  ~ProbEval() {
  //  }

  /// receives the map jets-dzs
    map<int, int> evalProb(SimpleCaloJetCollection vec_CaloJet, int jetNumCut = 20 );
};


#endif
