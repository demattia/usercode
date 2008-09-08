/**
 * 
 */


#ifndef JETVERTEXASSOCIATOR_H
#define JETVERTEXASSOCIATOR_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"

//Implemented functions
//-----------------------
#include "AnalysisExamples/AnalysisClasses/interface/SortingDef.h"
#include "AnalysisExamples/AnalysisClasses/interface/DzAssociator.h"
#include "AnalysisExamples/AnalysisClasses/interface/WAverager.h"

//Implemented Classes
//------------------------
#include "AnalysisExamples/AnalysisClasses/interface/RejectionAlg.h"
#include "AnalysisExamples/AnalysisClasses/interface/VerticesJetCounter.h"
#include "AnalysisExamples/AnalysisClasses/interface/ProbCal.h"
#include "AnalysisExamples/AnalysisClasses/interface/ProbEval.h"

using namespace std;  
using namespace anaobj;

class JetVertexAssociator {

  std::vector<double> pMatrixSing25_;
  std::vector<double> pMatrix25_;
  
  double jetEtCut_; 
  double  jetEtaCut_;
  int sigmaCut_;
  double sigmaVtx_;

 public:
  JetVertexAssociator( double jetEtCut, double jetEtaCut, int sigmaCut = 13, double SigmaVtx = 0.02 );

  OfflineJetCollection associate(const OfflineJetCollection & caloJets, const BaseVertexCollection & recoVertexes);


};


#endif


