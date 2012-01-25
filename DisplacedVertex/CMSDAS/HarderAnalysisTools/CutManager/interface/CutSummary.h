#ifndef CutSummary_h
#define CutSummary_h 1

#include "HarderAnalysisTools/CutManager/interface/CutFlow.h"
#include "HarderAnalysisTools/Histograms/interface/HistMap.h"

#include "TTree.h"

#include <string>
#include <vector>
#include <map>

class CutSummary {
  
 public:
 
  CutSummary(std::string title);
  ~CutSummary() {};

  void addEntry(CutFlow& cuts,
		std::map<std::string,double>& weightMap,
		const int numDecays=0,
		const double ctau1=-1,
		const double ctau2=-1,
		const double leptonD01=-1,
		const double leptonD02=-1,
		const double mass=-1);

 private:

  unsigned numCutNames_;
  std::map<std::string,unsigned> cutNameMap_;

  TTree* bigTree_;
  std::vector<Float_t> bigTreeValues_;
  std::vector<Float_t> bigTreeBools_; // cannot use actual bools here due to some implicit unpacking that messes up variable addresses for ROOT tree
  std::vector<std::string> bigTreeNames_;

  HistMap* histos_;

  typedef struct {
    Float_t mass;
    Float_t weight;
    Float_t weight_up;
    Float_t weight_down;
    Float_t numDecays;
    Float_t ctau1;
    Float_t ctau2;
    Float_t leptonD01;
    Float_t leptonD02;
    bool passesAllCuts;
    bool passesAllCutsIgnoreLifetime;
  } bigtree_type;
  bigtree_type bigTreeVars_;

  std::vector<TTree*> trees_;
  typedef struct {
    Float_t value;
    Float_t mass;
    Float_t weight;
    Float_t weight_up;
    Float_t weight_down;
    Float_t numDecays;
    Float_t ctau1;
    Float_t ctau2;
    Float_t leptonD01;
    Float_t leptonD02;
    bool passesThisCut;
    bool passesAllOtherCuts;
    bool passesAllOtherCutsIgnoreLifetime;
  } tree_type;
  std::vector<tree_type*> treeVariables_;

} ;

#endif
