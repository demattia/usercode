//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 19 13:21:45 2012 by ROOT version 5.32/00
// from TTree outputTree/outputTree
// found on file: histograms.root
//////////////////////////////////////////////////////////

#ifndef treeAnalyzer_h
#define treeAnalyzer_h

#include "treeAnalyzerBase.h"
#include "commonTools.C"

class treeAnalyzer : public treeAnalyzerBase {
  public :

    treeAnalyzer(TString fileName = "", const double & weight = 1., const bool electrons = false);
    ~treeAnalyzer();
        
   void     Loop();
   
    // Analysis cuts
    bool acceptanceCuts( std::vector<TreeCandidate>::const_iterator cand );
    bool trackSelectionCuts( std::vector<TreeCandidate>::const_iterator cand, const bool removeIsolationCut = false );
    bool lifetimeRelatedCuts( std::vector<TreeCandidate>::const_iterator candconst, bool removeLifetimeRelatedCuts = false, const bool decayLengthSigniInverted = false);
    bool dileptonSelectionCuts( std::vector<TreeCandidate>::const_iterator cand);
    bool analysisCuts( std::vector<TreeCandidate>::const_iterator cand, const bool removeIsolationCut = false, const bool removeLifetimeRelatedCuts = false, const bool decayLengthSigniInverted = false);

    bool triggerMatching( std::vector<TreeCandidate>::const_iterator cand );
    void initializeCuts();
    double ptCut_;
    double d0SignificanceCut_;
    double dPhiCorrCut_;
    double decayLengthSignificance2DCut_;
};

#endif

#ifdef treeAnalyzer_cxx
treeAnalyzer::treeAnalyzer(TString fileName, const double & weight, const bool electrons) :
    treeAnalyzerBase(fileName, weight, electrons)
{

}

treeAnalyzer::~treeAnalyzer()
{
}
#endif
