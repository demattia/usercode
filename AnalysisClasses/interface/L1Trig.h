#ifndef L1TRIG_H
#define L1TRIG_H
//#define DEBUG
/**
 *
 * This class gives the L1 trigger response. It requires a vector of Jet
 * objects (Jet objects are included) and can also get the L1MET.
 * the L1 response is given by 
 *
 */
#include <vector>
#include <algorithm>
#include "SimpleJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
// used to pass calibrated jets to this class, see Init_L1Trig

using namespace std;
using namespace anaobj;

class L1Trig {
  float _calib_jet1_et;
  float _calib_jet2_et;
  float _calib_jet3_et;
  float _calib_jet4_et;
  float _L1_met;
  bool et1_cut_;
  bool et2_cut_;
  bool et3_cut_;
  bool et4_cut_;
  bool response_;
  bool responseNoFor_;
  bool responseMEtJet_;
  bool responseMEtJetNoFor_;

 public:
  L1Trig () {
    _calib_jet1_et = 0;
    _calib_jet2_et = 0;
    _calib_jet3_et = 0;
    _calib_jet4_et = 0;
    _L1_met = 0;
    et1_cut_ = false;
    et2_cut_ = false;
    et3_cut_ = false;
    et4_cut_ = false;
    response_ = false;
    responseNoFor_ = false;
    responseMEtJet_ = false;
    responseMEtJetNoFor_ = false;

  }

  /// Function to sort the vector in Pt
  void Sort (vector<SimpleJet> &);
  /// Requires a vector of Jet objects and a float for the MEt
  void Fill (vector<SimpleJet> &, float L1MEt = 0);
  /// Requires three vectors of Jet objects and a bool default true asking if the forward jets must be used
  void Fill (const BaseJetCollection &, const BaseJetCollection &, const BaseJetCollection &, const double & MEt = 0);
  /** Get the response of the trigger. The bool is required: true = use L1FJets
   * false = do not use L1FJets.
   * Returns a pair(multijet response, MEt+Jet response).
   */
  pair<bool,bool> Response (bool useForwardL1Jets);
  /// Get the standard L1 response for the Multijet Trigger
  bool Response ();
  /// Response to the L1 Multijet trigger for the assigned cuts
  bool Response (float, float, float, float);
  /// Response to the L1 MET+Jet trigger for the assigned cuts (this order)
  bool Response (float, float);
  bool Et1Cut () {
    return et1_cut_;
  }
  bool Et2Cut () {
    return et2_cut_;
  }
  bool Et3Cut () {
    return et3_cut_;
  }
  bool Et4Cut () {
    return et4_cut_;
  }
};

#endif // L1TRIG_H
