#ifndef L1TRIG_H
#define L1TRIG_H
//#define DEBUG
////////////////////////////////////////////////////////////////////////////
//
// This class gives the L1 trigger response. It requires a vector of Jet
// objects (Jet objects are included) and can also get the L1MET.
// the L1 response is given by 
//
////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <algorithm>
#include "SimpleJet.h"
// used to pass calibrated jets to this class, see Init_L1Trig

using namespace std;

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
  }

  // Function to sort the vector in Pt
  void Sort (vector<SimpleJet> &);
  // Requires a vector of Jet objects and a float for the MET
  void Fill (vector<SimpleJet> &, float);
  // alternate method requiring only the vector of jet, MET is set = 0
  void Fill (vector<SimpleJet> & );
  // Get the standard L1 response for the Multijet Trigger
  bool Response ();
  // Response to the L1 Multijet trigger for the assigned cuts
  bool Response (float, float, float, float);
  // Response to the L1 MET+Jet trigger for the assigned cuts (this order)
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
