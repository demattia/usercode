#ifndef L1TRIG_CC
#define L1TRIG_CC

#include "AnalysisExamples/AnalysisClasses/interface/L1Trig.h"

void L1Trig::Sort ( vector<SimpleJet> & vec_L1jet ) {
  // Overload the < operator to sort Jet objects
  // -------------------------------------------

  sort( vec_L1jet.begin(), vec_L1jet.end() );
  reverse( vec_L1jet.begin(), vec_L1jet.end() );
}

void L1Trig::Fill (vector<SimpleJet> & vec_L1jet, float L1_met) {
  // sort the vector in Pt
  Sort(vec_L1jet);

  _calib_jet1_et = 0;
  _calib_jet2_et = 0;
  _calib_jet3_et = 0;
  _calib_jet4_et = 0;
  _L1_met = 0;

  if (vec_L1jet.size() > 0) {
    _calib_jet1_et = vec_L1jet.at(0).pt();
  }
  if (vec_L1jet.size() > 1) {
    _calib_jet2_et = vec_L1jet.at(1).pt();
  }
  if (vec_L1jet.size() > 2) {
    _calib_jet3_et = vec_L1jet.at(2).pt();
  }
  if (vec_L1jet.size() > 3) {
    _calib_jet4_et = vec_L1jet.at(3).pt();
  }
  _L1_met = L1_met;
}
void L1Trig::Fill (vector<SimpleJet> & vec_L1jet) {
  // sort the vector in Pt
  Sort(vec_L1jet);

  _calib_jet1_et = 0;
  _calib_jet2_et = 0;
  _calib_jet3_et = 0;
  _calib_jet4_et = 0;
  _L1_met = 0;

  if (vec_L1jet.size() > 0) {
    _calib_jet1_et = vec_L1jet.at(0).pt();
  }
  if (vec_L1jet.size() > 1) {
    _calib_jet2_et = vec_L1jet.at(1).pt();
  }
  if (vec_L1jet.size() > 2) {
    _calib_jet3_et = vec_L1jet.at(2).pt();
  }
  if (vec_L1jet.size() > 3) {
    _calib_jet4_et = vec_L1jet.at(3).pt();
  }
}

bool L1Trig::Response () {
#ifdef DEBUG
  std::cout << "_calib_jet1_et = " << _calib_jet1_et << std::endl;
  std::cout << "_calib_jet2_et = " << _calib_jet2_et << std::endl;
  std::cout << "_calib_jet3_et = " << _calib_jet3_et << std::endl;
  std::cout << "_calib_jet4_et = " << _calib_jet4_et << std::endl;
#endif
  et1_cut_ = false;
  if (_calib_jet1_et >= 250.) et1_cut_ = true;
  et2_cut_ = false;
  if (_calib_jet2_et >= 200.) et2_cut_ = true;
  et3_cut_ = false;
  if (_calib_jet3_et >= 100.) et3_cut_ = true;
  et4_cut_ = false;
  if (_calib_jet4_et >= 80.) et4_cut_ = true;

  // High luminosity thresholds (taken from the trigger TDR)
  if (_calib_jet1_et >= 250. || _calib_jet2_et >= 200. || _calib_jet3_et >= 100. || _calib_jet4_et >= 80.)
    return true;
  return false;
}
bool L1Trig::Response (float ET1, float ET2, float ET3, float ET4) {
  if (_calib_jet1_et >= ET1 || _calib_jet3_et >= ET2 || _calib_jet3_et >= ET3 || _calib_jet4_et >= ET4)
    return true;
  return false;
}
bool L1Trig::Response (float ET1, float L1MET) {
  if (_calib_jet1_et >= ET1 && _L1_met >= L1MET)
    return true;
  return false;
}

#endif // L1TRIG_CC
