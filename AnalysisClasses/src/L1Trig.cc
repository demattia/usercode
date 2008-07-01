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

// Requires three vectors of Jet objects and a float for MEt
void L1Trig::Fill (const BaseJetCollection & vec_L1CenJet, const BaseJetCollection & vec_L1ForJet, const BaseJetCollection & vec_L1TauJet, const double & L1MEt) {

  response_ = false;
  responseNoFor_ = false;
  responseMEtJet_ = false;
  responseMEtJetNoFor_ = false;

  vector<SimpleJet> vec_TriggerCenJet;
  vector<SimpleJet> vec_TriggerForJet;
  vector<SimpleJet> vec_TriggerTauJet;
  for ( BaseJetCollection::const_iterator tcj = vec_L1CenJet.begin(); tcj != vec_L1CenJet.end(); ++tcj ) {
    vec_TriggerCenJet.push_back( SimpleJet( tcj->et(), tcj->eta(), tcj->phi() ) );
  }
  for ( BaseJetCollection::const_iterator tfj = vec_L1ForJet.begin(); tfj != vec_L1ForJet.end(); ++tfj ) {
    vec_TriggerForJet.push_back( SimpleJet( tfj->et(), tfj->eta(), tfj->phi() ) );
  }
  for ( BaseJetCollection::const_iterator ttj = vec_L1TauJet.begin(); ttj != vec_L1TauJet.end(); ++ttj ) {
    vec_TriggerTauJet.push_back( SimpleJet( ttj->et(), ttj->eta(), ttj->phi() ) );
  }

  // Multijet
  // --------

  // Central
  // -------
  bool response_cen = false;
  Fill( vec_TriggerCenJet );
  response_cen = Response();  
  // Forward
  // -------
  bool response_for = false;
  Fill( vec_TriggerForJet );
  response_for = Response();
  // Tau
  // ---
  bool response_tau = false;
  Fill( vec_TriggerTauJet );
  response_tau = Response();

  // Full and no-forward
  // -------------------
  response_ = ( response_cen || response_tau || response_for );
  responseNoFor_ = ( response_cen || response_tau );

  // MEt + Jet
  // ---------
  if ( L1MEt > 0. ) {
    // Central
    // -------
    bool response_MEtJet_cen = false;
    sort( vec_TriggerCenJet.rbegin(), vec_TriggerCenJet.rend() );
    if ( vec_TriggerCenJet.size() != 0 ) {
      if ( (vec_TriggerCenJet[0].pt() >= 80.) && (L1MEt >= 100.) ) {
        response_MEtJet_cen = true;
      }
    }
    // Tau
    // ---
    bool response_MEtJet_tau = false;
    sort( vec_TriggerTauJet.rbegin(), vec_TriggerTauJet.rend() );
    if ( vec_TriggerTauJet.size() != 0 ) {
      if ( (vec_TriggerTauJet[0].pt() >= 80.) && (L1MEt >= 100.) ) {
        response_MEtJet_tau = true;
      }
    }
    // Forward
    // -------
    bool response_MEtJet_for = false;
    sort( vec_TriggerForJet.rbegin(), vec_TriggerForJet.rend() );
    if ( vec_TriggerForJet.size() != 0 ) {
      if ( (vec_TriggerForJet[0].pt() >= 80.) && (L1MEt >= 100.) ) {
        response_MEtJet_for = true;
      }
    }
    // Full and no-forward
    // -------------------
    responseMEtJet_ = ( response_MEtJet_cen || response_MEtJet_tau || response_MEtJet_for );
    responseMEtJetNoFor_ = ( response_MEtJet_cen || response_MEtJet_tau );
  }
}

pair<bool,bool> L1Trig::Response(bool useForwardL1Jets) {
  if (useForwardL1Jets) return (make_pair(response_, responseMEtJet_));
  else return (make_pair(responseNoFor_, responseMEtJetNoFor_));
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
