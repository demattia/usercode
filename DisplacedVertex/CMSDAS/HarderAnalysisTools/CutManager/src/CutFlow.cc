#include <iostream>
#include <cstdlib>
#include "HarderAnalysisTools/CutManager/interface/CutFlow.h"

CutFlow::CutFlow() : passAll_(true), passAllIgnoreLifetime_(true) {

  cutResults_.clear();
}


void CutFlow::applyCut(std::string cutName, double value, double minValue, double maxValue) {

  if (values_.find(cutName)!=values_.end()) {
    std::cout << "CutFlow ERROR: double use of cut name " << cutName << std::endl;
    std::exit(100);
  }
  values_[cutName]=value;
  bool passThisCut=(value>=minValue && value<=maxValue);
  cutResults_[cutName]=passThisCut;
  passAll_&=passThisCut;
  passAllIgnoreLifetime_&=passThisCut;
}


void CutFlow::applyCut(std::string cutName, bool passThisCut) {

  if (values_.find(cutName)!=values_.end()) {
    std::cout << "CutFlow ERROR: double use of cut name " << cutName << std::endl;
    std::exit(100);
  }
  values_[cutName]=double(passThisCut);
  cutResults_[cutName]=passThisCut;
  passAll_&=passThisCut;
  passAllIgnoreLifetime_&=passThisCut;
}


void CutFlow::applyLifetimeCut(std::string cutName, double value, double minValue, double maxValue) {

  if (values_.find(cutName)!=values_.end()) {
    std::cout << "CutFlow ERROR: double use of cut name " << cutName << std::endl;
    std::exit(100);
  }
  values_[cutName]=value;
  bool passThisCut=(value>=minValue && value<=maxValue);
  lifetimeCutResults_[cutName]=passThisCut;
  passAll_&=passThisCut;
}


void CutFlow::applyLifetimeCut(std::string cutName, bool passThisCut) {

  if (values_.find(cutName)!=values_.end()) {
    std::cout << "CutFlow ERROR: double use of cut name " << cutName << std::endl;
    std::exit(100);
  }
  values_[cutName]=double(passThisCut);
  cutResults_[cutName]=passThisCut;
  passAll_&=passThisCut;
}


bool CutFlow::passCut(std::string cutName) {

  if (cutResults_.find(cutName)!=cutResults_.end()) {
    return cutResults_[cutName];
  } else if (lifetimeCutResults_.find(cutName)!=lifetimeCutResults_.end()) {
    return lifetimeCutResults_[cutName];
  } else {
    std::cout << "CutFlow ERROR: cut name " << cutName << " not found" << std::endl;
    std::exit(100);
  }
}


double CutFlow::getValue(std::string cutName) {
  if (values_.find(cutName)!=values_.end()) {
    return values_[cutName];
  } else {
    std::cout << "CutFlow ERROR: cut name " << cutName << " not found" << std::endl;
    std::exit(100);
  }
}


bool CutFlow::passAll() {

  return passAll_;
}


bool CutFlow::passAllIgnoreLifetime() {

  return passAllIgnoreLifetime_;
}


bool CutFlow::passAllFailOne(std::string cutName) {

  if (values_.find(cutName)==values_.end()) {
    std::cout << "CutFlow ERROR: cut name " << cutName << " not found" << std::endl;
    std::exit(100);
  }

  // simple case: all cuts passed anyway
  if (passAll_) return false;

  // difficult case: at least one cut failed. check whether it is the requested one.
  bool passAllOthers=true;
  for (std::map<std::string,bool>::const_iterator it=cutResults_.begin();
       it!=cutResults_.end(); it++) {
    if (it->first!=cutName) {
      passAllOthers&=it->second;
    }
  }
  for (std::map<std::string,bool>::const_iterator it=lifetimeCutResults_.begin();
       it!=lifetimeCutResults_.end(); it++) {
    if (it->first!=cutName) {
      passAllOthers&=it->second;
    }
  }
  return passAllOthers;
}


bool CutFlow::passAllIgnoreOne(std::string cutName) {

  if (values_.find(cutName)==values_.end()) {
    std::cout << "CutFlow ERROR: cut name " << cutName << " not found" << std::endl;
    std::exit(100);
  }

  // simple case: all cuts passed anyway
  if (passAll_) return true;

  // difficult case: at least one cut failed. check whether it is the requested one.
  return passAllFailOne(cutName);
}


bool CutFlow::passAllIgnoreOneIgnoreLifetime(std::string cutName) {

  if (values_.find(cutName)==values_.end()) {
    std::cout << "CutFlow ERROR: cut name " << cutName << " not found" << std::endl;
    std::exit(100);
  }

  // simple case: all cuts passed anyway
  if (passAll_) return true;

  // difficult case: at least one cut failed. check whether it is the requested one.
  // loop only over regular cuts, not over lifetime cuts
  bool passAllOthers=true;
  for (std::map<std::string,bool>::const_iterator it=cutResults_.begin();
       it!=cutResults_.end(); it++) {
    if (it->first!=cutName) {
      passAllOthers&=it->second;
    }
  }
  return passAllOthers;
}
