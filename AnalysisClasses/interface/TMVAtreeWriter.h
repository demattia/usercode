#ifndef TMVATREEWRITER_H
#define TMVATREEWRITER_H

/**
 * Template class, with parameter T, used to write the variables in a root tree for TMVA.
 * The constructor requires a vector<TString> with the names of the variables.
 * The size of this vector is used to create an array for the variables, which are all T.
 * It optionally can get a suffix which it appends to
 * as tmva_suffix.root for the output file name.
 * The fill method gets a vector<T> for the variables to store. The recommended type is Float_t
 * for compatibility with TMVA.
 * It fills the array and then fills the TTree. The branches are of type Float_t (option /F).
 * This two vectors must be ordered accordingly to maintain correspondence name-variable.
 * The size of the vectors is not checked. 
 * The loop is made on the variables, 
 * a shorter vector of names would crash the program.
 * ATTENTION:
 * the methods writing to file are invoked in the destructor. 
 * If an object of this class is created with new, it must be deleted,
 * otherwise the tree will never be written. In this case an auto_ptr could be used.
 */

#include <vector>
#include <memory>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"

using namespace std;

template <class T>
class TMVAtreeWriter {

public:
  TMVAtreeWriter(const vector<TString> & variablesNames, const TString & suffix = "");
  ~TMVAtreeWriter();

  /// Added fot better compatibility with TMVA
  void fill(const vector<T> & variables, const T & weight = 1.);
private:
  TFile * tmvaRootFile_;
  TTree * rootTree_;

  unsigned int numberOfVariables_;
  auto_ptr<T> variablesPtr_;
  T weight_;
};

template <class T>
TMVAtreeWriter<T>::TMVAtreeWriter(const vector<TString> & variablesNames, const TString & suffix) {

  tmvaRootFile_ = new TFile("tmva"+suffix+".root", "RECREATE");
  tmvaRootFile_->cd();

  TString tempSuffix = "_";
  tempSuffix += suffix;
  rootTree_ = new TTree("tree"+tempSuffix, "tmva tree "+suffix);

  // This is used to store the per-event weights (Float_t branch)
  rootTree_->Branch("__WEIGHT__", &weight_, "__WEIGHT__/F");

  numberOfVariables_ = variablesNames.size();
  T * varPtr = new T[numberOfVariables_];
  variablesPtr_.reset(varPtr);

  vector<TString>::const_iterator name_itr = variablesNames.begin();
  int var = 0;
  for( ; name_itr != variablesNames.end(); ++name_itr, ++var ) {
    // Float_t branch
    rootTree_->Branch(*name_itr, &(varPtr[var]), *name_itr+"/F");
  }
}

template<class T>
void TMVAtreeWriter<T>::fill(const vector<T> & variables, const T & weight) {
  typename vector<T>::const_iterator variable_itr = variables.begin();
  int var = 0;
  // Write all the vector variables in the array
  // whose memory addresses are set for the TTree.
  T * varPtr = variablesPtr_.get();
  for( ; variable_itr != variables.end(); ++variable_itr, ++var ) {
    varPtr[var] = *variable_itr;
  }

  // Set the weight for this event
  weight_ = weight;

  rootTree_->Fill();
}

template<class T>
TMVAtreeWriter<T>::~TMVAtreeWriter() {

  tmvaRootFile_->cd();
  rootTree_->Write();
  tmvaRootFile_->Close();

}
#endif // TMVATREEWRITER_H
