#ifndef TMVATREEWRITER_H
#define TMVATREEWRITER_H

/**
 * Class used to write the variables in a root tree for TMVA.
 * The constructor requires a vector<TString> with the names of the variables. The size of this vector is used
 * to create an array for the variables, which are all double. It optionally can get a suffix which it appends to
 * as tmva_suffix.root for the output file name.
 * The fill method gets a vector<double> for the variables to store. It fills the array and then fills the TTree.
 * This too vectors must be ordered accordingly to maintain correspondence name-variable.
 * The size of the vectors is not checked. The loop is made on the variables, a shorter vector of names
 * would crash the program.
 * ATTENTION:
 * the methods writing to file are invoked in the destructor. If an object of this class is created with new,
 * it must be deleted, otherwise the tree will never be written. In this case an auto_ptr could be used.
 */

#include <vector>
#include <memory>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"

using namespace std;

class TMVAtreeWriter {

public:
  TMVAtreeWriter(const vector<TString> & variablesNames, const TString & suffix = "");
  ~TMVAtreeWriter();

  void fill(const vector<double> & variables, const double & weight = 1.);
private:
  TFile * tmvaRootFile_;
  TTree * rootTree_;

  unsigned int numberOfVariables_;
  auto_ptr<double> variablesPtr_;
  double weight_;
};

#endif // TMVATREEWRITER_H
