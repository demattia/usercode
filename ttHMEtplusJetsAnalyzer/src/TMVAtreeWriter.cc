#ifndef TMVATREEWRITER_CC
#define TMVATREEWRITER_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/TMVAtreeWriter.h"

TMVAtreeWriter::TMVAtreeWriter(const vector<TString> & variablesNames) {

  tmvaRootFile_ = new TFile("tmva.root", "RECREATE");
  tmvaRootFile_->cd();

  rootTree_ = new TTree("tree", "tmva tree");

  // This is used to store the per-event weights
  rootTree_->Branch("__WEIGHT__", &weight_, "__WEIGHT__/D");

  numberOfVariables_ = variablesNames.size();
  double * varPtr = new double[numberOfVariables_];
  variablesPtr_.reset(varPtr);

  vector<TString>::const_iterator name = variablesNames.begin();
  int var = 0;
  for( ; name != variablesNames.end(); ++name, ++var ) {
    rootTree_->Branch(*name, &(varPtr[var]), *name+"/D");
  }
}

void TMVAtreeWriter::fill(const vector<double> & variables, const double & weight) {
  vector<double>::const_iterator variable = variables.begin();
  int var = 0;
  // Write all the vector variables in the array
  // whose memory addresses are set for the TTree.
  double * varPtr = variablesPtr_.get();
  for( ; variable != variables.end(); ++variable, ++var ) {
    varPtr[var] = *variable;
  }

  // Set the weight for this event
  weight_ = weight;

  rootTree_->Fill();
}

TMVAtreeWriter::~TMVAtreeWriter() {

  tmvaRootFile_->cd();
  rootTree_->Write();
  tmvaRootFile_->Close();

}

#endif // TMVATREEWRITER_CC
