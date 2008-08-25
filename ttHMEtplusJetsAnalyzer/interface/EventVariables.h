#ifndef EVENTVARIABLES_H
#define EVENTVARIABLES_H

/**
 * Evaluates variables that characterize the event.
 * This variables are stored in histograms which names are
 * built appending a suffix passed as argument to the constructor.
 * It can optionally write root files in TMVA format with all
 * the variables.
 *
 * Author: M. De Mattia
 * Date: 7/7/2008
 */

#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/Particle.h"
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace anaobj;

class EventVariables {
public:
  EventVariables( const string & higgsFileName, const string & hadronicTopFileName, const string & qcdFileName, TString suffix, TFile * outputFile, bool fillHistograms = true );
  ~EventVariables();
  /// Used to pass the collections. Takes the jetCollection by value since it modifies it removing the jets associated to the Higgs.
  vector<double> fill( vector<const OfflineJet *> jetCollection, const vector<const OfflineJet *> & bTaggedJetCollection, const OfflineMEt * offlineMEt );

private:

  /** Used to fill the matrices for the probability of Higgs jet pairs, top jet triplets and also the b-tag probability matrix for qcd
   * takes by reference the *** because it creates them with new inside and assign the pointers to new values
   */
  void fillProbabilityMatrices(const string & probabilityFileName, unsigned int * binNum, double * binSize, unsigned int ***& trueArray, unsigned int ***& falseArray);
  /// Used to evaluate a vector of the first N jets in the event
  Particle<const OfflineJet> firstNjetsParticle( const vector<const OfflineJet *> & jetCollection, const int N );
  /// Used to evaluate variables on all the selected jets (Ht, SumEt, ...)
  void allGoodJetsVariables( const vector<const OfflineJet *> & offlineJets, const OfflineMEt * offlineMEt );
  /// Used to evaluate the ratios and select the combinations of jets for the Higgs candidates
  double evalHiggsPairProbability(const Particle<const OfflineJet> & higgsCandidate) const;
  /// Used to evaluate the ratios and select the combinations of jets for the hadronic top candidates
  double evalTopTripletProbability(const Particle<const OfflineJet> & hadronicTopCandidate, const Particle<const OfflineJet> & selectedHiggs) const;
  /// Used to evaluate the hadronic W from the selected Hadronic Top
  Particle<const OfflineJet> getWfromHadronicTop(const Particle<const OfflineJet> & selectedHadronicTop ) const;

  bool fillHistograms_;

  // This will be the multidimensional array
  unsigned int *** trueH_;
  unsigned int *** falseH_;
  unsigned int *** trueHadronicTop_;
  unsigned int *** falseHadronicTop_;
  unsigned int *** taggedJet_;
  unsigned int *** notTaggedJet_;

  unsigned int higgsBinNum_[3];
  double higgsBinSize_[3];
  unsigned int hadronicTopBinNum_[3];
  double hadronicTopBinSize_[3];
  unsigned int qcdBinNum_[3];
  double qcdBinSize_[3];

  // Reference masses for Higgs, Top and W
  double referenceHiggsMass_;
  double referenceTopMass_;
  double referenceWmass_;

  vector<double> eventVariablesVector_;

  // Histograms
  TH1D * higgsMass_;
  TH1D * hadronicTopMass_;
  TH1D * hadronicWmass_;
  TH1D * chi2ofMasses_;
  TH1D * firstNjetsMass_[2];
  TH1D * firstNjetsCentrality_[2];
  TH1D * hadronicTopProjectionAlongHiggsDirection_;
  TH1D * deltaEtaHadronicTopHiggs_;
  TH1D * goodHt_;
  TH1D * mEtSig_;
  TH1D * deltaPhiMEtNthLeadingJet_[3];
  TH1D * hadronicTopPlusHiggsMass_;
  TH1D * sumHighEffDiscriminantFirst4Jets_;
  TH1D * sumHighEffDiscriminantFirst6Jets_;
  TH1D * remainingJetsMass_;
  TH1D * sixthJetEt_;
  TH1D * bTagTkInvMass_;

  TFile * outputFile_;
  TDirectory * outputDir_;
};

#endif // EVENTVARIABLES_H
