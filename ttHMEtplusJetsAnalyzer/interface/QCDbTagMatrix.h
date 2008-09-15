#ifndef QCDBTAGMATRIX_HH
#define QCDBTAGMATRIX_HH

/**
 * Class used for the multiplication of QCD with the probability matrix method.
 * It inherits from EventVariables.
 * the multiply method receives a vector<OfflineJets> (of at most 8 jets by default) and
 * it evaluates all possible combinations giving at least N tags. For each combination,
 * the fill method of EventVariables is called with the probability of that combination
 * as weight factor, effectively filling an event multiple times. Each jet has its discriminant
 * variable taken randomly from the histogram (b or no-b depending if in that combination it
 * is considered tagged or not).
 * The number of required b-jets (greatel-equal) is passed in the constructor.
 *
 * Author: M. De Mattia
 * Date: 12/9/2008
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/EventVariables.h"

class QCDbTagMatrix : public EventVariables {
public:
  /// This constructor calls the EventVariables constructor
  QCDbTagMatrix( const string & higgsFileName, const string & hadronicTopFileName, const string & qcdFileName, TString suffix, TFile * outputFile, bool fillHistograms, const string & qcdHistoFileName, int bJetNumCut );
  ~QCDbTagMatrix();
  void multiply( const vector<const OfflineJet *> jetCollection, const OfflineMEt * offlineMEt );
protected:
  int bJetNumCut_;
  TFile * inputFileSignal_;
  // TH1F * taggedJetEt_;
  // TH1F * taggedJetEta_;
  // TH1F * taggedJetS1_;
  TH1F * taggedJetTagMass_;
  TH1F * taggedJetDiscriminatorHighEff_;
  // TH1F * notTaggedJetEt_;
  // TH1F * notTaggedJetEta_;
  // TH1F * notTaggedJetS1_;
  TH1F * notTaggedJetTagMass_;
  TH1F * notTaggedJetDiscriminatorHighEff_;
};

#endif // QCDBTAGMATRIX_HH
