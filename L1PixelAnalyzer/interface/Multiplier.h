#ifndef MULTIPLIER_HH
#define MULTIPLIER_HH

/**
 * Takes the OfflineJet and L1Jet collections and generates new collections varying the 
 * et of the jets according to the resolution functions produced by JetResMaker. <br>
 *
 * The constructor takes the resolution files produced by JetResMaker and loads the histograms. <br>
 *
 * The initialize method builds the collections of ModL1Jets and ModOfflineJets which are used
 * for the multiplication. These are the copy of the initial collections and are made of ModJets. <br>
 * It also stores the original values of sMEt, MEt and MEtPhi, so it is not needed to use a new object
 * also for the MEt. The type given back is still an OfflineMEt object.
 *
 * The l1Jet collections must be passed in the following order:
 * central, tau and forward.
 *
 * The fillGen() method receives the GenJet collection (if exists). This method must be called
 * only if the multiplication is needed and only once per event. <br>
 * The same method makes the associations (in the same way as JetResMaker) and writes the refEt values
 * in the corresponding collections, only for the associated jets (the remaining jets will have refEt = 0). <br>
 *
 * The generate() method varies the et of the ModL1Jet and ModOfflineJet with refEt != 0, using
 * the histograms produced by JetResMaker.
 *
 * Author M. De Mattia - 6/1/2008
 */

#include <algorithm>
// L1Jet collection
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
// OfflineJet collection
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
// ModJet collection used for the multiplication
#include "AnalysisExamples/AnalysisObjects/interface/ModOfflineJet.h"

#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"
#include "AnalysisExamples/AnalysisClasses/interface/EtSort.h"

#include "TFile.h"
#include "TH1F.h"
//#include "TF1.h"

using namespace std;
using namespace anaobj;
using namespace reco;

class Multiplier {
public:
  /// Loads the files with the resolution histograms (once)
  Multiplier( const double & DR, const double & ETMIN, const double & ALPHA );
  /// Loads the collection of OfflineJets and OfflineMEt (once per event)
  vector<ModOfflineJet> * initialize( const OfflineJetCollection & OFFLINEJETS, const OfflineMEt & OFFLINEMET );
  /// Returns a pointer to the internal OfflineMEt object. Should be called only after initialize.
  OfflineMEt * offlineMEt() {
    return (& offlineMEt_);
  }
  /// Gets the GenJetCollection and makes the associations
  void fillGen( const GenJetCollection & GENJETS );
  /// Change the alpha factor
  void setAlpha( const double & ALPHA ) {
    alpha_ = ALPHA;
  }
  /// modifyes the et of associated offline jets and changes offline MEt
  void generate();
  /// Returns the number of jets that have been changed
  int numChanged() const {
    return numChanged_;
  }

private:
  double dR_;
  double etMin_;
  double alpha_;
  // To store the original MEt values
  double origMEt_;
  double origPhi_;
  double origSumEt_;

  vector<ModOfflineJet> modOfflineJets_;
  OfflineMEt offlineMEt_;
  TFile * JetRes3_histograms_;
  TH1F * Ptres_[40];
  char hname_[9];

  int numChanged_;
};

#endif // MULTIPLIER_HH
