#ifndef MULTIPLIER_HH
#define MULTIPLIER_HH

/**
 * Takes the OfflineJet and L1Jet collections and generates new collections varying the 
 * et of the jets according to the resolution functions produced by JetResMaker. <br>
 *
 * The constructor builds the collections of ModL1Jets and ModOfflineJets which are used
 * for the multiplication. These are the copy of the initial collections and are made of ModJets. <br>
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
// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
// ModJet collection used for the multiplication
#include "AnalysisExamples/AnalysisObjects/interface/ModJet.h"

#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"
#include "AnalysisExamples/AnalysisClasses/interface/EtSort.h"

using namespace std;
using namespace anaobj;
using namespace reco;

class Multiplier {
public:
  Multiplier( const BaseJetCollection & CENJETS, const BaseJetCollection & TAUJETS,
              const BaseJetCollection & FORJETS, const OfflineJetCollection & OFFLINEJETS );
  /// Gets the GenJetCollection and makes the associations
  void fillGen( const GenJetCollection & GENJETS );
  /// modifyes the et of associated offline and l1 jets and changes offline MEt
  void generate();

private:
  vector<ModJet<BaseJet> > modCenJets;
  vector<ModJet<BaseJet> > modTauJets;
  vector<ModJet<BaseJet> > modForJets;
  vector<ModJet<OfflineJet> > modOfflineJets;
};

#endif // MULTIPLIER_HH
