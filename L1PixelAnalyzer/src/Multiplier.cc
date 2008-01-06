#ifndef MULTIPLIER_CC
#define MULTIPLEIR_CC

#include "AnalysisExamples/L1PixelAnalyzer/interface/Multiplier.h"

using namespace std;
using namespace anaobj;
using namespace reco;

Multiplier::Multiplier( const BaseJetCollection & CENJETS, const BaseJetCollection & TAUJETS,
                        const BaseJetCollection & FORJETS, const OfflineJetCollection & OFFLINEJETS ) {

  // Create the new ModJet collections holding the original values and the modified ones
  BaseJetCollection::const_iterator base_it = CENJETS.begin();
  for( ; base_it != CENJETS.end(); ++base_it ) {
    modCenJets.push_back( ModJet<BaseJet>(*base_it) );
  }
  for( base_it = TAUJETS.begin(); base_it != TAUJETS.end(); ++base_it ) {
    modTauJets.push_back( ModJet<BaseJet>(*base_it) );
  }
  for( base_it = FORJETS.begin(); base_it != FORJETS.end(); ++base_it ) {
    modForJets.push_back( ModJet<BaseJet>(*base_it) );
  }
  OfflineJetCollection::const_iterator off_it = OFFLINEJETS.begin();
  for( ; off_it != OFFLINEJETS.end(); ++off_it ) {
    modOfflineJets.push_back( ModJet<OfflineJet>(*off_it) );
  }
}

void Multiplier::fillGen( const GenJetCollection & GENJETS ) {
  // Create vectors of pointers to use in the association
  vector<const ModJet<BaseJet> *> modL1JetPtrs;
  vector<ModJet<BaseJet> >::const_iterator base_it = modCenJets.begin();
  for( ; base_it != modCenJets.end(); ++base_it ) {
    modL1JetPtrs.push_back( &(*base_it) );
  }
  for( base_it = modTauJets.begin(); base_it != modTauJets.end(); ++base_it ) {
    modL1JetPtrs.push_back( &(*base_it) );
  }
  for( base_it = modForJets.begin(); base_it != modForJets.end(); ++base_it ) {
    modL1JetPtrs.push_back( &(*base_it) );
  }
  vector<const ModJet<OfflineJet> *> modOfflineJetPtrs;
  vector<ModJet<OfflineJet> >::const_iterator off_it = modOfflineJets.begin();
  for( ; off_it != modOfflineJets.end(); ++off_it ) {
    modOfflineJetPtrs.push_back( &(*off_it) );
  }
  vector<const GenJet *> genJetPtrs;
  GenJetCollection::const_iterator gen_it = GENJETS.begin();
  for( ; gen_it != GENJETS.end(); ++gen_it ) {
    genJetPtrs.push_back( &(*gen_it) );
  }

  // Sort all the collections
  sort( modL1JetPtrs.rbegin(), modL1JetPtrs.rend(), EtSort<ModJet<BaseJet> >() );
  sort( modOfflineJetPtrs.rbegin(), modOfflineJetPtrs.rend(), EtSort<ModJet<OfflineJet> >() );
  sort( genJetPtrs.rbegin(), genJetPtrs.rend(), EtSort<GenJet>() );

  const float drmin = 0.3;

  // Associate L1Jets with offlineJets
  AssociatorEt<ModJet<OfflineJet>, ModJet<BaseJet> > l1OffAssociator( drmin );
  auto_ptr<map<const ModJet<OfflineJet>*, const ModJet<BaseJet>*> > l1JetOfflineJetMap( l1OffAssociator.Associate( modOfflineJetPtrs, modL1JetPtrs ) );

  // Fill the vector of offlineJets associated to L1Jets
  vector<const ModJet<OfflineJet> *> assocModOfflineJetPtrs;
  map<const ModJet<OfflineJet>*, const ModJet<BaseJet>*>::const_iterator l1OffMap_it = l1JetOfflineJetMap->begin();
  for( ; l1OffMap_it != l1JetOfflineJetMap->end(); ++l1OffMap_it ) {
    assocModOfflineJetPtrs.push_back( l1OffMap_it->first );
  }
  sort( assocModOfflineJetPtrs.begin(), assocModOfflineJetPtrs.end(), EtSort<ModJet<OfflineJet> >() );

  // Associate this new collection with the goodGenJets
  AssociatorEt<ModJet<OfflineJet>, GenJet > genOffAssociator( drmin );
  auto_ptr<map<const ModJet<OfflineJet>*, const GenJet*> > genJetOfflineJetMap( genOffAssociator.Associate( assocModOfflineJetPtrs, genJetPtrs ) );

  // Loop on the map and write the refEt variables for the associated jets
  map<const ModJet<OfflineJet>*, const GenJet*>::const_iterator genOffMap_it = genJetOfflineJetMap->begin();
  for( ; genOffMap_it != genJetOfflineJetMap->end(); ++genOffMap_it ) {

    // Remove const from the pointer, needed to access non const method
    ModJet<OfflineJet> * nonConstOffPtr = const_cast<ModJet<OfflineJet>* >(genOffMap_it->first);
    nonConstOffPtr->setRefEt( genOffMap_it->second->et() );

    ModJet<BaseJet> * nonConstL1Ptr = const_cast<ModJet<BaseJet>* >( ( l1JetOfflineJetMap->find( genOffMap_it->first ) )->second );
    nonConstL1Ptr->setRefEt( genOffMap_it->first->et() );
  }
}

void generate() {

}

#endif // MULTIPLIER_CC
