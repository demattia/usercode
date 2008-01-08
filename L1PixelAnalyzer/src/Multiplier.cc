#ifndef MULTIPLIER_CC
#define MULTIPLEIR_CC

#include "AnalysisExamples/L1PixelAnalyzer/interface/Multiplier.h"

using namespace std;
using namespace anaobj;
using namespace reco;

Multiplier::Multiplier( const double & DR, const double & ETMIN, const double & ALPHA ) {
  dR_ = DR;
  etMin_ = ETMIN;
  alpha_ = ALPHA;
  numChanged_ = 0;

  // Get resolution histograms from original files 
  // ---------------------------------------------
  JetRes3_histograms_ = new TFile("JetRes3.hist");
  
  for ( int i=0; i<40; ++i ) {
    sprintf ( hname_, "Ptres%d", i );
    Ptres_[i] = dynamic_cast<TH1F*>(JetRes3_histograms_->Get(hname_)); 
  }  
}

vector<ModOfflineJet> * Multiplier::initialize( const OfflineJetCollection & OFFLINEJETS, const OfflineMEt & OFFLINEMET ) {

  // Save the original values
  origMEt_ = OFFLINEMET.et();
  origPhi_ = OFFLINEMET.phi();
  origSumEt_ = OFFLINEMET.sumEt();
  //  origMEtSig_ = OFFLINEMET.mEtSig();
  // Copy the values in the local offlineMEt
  offlineMEt_.setEt(origMEt_);
  offlineMEt_.setPhi(origPhi_);
  offlineMEt_.setSumEt(origSumEt_);
  offlineMEt_.setMEtSig(OFFLINEMET.mEtSig());

  // Free the vector of offlineJets, if this is not the first event it still contains the jets from the previous event
  modOfflineJets_.clear();

  // Create the new ModOfflineJet collection holding the original values and the modified ones
  OfflineJetCollection::const_iterator off_it = OFFLINEJETS.begin();
  for( ; off_it != OFFLINEJETS.end(); ++off_it ) {
    modOfflineJets_.push_back( ModOfflineJet(*off_it) );
  }
  return (& modOfflineJets_);
}

void Multiplier::fillGen( const GenJetCollection & GENJETS ) {
  // Create vectors of pointers to use in the association
  // Use only offlineJets with Et > etMin_ GeV and eta < 2.5
  vector<const ModOfflineJet*> modOfflineJetPtrs;
  vector<ModOfflineJet>::const_iterator off_it = modOfflineJets_.begin();
  for( ; off_it != modOfflineJets_.end(); ++off_it ) {
    if ( off_it->et() > etMin_ && fabs(off_it->eta()) < 2.5 ) {
      modOfflineJetPtrs.push_back( &(*off_it) );
    }
  }
  vector<const GenJet *> genJetPtrs;
  GenJetCollection::const_iterator gen_it = GENJETS.begin();
  for( ; gen_it != GENJETS.end(); ++gen_it ) {
    if ( gen_it->et() > 20 ) {
      genJetPtrs.push_back( &(*gen_it) );
    }
  }

  // Sort all the collections
  sort( modOfflineJetPtrs.rbegin(), modOfflineJetPtrs.rend(), EtSort<ModOfflineJet>() );
  sort( genJetPtrs.rbegin(), genJetPtrs.rend(), EtSort<GenJet>() );

  // Associate this new collection with the goodGenJets
  AssociatorEt<ModOfflineJet, GenJet > genOffAssociator( dR_ );
  auto_ptr<map<const ModOfflineJet*, const GenJet*> > genJetOfflineJetMap( genOffAssociator.Associate( modOfflineJetPtrs, genJetPtrs ) );

  // Loop on the map and write the refEt variables for the associated jets
  numChanged_ = 0;
  map<const ModOfflineJet*, const GenJet*>::const_iterator genOffMap_it = genJetOfflineJetMap->begin();
  for( ; genOffMap_it != genJetOfflineJetMap->end(); ++genOffMap_it ) {
    // Remove const from the pointer, needed to access non const method
    ModOfflineJet * nonConstOffPtr = const_cast<ModOfflineJet* >(genOffMap_it->first);
    nonConstOffPtr->setRefEt( genOffMap_it->second->et() );
    nonConstOffPtr->setRefEta( genOffMap_it->second->eta() );
    ++numChanged_;
  }
}

void Multiplier::generate() {

  // Store initial values for the MEt
  float pi = 3.1415926;
  //  float metc_dpmin = pi;
  float metc_met = 0.;
  //  float metc_smet = 0.;
  // Set initial values
  float metc_metx = origMEt_*cos(origPhi_);
  float metc_mety = origMEt_*sin(origPhi_);
  float metc_sumet = origSumEt_;

  // Loop on the modOfflineJet vector and modify the et of associated jets (those with refEt != 0)
  vector<ModOfflineJet>::iterator off_it = modOfflineJets_.begin();
  int countJet = 0;
  for( ; off_it != modOfflineJets_.end(); ++off_it, ++countJet ) {
    if ( off_it->refEt() != 0. ) {
      // Take the refEt and refEta (those of the associated GenJet)
      int ipt = (int)((off_it->refEt())/20.);
      if ( ipt>9 ) ipt=9; 
      int ieta = (int)(fabs(off_it->refEta())/1.2);
      if ( ieta>3 ) ieta=3;

      int kk=4*ipt+ieta;
      Double_t a = Ptres_[kk]->GetRandom();
      // Change the et of the jet with the new value
      off_it->setEt( off_it->refEt()*(1-a) );

      // Compute new missing Et, new number of jets passing cuts, 
      // new sumet, new missing Et significance, new delta phi min
      // with corrected IC5C jets
      // ---------------------------------------------------------
      float deltaPt = off_it->et() - off_it->origEt();
      // instead correct with an alpha factor
      metc_metx -= deltaPt*cos(off_it->phi())*alpha_;
      metc_mety -= deltaPt*sin(off_it->phi())*alpha_;
      metc_sumet += deltaPt*alpha_; // beware, we might need a correction factor (for jet calibration)
    }
  }
  // Sort the jets before returning the vector
  sort( modOfflineJets_.rbegin(), modOfflineJets_.rend(), EtSort<ModOfflineJet>() );

  // Set the new values in offlineMEt
  metc_met = sqrt(metc_metx*metc_metx+metc_mety*metc_mety);
  offlineMEt_.setEt( metc_met );
  offlineMEt_.setMEtSig( metc_met/sqrt(metc_sumet) );
  offlineMEt_.setPhi( pi+atan2(-metc_mety,-metc_metx) );
  offlineMEt_.setSumEt( metc_sumet );
}

#endif // MULTIPLIER_CC
