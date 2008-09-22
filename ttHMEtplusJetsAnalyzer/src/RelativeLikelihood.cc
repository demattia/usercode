#ifndef RELATIVELIKELIHOOD_CC
#define RELATIVELIKELIHOOD_CC

/**
 * Original Author: Marco De Mattia
 * Creation date: 25/8/2008
 * Mail: demattia@pd.infn.it
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/RelativeLikelihood.h"
#include <fstream>
#include <sstream>

using namespace edm;
using namespace std;
using namespace anaobj;

L1Trig RelativeLikelihood::L1Trigger;

// Constructors and destructor
// ---------------------------
RelativeLikelihood::RelativeLikelihood(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  cenJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "CenJets" ) ),
  forJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "ForJets" ) ),
  tauJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "TauJets" ) ),
  l1MEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "L1MEt" ) ),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  offlineMEtLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  globalMuonLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons" ) ),
  simpleElectronLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons" ) ),
  simpleTauLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus" ) ),
  summaryLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "Summary" ) ),
  vtxLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "vtxLabel" ) ),
  withL1ForwardJets_( iConfig.getUntrackedParameter<bool>("withL1ForwardJets") ),
  vtxAssoc_( iConfig.getUntrackedParameter<bool>("VtxAssoc") ),
  higgsFileName_(iConfig.getUntrackedParameter<string>("HiggsFileName") ),
  hadronicTopFileName_(iConfig.getUntrackedParameter<string>("HadronicTopFileName") ),
  qcdFileName_(iConfig.getUntrackedParameter<string>("QCDfileName") ),
  qcdHistoFileName_(iConfig.getUntrackedParameter<string>("QCDhistoFileName") ),
  jetEtCut_(iConfig.getUntrackedParameter<double>("JetEtCut") ),
  jetEtaCut_(iConfig.getUntrackedParameter<double>("JetEtaCut") ),

  countTTHdecaysFileName_(iConfig.getUntrackedParameter<string>("CountTTHdecaysFileName") ),
  countTTHdecays2tagsFileName_(iConfig.getUntrackedParameter<string>("CountTTHdecays2tagsFileName") ),

  inputFileNameSignal_(iConfig.getUntrackedParameter<string>("InputFileNameSignal") ),
  inputFileNameBackground_(iConfig.getUntrackedParameter<string>("InputFileNameBackground") ),
  outputFileName_(iConfig.getUntrackedParameter<string>("OutputFileName") ),
  tmvaSuffix_(iConfig.getUntrackedParameter<string>("TMVAsuffix") ),
  useTagMatrixForQCD_(iConfig.getUntrackedParameter<bool>("UseTagMatrixForQCD") ),
  eventCounter_(0),
  l1Eff_(0)
{
  countTTHdecays_ = new ttHdecaysCounter(countTTHdecaysFileName_);
  countTTHdecays2tags_ = new ttHdecaysCounter(countTTHdecays2tagsFileName_);

  inputFileSignal_  = new TFile(inputFileNameSignal_.c_str());
  inputFileBackground_  = new TFile(inputFileNameBackground_.c_str());

  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");

  // Directory names and suffixes
  TString suffixSignal = "2tags";
  TString suffixBackground = "2tags_tagMatrix";
  // if ( useTagMatrixForQCD_ ) suffixBackground += "_tagMatrix";
  TString dirNameSignal = "EventVariables";
  TString dirNameBackground = "EventVariables";
  if ( suffixSignal != "" ) {
    suffixSignal.Prepend("_");
    dirNameSignal.Append(suffixSignal);
  }
  if ( suffixBackground != "" ) {
    suffixBackground.Prepend("_");
    dirNameBackground.Append(suffixBackground);
  }

  // Production of pseudo-events for qcd with 2 b-tags.
  if ( useTagMatrixForQCD_ ) {
    qcdbTagMatrixMultiplier_ = new QCDbTagMatrix(higgsFileName_, hadronicTopFileName_, qcdFileName_, suffixSignal, outputFile_, false, qcdHistoFileName_, 2, tmvaSuffix_);
  }
  else {
    eventVariables2Tags_ = new EventVariables(higgsFileName_, hadronicTopFileName_, qcdFileName_, suffixSignal, outputFile_, false);
  }

  jetVertexAssociator_ = new JetVertexAssociator(jetEtCut_,jetEtaCut_);

  // Load selected histograms from the input file
  TDirectory * dirSignal = dynamic_cast<TDirectory*>(inputFileSignal_->Get(dirNameSignal));
  TDirectory * dirBackground = dynamic_cast<TDirectory*>(inputFileBackground_->Get(dirNameBackground));
  
  const int totHistogramNum = 21;
  // Variables names: ATTENTION, they must be ordered as those returned by EventVariables

  TString histogramNames[totHistogramNum] = { "higgsMass","hadronicTopMass", "hadronicWmass", "chi2ofMasses",
                                              "hadronicTopProjectionAlongHiggsDirection", "deltaEtaHadronicTopHiggs", "hadronicTopPlusHiggsMass", "remainingJetsMass",
                                              "first6jetsMass",
					      "first6jetsCentrality", "first8jetsMass", "first8jetsCentrality",
                                              "deltaPhiMEt1stLeadingJet",
					      "deltaPhiMEt2ndLeadingJet",
					      "deltaPhiMEt3rdLeadingJet",
					      "sixthJetEt",
                                              "goodHt", "mEtSig", "sumHighEffDiscriminantFirst4Jets", "sumHighEffDiscriminantFirst6Jets",
                                              "bTagTkInvMass" };

  for ( int histogramNum = 0; histogramNum < totHistogramNum; ++histogramNum ) {
    histogramVariableSignal_.push_back(*(dynamic_cast<TH1D*>(dirSignal->Get(histogramNames[histogramNum] + suffixSignal))));
    TH1D * tempHistoPtr = &(histogramVariableSignal_.back());
    tempHistoPtr->Scale(1./(tempHistoPtr->Integral()));
    histogramVariableBackground_.push_back(*(dynamic_cast<TH1D*>(dirBackground->Get(histogramNames[histogramNum] + suffixBackground))));
    tempHistoPtr = &(histogramVariableBackground_.back());
    tempHistoPtr->Scale(1./(tempHistoPtr->Integral()));
  }

  // This directory is created by the eventVariables2Tags class
  outputDir_ = dynamic_cast<TDirectory*>(outputFile_->Get(dirNameSignal));
  outputDir_->cd();
  relativeLikelihood_ = new TH1D( "relativeLikelihood" + suffixSignal, "relative likelihood " + suffixSignal, 50, -10., 10. );
}

RelativeLikelihood::~RelativeLikelihood()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void RelativeLikelihood::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // L1 Calo
  // -------
  edm::Handle < BaseJetCollection > l1eCenJets;
  edm::Handle < BaseJetCollection > l1eForJets;
  edm::Handle < BaseJetCollection > l1eTauJets;
  edm::Handle < BaseMEt > l1eEtMiss;

  try {
    iEvent.getByLabel(cenJetLabel_, l1eCenJets);
    iEvent.getByLabel(forJetLabel_, l1eForJets);
    iEvent.getByLabel(tauJetLabel_, l1eTauJets);
    iEvent.getByLabel(l1MEtLabel_, l1eEtMiss);
  }
  catch (...) {
    std::cerr << "L1TGCT: could not find one of the classes?" << std::endl;
    return;
  }

  // Global muons
  // SimpleElectrons
  // SimpleTaus
  // Summary
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  edm::Handle < SimpleTauCollection > simpleTaus;
  edm::Handle < Summary > summary;
  try {
    iEvent.getByLabel( globalMuonLabel_, globalMuons );
    iEvent.getByLabel( simpleElectronLabel_, simpleElectrons );
    iEvent.getByLabel( simpleTauLabel_, simpleTaus );
    iEvent.getByLabel( summaryLabel_, summary );
  }
  catch (...) {
    std::cerr << "One of the remaining collections cannot be found" << endl;
  }
  // --------------------------------------------
  // --- end loading handles with collections ---
  // --------------------------------------------

  eventCounter_++;
  if ( eventCounter_/100 == float(eventCounter_)/100. ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
  }
  // ------------------- //
  // L1 Trigger response //
  // ------------------- //

  // Must be fed level1: cen;for;tau;MEt in this order.
  L1Trigger.Fill(*l1eCenJets, *l1eForJets, *l1eTauJets, l1eEtMiss->et());
  // The bool asks to use or not the level 1 forward jets in the trigger.
  // The first bool is multijet, the second met+jet response.
  pair<bool, bool> responsePair(L1Trigger.Response(withL1ForwardJets_));

  // Evaluate level 1 efficiency for multijet and met+jet
  if (responsePair.first || responsePair.second) ++l1Eff_;

  // Take the offlineMEt
  edm::Handle<OfflineMEt> offlineMEt;
  iEvent.getByLabel( offlineMEtLabel_, offlineMEt );

  // ----------------------------------- //
  // Determine the type of event from MC //
  // ----------------------------------- //
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );
  // The first index is for the Higgs decay, the second for the ttbar.
  // Check ttHdecaysCounter for the decay codes.
  pair<int,int> decayType(countTTHdecays_->countDecays(*MCpartons));
  // ----------------------------------- //

  // ------------------------------------------------ //
  // Check if a good electron or muon is in the event //
  // ------------------------------------------------ //
  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );
  //  bool goodElectronFound = goodElectron( *simpleElectrons, *caloJets);
  //  bool goodMuonFound = goodMuon( *globalMuons, *caloJets);
  // ------------------------------------------------ //

  //--------------------------
  // Reco vertex collection
  //--------------------------
  Handle<BaseVertexCollection> recVtxs;
  BaseVertexCollection recoVertexes;
  try {
    iEvent.getByLabel( vtxLabel_ , recVtxs);
  } 
  catch (...) {
    std::cerr << "Could not find the recovertex collection" << std::endl;
    return;
  } 
  recoVertexes = *(recVtxs.product());

  // Select only jets with Et>(JetEtCut_)GeV and Eta<(JetEtaCut_) and write them in the goodJets collection
  vector<const OfflineJet *> goodJets;
  vector<const OfflineJet *> goodbTaggedJets;

  // Consider only the first 8 most energetic (Et) jets
  // --------------------------------------------------
  int numGoodJets = 0;

  // Jet-Vertex Algorithm
  // Need to define this outside the if clause, otherwise the objects will not survive it and the
  // pointers stored in goodJets will be invalid.
  OfflineJetCollection associatedJet;
  if(vtxAssoc_){
    associatedJet = jetVertexAssociator_->associate(*(caloJets.product()),recoVertexes);
    for ( OfflineJetCollection::const_iterator assocJetIt = associatedJet.begin(); assocJetIt != associatedJet.end() && numGoodJets < 8; ++assocJetIt) {
      //non metto i tagli in et ed eta sono  fatti internamente dall'algoritmo di associazione
      goodJets.push_back(&(*assocJetIt));
      ++numGoodJets;
      // -------------------------------------------------------------- //
      // -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
      // -------------------------------------------------------------- //
      // Consider as tagged those jets with highEff > 5.3.
      // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
      // Set the b-tag cut value
      float medium = 5.3;
	if ( assocJetIt->discriminatorHighEff()>medium ) goodbTaggedJets.push_back(&(*assocJetIt));
    }
  } else {
    for ( OfflineJetCollection::const_iterator allJetIt = caloJets->begin(); allJetIt != caloJets->end() && numGoodJets < 8; ++allJetIt ) {
      if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) {
	goodJets.push_back(&(*allJetIt));
        ++numGoodJets;
        
	// -------------------------------------------------------------- //
	// -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
	// -------------------------------------------------------------- //
	// Consider as tagged those jets with highEff > 5.3.
	// high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
	// Set the b-tag cut value
	float medium = 5.3;
	if ( allJetIt->discriminatorHighEff()>medium ) goodbTaggedJets.push_back(&(*allJetIt));
      }
    }
  }

  // Preliminary cuts
  // ----------------
  if ( goodJets.size() >= 5 && offlineMEt->corrL3mEtSig() > 3. ){ // && goodbTaggedJets.size() >= 2 ) {

    // If it is QCD, do not cut on the b-tags, but loop on all the combinations of jets
    // --------------------------------------------------------------------------------
    if ( useTagMatrixForQCD_ ) {
      cout << "calling multiply for event = " << eventCounter_ << endl;
      vector<vector<double> > eventVariablesVectors( qcdbTagMatrixMultiplier_->multiply( goodJets, &(*offlineMEt) ) );
      cout << "number of vectors = " << eventVariablesVectors.size() << endl;
      vector<vector<double> >::const_iterator eventVariablesVec = eventVariablesVectors.begin();
      for( ; eventVariablesVec != eventVariablesVectors.end(); ++eventVariablesVec ) {
        cout << "number of variable histograms for signal = " << histogramVariableSignal_.size() << endl;
        cout << "number of variable histograms for background = " << histogramVariableBackground_.size() << endl;
        cout << "number of variables in vector = " << eventVariablesVec->size() << endl;
        if (histogramVariableSignal_.size() != eventVariablesVec->size() || histogramVariableBackground_.size() != eventVariablesVec->size()) cout << "ATTENTION: number of variables does not match number of loaded histograms." << endl;
        evaluateLikelihood( *eventVariablesVec );
      }
    }
    // For normal events evaluate likelihood
    // -------------------------------------
    else {
      // Require at least two b-tags
      if ( goodbTaggedJets.size() >= 2 ) {
        pair<int,int> decayType2tags(countTTHdecays2tags_->countDecays(*MCpartons));
        vector<double> eventVariablesVector( eventVariables2Tags_->fill( goodJets, goodbTaggedJets, &(*offlineMEt) ) );
        if (histogramVariableSignal_.size() != eventVariablesVector.size() || histogramVariableBackground_.size() != eventVariablesVector.size()) cout << "ATTENTION: number of variables does not match number of loaded histograms." << endl;
        evaluateLikelihood( eventVariablesVector );
      } // end if at least two b-tags
    }
  } // end preliminary cuts

}

void RelativeLikelihood::evaluateLikelihood( const vector<double> & eventVariablesVector ) {

  // ATTENTION: the exclusion is done by hand here, to skip the variables we do not want to use in the likelihood.

  double ratiosProduct = 1.;
  int varCounter = 0;
  for ( vector<double>::const_iterator varIter = eventVariablesVector.begin(); varIter != eventVariablesVector.end(); ++varIter, ++varCounter ) {

    // Skip this variable if it is in the blacklist
    TString histoName = histogramVariableSignal_[varCounter].GetName();
    if( histoName == "first6jetsMass" ||
        histoName == "deltaPhiMEt1stLeadingJet" ||
        histoName == "deltaPhiMEt3rdLeadingJet" ) continue;

    cout << "histoName["<<varCounter<<"] = " << histogramVariableSignal_[varCounter].GetName() << endl;
    // cout << "varIter["<<varCounter<<"] = " << *varIter << endl;
    // cout << "BinWidth["<<varCounter<<"] = " << histogramVariableSignal_[varCounter]->GetBinWidth(1) << endl;
    // cout << "BinLowEdge["<<varCounter<<"] = " << histogramVariableSignal_[varCounter]->GetBinLowEdge(1) << endl;

    double binWidthSignal = histogramVariableSignal_[varCounter].GetBinWidth(1);
    double binWidthBackground = histogramVariableBackground_[varCounter].GetBinWidth(1);

    int binNumberSignal = histogramVariableSignal_[varCounter].GetNbinsX();
    int binNumberBackground = histogramVariableBackground_[varCounter].GetNbinsX();

    // double upperEdge = histogramVariableSignal_[varCounter]->GetBinLowEdge(binNumber) + binWidthSignal << endl;
    // cout << "UpperEdge["<<varCounter<<"] = " << lastBin << endl;

    // Determine the bin index
    // First translate the variable of the negative edge of the histogram: (-3, 3) -> translate the variable by 3.
    double shiftedVar = *varIter + histogramVariableSignal_[varCounter].GetBinLowEdge(1);
    // Now divide by the bin width
    int binIndexSignal = int(shiftedVar/binWidthSignal)+1;
    int binIndexBackground = int(shiftedVar/binWidthBackground)+1;
    cout << "binIndexSignal = " << binIndexSignal << endl;
    cout << "binIndexBackground = " << binIndexBackground << endl;
    double signalVar = 1.;
    double backgroundVar = 1.;
    if ( binIndexSignal >= 0 && binIndexSignal <= binNumberSignal ) signalVar = histogramVariableSignal_[varCounter].GetBinContent( binIndexSignal );
    if ( binIndexBackground > 0 && binIndexBackground < binNumberBackground ) backgroundVar = histogramVariableBackground_[varCounter].GetBinContent( binIndexBackground );

    if ( backgroundVar == 0 ) {
      if ( signalVar == 0 ) ratiosProduct *= 1.;
      // ATTENTION: this is arbitrarily put to 10, consider a more appropriate value
      else ratiosProduct *= 100.;
    }
    else if ( signalVar == 0 ) ratiosProduct *= 0.01;
    else ratiosProduct *= signalVar/backgroundVar;

    cout << "signalVar = " << signalVar << endl;
    cout << "backgroundVar = " << backgroundVar << endl;

    cout << "ratiosProduct = " << ratiosProduct << endl;
  }

  double relativeLikelihood = 0.;

  if ( ratiosProduct >= 0. ) relativeLikelihood = log(ratiosProduct);
  else relativeLikelihood = -9.99;
  if ( relativeLikelihood < -10. ) relativeLikelihood = -9.99;
  if ( relativeLikelihood > 10. ) relativeLikelihood = 9.99;

  cout << "relativeLikelihood = ";

  relativeLikelihood_->Fill(relativeLikelihood);

  cout << relativeLikelihood << endl;
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void RelativeLikelihood::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void RelativeLikelihood::endJob() {

  countTTHdecays_->writeDecays();
  countTTHdecays2tags_->writeDecays();
  delete countTTHdecays_;
  delete countTTHdecays2tags_;
  if ( useTagMatrixForQCD_ ) { delete qcdbTagMatrixMultiplier_; }
  else { delete eventVariables2Tags_; }
  outputDir_->cd();
  relativeLikelihood_->Write();
  outputFile_->Write();
}

// See if we have a good high-Pt electron or muon
// ----------------------------------------------
bool RelativeLikelihood::goodMuon( const GlobalMuonCollection & globalMuons, const OfflineJetCollection & caloJets ) {
  bool muonEvent = false;
  for ( GlobalMuonCollection::const_iterator muon = globalMuons.begin(); 
	muon != globalMuons.end() && !muonEvent; ++muon ) {
    if ( muon->pt()>25 && fabs(muon->eta())<2.5 ) { 
      // See if there are jets closer than 0.5 from this one
      // ---------------------------------------------------
      double dRmin=25.;
      for ( OfflineJetCollection::const_iterator cal = caloJets.begin(); 
	    cal != caloJets.end(); ++cal ) {
	if ( cal->et()>25. ) { 
	  double dR = DeltaR(cal->eta(), cal->phi(), muon->eta(), muon->phi());
	  if ( dR<dRmin ) dRmin=dR;
	}
      }
      if ( dRmin>0.25 ) muonEvent = true;
    }
  }
  return muonEvent;
}
bool RelativeLikelihood::goodElectron( const SimpleElectronCollection & simpleElectrons, const OfflineJetCollection & caloJets ) {
  bool elecEvent = false;
  for ( SimpleElectronCollection::const_iterator elec = simpleElectrons.begin(); 
	elec != simpleElectrons.end() && !elecEvent; ++elec ) {
    if ( elec->et()>30 && fabs(elec->eta())<2.5 && elec->hadOverEm()<0.05 ) {
      // See if there are jets closer than 0.5 from this one
      // ---------------------------------------------------
      double dRmin=25.;
      for ( OfflineJetCollection::const_iterator cal = caloJets.begin(); 
	    cal != caloJets.end(); ++cal ) {
	if ( cal->et()>25. && cal->emEnergyFraction()<0.95 ) { 
	  double dR = DeltaR(cal->eta(), cal->phi(), elec->eta(), elec->phi());
	  if ( dR<dRmin ) dRmin=dR;
	}
      }
      if ( dRmin>0.25 ) elecEvent = true;
    }
  }
  return elecEvent;
}

#endif // RELATIVELIKELIHOOD_CC

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RelativeLikelihood);
