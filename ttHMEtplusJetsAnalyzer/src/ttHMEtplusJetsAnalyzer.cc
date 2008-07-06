#ifndef TTHMETPLUSJETSANALYZER_CC
#define TTHMETPLUSJETSANALYZER_CC

/**
 * Original Author: Marco De Mattia
 * Creation date: 29/6/2008
 * Mail: demattia@pd.infn.it
 */

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/ttHMEtplusJetsAnalyzer.h"
#include <fstream>
#include <sstream>

using namespace edm;
using namespace std;
using namespace anaobj;

// Small structs used in the analyze method. Defined in a nameless namespace so that they only exist in this file.
namespace {
  struct pairStruct {
    pairStruct( const OfflineJet * JET1, const OfflineJet * JET2 ) {
      //, const double & PROB ) {
      jet1 = JET1;
      jet2 = JET2;
    }
    double prob;
    const OfflineJet * jet1;
    const OfflineJet * jet2;
    /// To sort the structs
    bool operator< ( const pairStruct& b ) const {
      return prob < b.prob;
    }
  };
  struct tripletStruct : public pairStruct {
    tripletStruct( const OfflineJet * JET1, const OfflineJet * JET2, const OfflineJet * JET3 ) : pairStruct(JET1, JET2) {
      jet3 = JET3;
    }
    const OfflineJet * jet3;
    bool operator< ( const tripletStruct& b ) const {
      return prob < b.prob;
    }
  };
}

L1Trig ttHMEtplusJetsAnalyzer::L1Trigger;

// Constructors and destructor
// ---------------------------
ttHMEtplusJetsAnalyzer::ttHMEtplusJetsAnalyzer(const edm::ParameterSet& iConfig) :
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
  withL1ForwardJets_( iConfig.getUntrackedParameter<bool>("withL1ForwardJets") ),
  higgsFileName_(iConfig.getUntrackedParameter<string>("HiggsFileName") ),
  hadronicTopFileName_(iConfig.getUntrackedParameter<string>("HadronicTopFileName") ),
  qcdFileName_(iConfig.getUntrackedParameter<string>("QCDfileName") ),
  jetEtCut_(iConfig.getUntrackedParameter<double>("JetEtCut") ),
  jetEtaCut_(iConfig.getUntrackedParameter<double>("JetEtaCut") ),
  eventCounter_(0),
  l1Eff_(0)
{
  countTTHdecays_ = new ttHdecaysCounter("ttHdecays.txt");

  // Load the matrices for Higgs, hadronic top and QCD b-tag probabiliy from the files
  // ---------------------------------------------------------------------------------
  fillProbabilityMatrices( higgsFileName_, higgsBinNum_, higgsBinSize_, trueH_, falseH_ );
  fillProbabilityMatrices( hadronicTopFileName_, hadronicTopBinNum_, hadronicTopBinSize_, trueHadronicTop_, falseHadronicTop_ );
  fillProbabilityMatrices( qcdFileName_, qcdBinNum_, qcdBinSize_, taggedJet_, notTaggedJet_ );

}

ttHMEtplusJetsAnalyzer::~ttHMEtplusJetsAnalyzer()
{
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void ttHMEtplusJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  // Select only jets with Et>(JetEtCut_)GeV and Eta<(JetEtaCut_) and write them in the goodJets collection
  vector<const OfflineJet *> goodJets;
  vector<const OfflineJet *> goodbTaggedJets;
  for ( OfflineJetCollection::const_iterator allJetIt = caloJets->begin(); allJetIt != caloJets->end(); ++allJetIt ) {
    if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) {
      goodJets.push_back(&(*allJetIt));

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

  // Require at least two b-tags
  if ( goodbTaggedJets.size() >= 2 ) {

    // Create pairs of b-jets and evaluate their probability to come from the Higgs decay
    //  vector<pair<true/false ratio, candidate> >
    vector<pair<double, Particle<const OfflineJet> > > bTaggedPairs;
    vector<const OfflineJet *>::const_iterator bTaggedJetIt = goodbTaggedJets.begin();
    for ( ; bTaggedJetIt != goodbTaggedJets.end(); ++bTaggedJetIt ) {
      vector<const OfflineJet *>::const_iterator subbTaggedJetIt = bTaggedJetIt+1;
      for ( ; subbTaggedJetIt != goodbTaggedJets.end(); ++subbTaggedJetIt ) {
        // bTaggedPairs.push_back(pairStruct( *bTaggedJetIt, *subbTaggedJetIt ));
        Particle<const OfflineJet> higgsCandidate( *bTaggedJetIt );
        higgsCandidate.add( *subbTaggedJetIt );
        bTaggedPairs.push_back( make_pair( evalHiggsPairProbability(higgsCandidate), higgsCandidate ) ); 
        cout << "true/false ratio = " << bTaggedPairs.back().first << endl;
      }
    }






  } // end if at least two b-tags




  // Associate partons to offlineJets only for H->bb/cc and tt->leptons+4jets decays
//  int higgsDecayType = decayType.first;
//  int ttDecayType = decayType.second;
//  if ( (higgsDecayType == 0 || higgsDecayType == 1) && (ttDecayType == 11 || ttDecayType == 101 || ttDecayType == 1001) ) {
//    cout << "ttH->MEt+4Jets decay" << endl;
//  }

  // Select the Higgs using the true/falseHiggs matrices

}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void ttHMEtplusJetsAnalyzer::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void ttHMEtplusJetsAnalyzer::endJob() {

  countTTHdecays_->writeDecays();
  delete countTTHdecays_;
}

// See if we have a good high-Pt electron or muon
// ----------------------------------------------
bool ttHMEtplusJetsAnalyzer::goodMuon( const GlobalMuonCollection & globalMuons, const OfflineJetCollection & caloJets ) {
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
bool ttHMEtplusJetsAnalyzer::goodElectron( const SimpleElectronCollection & simpleElectrons, const OfflineJetCollection & caloJets ) {
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



void ttHMEtplusJetsAnalyzer::fillProbabilityMatrices(const string & probabilityFileName, unsigned int * binNum, double * binSize, unsigned int ***& trueArray, unsigned int ***&falseArray ) {

  ifstream probabilityFile(probabilityFileName.c_str());
  if ( !probabilityFile.is_open() ) {
    cout << "file: " << probabilityFileName << "not found or unable to open" << endl;
    exit(1);
  }

  string line;
  // Read the first line with the bin numbers and sizes
  getline(probabilityFile, line);
  stringstream probabilityCounts(line);
  // This loop is based on the text file structure
  string skip;
  int wordCount=0;
  int binNumCount=0;
  int binSizeCount=0;
  while ( !probabilityCounts.eof() ) {
    ++wordCount;
    // Take every three words
    if ( wordCount%3 == 0 ) {
      // Every six words there is the size of the bins, in the other case is the number of bins
      if ( wordCount%6 != 0 ) {
        probabilityCounts >> binNum[binNumCount];
        // cout << "binNum_["<<binNumCount<<"] = " << binNum[binNumCount] << endl;
        ++binNumCount;
      }
      else {
        probabilityCounts >> binSize[binSizeCount];
        // cout << "binSize_["<<binSizeCount<<"] = " << binSize[binSizeCount] << endl;
        ++binSizeCount;
      }
    }
    // Skip the rest of the words
    else {
      probabilityCounts >> skip;
    }
  }

  // Create and fill the arrays
  trueArray = new unsigned int**[binNum[0]];
  falseArray = new unsigned int**[binNum[0]];
  for(unsigned int i=0; i < binNum[0]; ++i) {
    trueArray[i] = new unsigned int*[binNum[1]];
    falseArray[i] = new unsigned int*[binNum[1]];
    for(unsigned int j=0; j < binNum[1]; ++j) {
      trueArray[i][j] = new unsigned int[binNum[2]];
      falseArray[i][j] = new unsigned int[binNum[2]];
      for (unsigned int k=0; k < binNum[2]; ++k) {
        if ( probabilityFile.eof() ) { 
          cout << "ERROR: not enough lines in the file for the required bins" << endl;
          exit(1);
        }

        getline(probabilityFile, line);
        stringstream probabilityCounts(line);
        for ( int w=0; w<6; ++w ) {
          // The second word is for the true case probability
          if( w == 2 ) {
            probabilityCounts >> trueArray[i][j][k];
            // cout << "trueArray["<<i<<"]["<<j<<"]["<<k<<"] = " << trueArray[i][j][k];
            // The fourth word is for the false case probability
          }
          if( w == 4 ) {
            probabilityCounts >> falseArray[i][j][k];
            // cout << " trueArray["<<i<<"]["<<j<<"]["<<k<<"] = " << falseArray[i][j][k] << endl;;
          }
          // Skip the other words
          else probabilityCounts >> skip;
        }
      }

    }
  }
}

double ttHMEtplusJetsAnalyzer::evalHiggsPairProbability(const Particle<const OfflineJet> & higgsCandidate) const {

  // Determine bin index
  int etId = int(higgsCandidate.pt()/higgsBinSize_[0]);
  int etaId = int(fabs(higgsCandidate.eta())/higgsBinSize_[1]);
  if ( higgsCandidate.components().size() != 2 ) {
    cout << "Error: higgsCandidate does not have 2 component jets" << endl;
    exit(1);
  }
  const OfflineJet * jet1 = (higgsCandidate.components())[0];
  const OfflineJet * jet2 = (higgsCandidate.components())[1];
  int dRId = int( DeltaR( jet1->eta(), jet1->phi(), jet2->eta(), jet2->phi())/higgsBinSize_[2] );
  if ( etId < 0 || etaId < 0 || dRId < 0 ) cout << "Error: index < 0, will crash..." << endl;
  if ( etId>int(higgsBinNum_[0]-1) ) etId = higgsBinNum_[0]-1;
  if ( etaId>int(higgsBinNum_[1]-1) ) etaId = higgsBinNum_[1]-1;
  if ( dRId>int(higgsBinNum_[2]-1) ) dRId = higgsBinNum_[2]-1;

  cout << "before crash" << endl;
  cout << "etId = " << etId << endl;
  cout << "etaId = " << etaId << endl;
  cout << "dRId = " << dRId << endl;
  cout << "trueH_["<<etId<<"]["<<etaId<<"]["<<dRId<<"] = " << trueH_[etId][etaId][dRId] << endl;
  cout << "falseH_["<<etId<<"]["<<etaId<<"]["<<dRId<<"] = " << falseH_[etId][etaId][dRId] << endl;

  cout << "before before crash" << endl;

  return (double(trueH_[etId][etaId][dRId])/double(falseH_[etId][etaId][dRId]));
}


#endif // TTHMETPLUSJETSANALYZER_CC

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ttHMEtplusJetsAnalyzer);
