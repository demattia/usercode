#ifndef QCDBTAGPROBABILITY_CC
#define QCDBTAGPROBABILITY_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/QCDbTagProbability.h"
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"

#include <fstream>
#include <memory>

QCDbTagProbability::QCDbTagProbability(const ParameterSet & iConfig) :
  conf_(iConfig),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  jetEtCut_(iConfig.getUntrackedParameter<double>( "JetEtCut" ) ),
  jetEtaCut_(iConfig.getUntrackedParameter<double>( "JetEtaCut" ) ),
  maxConsideredJets_(iConfig.getUntrackedParameter<unsigned int>( "MaxConsideredJets" ) ),
  etBinNum_( iConfig.getUntrackedParameter<unsigned int>( "EtBinNum" ) ),
  etBinSize_( iConfig.getUntrackedParameter<double>( "EtBinSize" ) ),
  etaBinNum_( iConfig.getUntrackedParameter<unsigned int>( "EtaBinNum" ) ),
  etaBinSize_( iConfig.getUntrackedParameter<double>( "EtaBinSize" ) ),
  s1BinNum_( iConfig.getUntrackedParameter<unsigned int>( "S1BinNum" ) ),
  s1BinSize_( iConfig.getUntrackedParameter<double>( "S1BinSize" ) )
{
  eventCounter_ = 0;

  // Name of the file for the probability counts
  // -------------------------------------------
  outputProbabilityFileName_ = conf_.getUntrackedParameter<std::string>("OutputProbabilityFileName");

  // File for output histograms
  // --------------------------
  outputFile_ = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str(), "RECREATE");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  outputFile_->cd(); 

  // Create histograms
  jetVsMCpEt_ = new TProfile("jetVsMCpEt", "jet-et vs parton-pt", 100, 0, 600, 0, 600);

  // True and false Higgs pair histograms
  taggedJetEt_ = new TH1F("taggedJetEt", "Et of b-tagged jets", 100, 0, 200);
  taggedJetEta_ = new TH1F("taggedJetEta", "Eta of b-tagged jets", 100, -3, 3);
  taggedJetS1_ = new TH1F("taggedJetS1", "Number of tracks with significance above 1 for b-tagged jets", 20, 0, 20);
  taggedJetTagMass_ = new TH1F("taggedJetTagMass", "Mass evaluated using the associated tracks for b-tagged jets", 100, 0, 100);
  notTaggedJetEt_ = new TH1F("notTaggedJetEt", "Et of non b-tagged jets", 100, 0, 200);
  notTaggedJetEta_ = new TH1F("notTaggedJetEta", "Eta of non b-tagged jets", 100, -3, 3);
  notTaggedJetS1_ = new TH1F("notTaggedJetS1", "Number of tracks with significance above 1 for non b-tagged jets", 20, 0, 20);
  notTaggedJetTagMass_ = new TH1F("notTaggedJetTagMass", "Mass evaluated using the associated tracks for non b-tagged jets", 100, 0, 100);

  // Initialized to 1, not 0, since they are used to evaluate fractions
  // ------------------------------------------------------------------
  taggedJet_ = new unsigned int**[etBinNum_];
  notTaggedJet_ = new unsigned int**[etBinNum_];
  for(unsigned int i=0; i != etBinNum_; ++i) {
    taggedJet_[i] = new unsigned int*[etaBinNum_];
    notTaggedJet_[i] = new unsigned int*[etaBinNum_];
    for(unsigned int j=0; j != etaBinNum_; ++j) {
      taggedJet_[i][j] = new unsigned int[s1BinNum_];
      notTaggedJet_[i][j] = new unsigned int[s1BinNum_];
      for(unsigned int k=0; k != s1BinNum_; ++k) {
        taggedJet_[i][j][k] = 1;
        notTaggedJet_[i][j][k] = 1;
      }
    }
  }
}

QCDbTagProbability::~QCDbTagProbability() {
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void QCDbTagProbability::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  eventCounter_++;
  if ( eventCounter_/100 == float(eventCounter_)/100. ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
  }

  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );

  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Select only jets with Et>(JetEtCut_)GeV and Eta<(JetEtaCut_) and write them in the goodJets collection
  vector<const OfflineJet *> goodJets;
  for ( OfflineJetCollection::const_iterator allJetIt = caloJets->begin(); allJetIt != caloJets->end(); ++allJetIt ) {
    if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) goodJets.push_back(&(*allJetIt));
  }

  typedef map<const MCParticle*, const OfflineJet*> mapMCpJet;

//  // jetEt vs partonEt histogram
//  AssociatorEt<MCParticle, OfflineJet> assocJetMCpartons(0.5);
//  auto_ptr<mapMCpJet> mapPtr(assocJetMCpartons.Associate(MCpartons, goodJets));
//  for ( mapMCpJet::const_iterator mapIt = mapPtr->begin(); mapIt != mapPtr->end(); ++mapIt ) {
//    // cout << "MC particle = " << mapIt->first->pid() << " of pt = " << mapIt->first->pt() << " associated to jet of et = " << mapIt->second->et() << endl;
//    // Resolution histogram for jets associated to quarks
//    jetVsMCpEt_->Fill(mapIt->first->pt(), mapIt->second->et());
//  }

  vector<const OfflineJet *>::const_iterator goodJetIt = goodJets.begin();
  unsigned int goodJetNum = 0;
  for ( ; goodJetIt != goodJets.end() && goodJetNum < maxConsideredJets_; ++goodJetIt, ++goodJetNum ) {
    int numTkS1 = 0;
    edm::RefVector<std::vector<SimpleTrack> > tkRefVec((*goodJetIt)->tkRefVec());
    edm::RefVector<std::vector<SimpleTrack> >::const_iterator tkRefIt = tkRefVec.begin();
    for ( ; tkRefIt != tkRefVec.end(); ++tkRefIt ) {
      if ( (*tkRefIt)->ip2Dsignificance() >= 1. ) ++numTkS1;
    }

    double jetEt = (*goodJetIt)->et();
    double jetEta = (*goodJetIt)->eta();
    // Determine bin index
    int etId = int(jetEt/etBinSize_);
    int etaId = int(fabs(jetEta)/etaBinSize_);
    int numTkS1Id = int(numTkS1/s1BinSize_);
    if ( etId < 0 || etaId < 0 || numTkS1Id < 0 ) cout << "Error: index < 0, will crash..." << endl;
    if ( etId>int(etBinNum_-1) ) etId = etBinNum_-1;
    if ( etaId>int(etaBinNum_-1) ) etaId = etaBinNum_-1;
    if ( numTkS1Id>int(s1BinNum_-1) ) numTkS1Id = s1BinNum_-1;

    // -------------------------------------------------------------- //
    // -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
    // -------------------------------------------------------------- //
    // Consider as tagged those jets with highEff > 5.3.
    // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
    // Set the b-tag cut value
    float medium = 5.3;
    if ( (*goodJetIt)->discriminatorHighEff()>medium ) {

      taggedJetEt_->Fill(jetEt);
      taggedJetEta_->Fill(jetEta);
      taggedJetS1_->Fill(numTkS1);
      taggedJetTagMass_->Fill( (*goodJetIt)->bTagTkInvMass() );
      (taggedJet_[etId][etaId][numTkS1Id])++;
    }
    // else fill HcombFalse.txt
    else {
      notTaggedJetEt_->Fill(jetEt);
      notTaggedJetEta_->Fill(jetEta);
      notTaggedJetS1_->Fill(numTkS1);
      notTaggedJetTagMass_->Fill( (*goodJetIt)->bTagTkInvMass() );
      (notTaggedJet_[etId][etaId][numTkS1Id])++;
    }
  }

  // Write the trueH_ and falseH_ to a txt file
  // ------------------------------------------
  ofstream bTagProbabilityFile(outputProbabilityFileName_.c_str());
  // The first line has informations on bin numbers and sizes
  bTagProbabilityFile <<  "etBinNum = "  << etBinNum_  << " etBinSize = "  << etBinSize_
                      << " etaBinNum = " << etaBinNum_ << " etaBinSize = " << etaBinSize_
                      << " s1BinNum = "  << s1BinNum_  << " s1BinSize = "  << s1BinSize_ 
                      << " totEvents = " << eventCounter_ << endl;

  // double taggedCount = 0.;
  // double notTaggedCount = 0.;
  // double norm = 0.;
  // The following lines have the counts for b-tagged and non b-tagged jets
  for(unsigned int i=0; i != etBinNum_; ++i) {
    for(unsigned int j=0; j != etaBinNum_; ++j) {
      for(unsigned int k=0; k != s1BinNum_; ++k) {
        // Evaluate the probability: tagged/(tagged+notTagged) and notTagged/(tagged+notTagged)
        // taggedCount = taggedJet_[i][j][k];
        // notTaggedCount = notTaggedJet_[i][j][k];
        // norm = taggedCount + notTaggedCount;

        // ATTENTION: the first line must not have space, while the second one must have it.
        bTagProbabilityFile << "taggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << taggedJet_[i][j][k]
                            << " notTaggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << notTaggedJet_[i][j][k] << endl;
        //                  <<  "taggedJet["<<i<<"]["<<j<<"]["<<k<<"] = "  << taggedCount/norm
        //                  << " notTaggedJet["<<i<<"]["<<j<<"]["<<k<<"] = " << notTaggedCount/norm
      }
    }
  }
  bTagProbabilityFile.close();
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void QCDbTagProbability::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void QCDbTagProbability::endJob() {

  // delete the multidimensional arrays
  for(unsigned int i=0; i != etBinNum_; ++i) {
    for(unsigned int j=0; j != etaBinNum_; ++j) {
      delete[] taggedJet_[i][j];
      delete[] notTaggedJet_[i][j];
    }
    delete[] taggedJet_[i];
    delete[] notTaggedJet_[i];
  }

  // Save histograms
  jetVsMCpEt_->Write();

  taggedJetEt_->Write();
  taggedJetEta_->Write();
  taggedJetS1_->Write();
  notTaggedJetEt_->Write();
  notTaggedJetEta_->Write();
  notTaggedJetS1_->Write();

  outputFile_->Write();
}

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(QCDbTagProbability);

#endif // QCDBTAGPROBABILITY_CC
