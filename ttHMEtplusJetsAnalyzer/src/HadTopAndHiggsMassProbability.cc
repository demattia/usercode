#ifndef HADTOPANDHIGGSMASSPROBABILITY_CC
#define HADTOPANDHIGGSMASSPROBABILITY_CC

#include "AnalysisExamples/ttHMEtplusJetsAnalyzer/interface/HadTopAndHiggsMassProbability.h"
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"
#include "AnalysisExamples/AnalysisClasses/interface/SumJets.h"

#include <fstream>
#include <memory>

HadTopAndHiggsMassProbability::HadTopAndHiggsMassProbability(const ParameterSet & iConfig) :
  conf_(iConfig),
  offlineJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets" ) ),
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "MCParticles" ) ),
  jetEtCut_(iConfig.getUntrackedParameter<double>( "JetEtCut" ) ),
  jetEtaCut_(iConfig.getUntrackedParameter<double>( "JetEtaCut" ) ),
  higgsEtBinNum_( iConfig.getUntrackedParameter<unsigned int>( "HiggsEtBinNum" ) ),
  higgsEtBinSize_( iConfig.getUntrackedParameter<double>( "HiggsEtBinSize" ) ),
  higgsEtaBinNum_( iConfig.getUntrackedParameter<unsigned int>( "HiggsEtaBinNum" ) ),
  higgsEtaBinSize_( iConfig.getUntrackedParameter<double>( "HiggsEtaBinSize" ) ),
  higgsDRbinNum_( iConfig.getUntrackedParameter<unsigned int>( "HiggsDRbinNum" ) ),
  higgsDRbinSize_( iConfig.getUntrackedParameter<double>( "HiggsDRbinSize" ) ),
  hadronicTopEtBinNum_( iConfig.getUntrackedParameter<unsigned int>( "HadronicTopEtBinNum" ) ),
  hadronicTopEtBinSize_( iConfig.getUntrackedParameter<double>( "HadronicTopEtBinSize" ) ),
  hadronicTopEtaBinNum_( iConfig.getUntrackedParameter<unsigned int>( "HadronicTopEtaBinNum" ) ),
  hadronicTopEtaBinSize_( iConfig.getUntrackedParameter<double>( "HadronicTopEtaBinSize" ) ),
  dPhiHiggsHadronicTopBinNum_( iConfig.getUntrackedParameter<unsigned int>( "DPhiHiggsHadronicTopBinNum" ) ),
  dPhiHiggsHadronicTopBinSize_( iConfig.getUntrackedParameter<double>( "DPhiHiggsHadronicTopBinSize" ) )
{
  countTTHdecays_ = new ttHdecaysCounter("ttHdecays.txt");

  eventCounter_ = 0;

  // File for output histograms
  // --------------------------
  outputFile_ = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str(), "RECREATE");
  // The file must be opened first, so that becomes the default position for all the histograms
  // ------------------------------------------------------------------------------------------
  outputFile_->cd(); 

  // Create histograms
  jetVsMCpEt_ = new TProfile("jetVsMCpEt", "jet-et vs parton-pt", 100, 0, 600, 0, 600);
  higgsMassTrue_ = new TH1F("higgsMassTrue", "mass of the two b-jets associated to the MC Higgs", 100, 0, 200);
  hadronicTopMassTrue_ = new TH1F("hadronicTopMassTrue", "mass of the three jets associated to the MC top", 100, 0, 300);

  // True and false Higgs pair histograms
  trueHiggsPairEt_ = new TH1F("trueHiggsPairEt", "Et of the two b-jets associated to the MC Higgs", 100, 0, 200);
  trueHiggsPairEta_ = new TH1F("trueHiggsPairEta", "Eta of the two b-jets associated to the MC Higgs", 100, -3, 3);
  trueHiggsPairDR_ = new TH1F("trueHiggsPairDR", "#DeltaR of the two b-jets associated to the MC Higgs", 100, 0, 10);
  falseHiggsPairEt_ = new TH1F("falseHiggsPairEt", "Et of two b-jets not associated to the MC Higgs", 100, 0, 200);
  falseHiggsPairEta_ = new TH1F("falseHiggsPairEta", "Eta of two b-jets not associated to the MC Higgs", 100, -3, 3);
  falseHiggsPairDR_ = new TH1F("falseHiggsPairDR", "#DeltaR of two b-jets not associated to the MC Higgs", 100, 0, 10);
  // True and false hadronic top triplet histograms
  trueHadronicTopTripletEt_ = new TH1F("trueHadronicTopTripletEt", "Et of the three jets associated to the MC hadronic top", 100, 0, 300);
  trueHadronicTopTripletEta_ = new TH1F("trueHadronicTopTripletEta", "Eta of the three jets associated to the MC hadronic top", 100, -3, 3);
  trueHadronicTopTripletDphiHiggsHadronicTop_ = new TH1F("trueHadronicTopTripletDphiHiggsHadronicTop",
                                                         "#Delta#Phi of the three jets associated to the MC hadronic top", 100, 0, 3.15);
  falseHadronicTopTripletEt_ = new TH1F("falseHadronicTopTripletEt", "Et of the three jets not associated to the MC hadronic top", 100, 0, 300);
  falseHadronicTopTripletEta_ = new TH1F("falseHadronicTopTripletEta", "Eta of the three jets not associated to the MC hadronic top", 100, -3, 3);
  falseHadronicTopTripletDphiHiggsHadronicTop_ = new TH1F("falseHadronicTopTripletDphiHiggsHadronicTop",
                                                          "#Delta#Phi of the three jets not associated to the MC hadronic top", 100, 0, 3.15);


  // Initialized to 1, not 0, since they are used to evaluate fractions
  // ------------------------------------------------------------------
  trueH_ = new unsigned int**[higgsEtBinNum_];
  falseH_ = new unsigned int**[higgsEtBinNum_];
  for(unsigned int i=0; i != higgsEtBinNum_; ++i) {
    trueH_[i] = new unsigned int*[higgsEtaBinNum_];
    falseH_[i] = new unsigned int*[higgsEtaBinNum_];
    for(unsigned int j=0; j != higgsEtaBinNum_; ++j) {
      trueH_[i][j] = new unsigned int[higgsDRbinNum_];
      falseH_[i][j] = new unsigned int[higgsDRbinNum_];
      for(unsigned int k=0; k != higgsDRbinNum_; ++k) {
        trueH_[i][j][k] = 1;
        falseH_[i][j][k] = 1;
      }
    }
  }

  trueHadronicTop_ = new unsigned int**[hadronicTopEtBinNum_];
  falseHadronicTop_ = new unsigned int**[hadronicTopEtBinNum_];
  for(unsigned int i=0; i != hadronicTopEtBinNum_; ++i) {
    trueHadronicTop_[i] = new unsigned int*[hadronicTopEtaBinNum_];
    falseHadronicTop_[i] = new unsigned int*[hadronicTopEtaBinNum_];
    for(unsigned int j=0; j != hadronicTopEtaBinNum_; ++j) {
      trueHadronicTop_[i][j] = new unsigned int[dPhiHiggsHadronicTopBinNum_];
      falseHadronicTop_[i][j] = new unsigned int[dPhiHiggsHadronicTopBinNum_];
      for(unsigned int k=0; k != dPhiHiggsHadronicTopBinNum_; ++k) {
        trueHadronicTop_[i][j][k] = 1;
        falseHadronicTop_[i][j][k] = 1;
      }
    }
  }
}

HadTopAndHiggsMassProbability::~HadTopAndHiggsMassProbability() {
  // Do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  // -----------------------------------------------------------
}

// Member functions
// ----------------

// ------------ method called to for each event  ------------
// ----------------------------------------------------------
void HadTopAndHiggsMassProbability::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace anaobj;

  eventCounter_++;
  if ( eventCounter_/100 == float(eventCounter_)/100. ) {
    std::cout << "Event number " << eventCounter_ << std::endl;
  }

  // ----------------------------------- //
  // Determine the type of event from MC //
  // ----------------------------------- //
  edm::Handle < MCParticleCollection > MCpartons;
  iEvent.getByLabel( MCParticleLabel_, MCpartons );
  // The first index is for the Higgs decay, the second for the ttbar.
  // Check ttHdecaysCounter for the decay codes.
  pair<int,int> decayType(countTTHdecays_->countDecays(*MCpartons));
  // Take the partons from the decay of the ttH
  vector<const MCParticle *> ttHpartons(countTTHdecays_->hadronicPartons());
  // ----------------------------------- //

  edm::Handle<OfflineJetCollection> caloJets;
  iEvent.getByLabel( offlineJetLabel_, caloJets );

  // Select only jets with Et>(JetEtCut_)GeV and Eta<(JetEtaCut_) and write them in the goodJets collection
  vector<const OfflineJet *> goodJets;
  for ( OfflineJetCollection::const_iterator allJetIt = caloJets->begin(); allJetIt != caloJets->end(); ++allJetIt ) {
    if ( allJetIt->et() >= jetEtCut_ && fabs(allJetIt->eta())< jetEtaCut_ ) goodJets.push_back(&(*allJetIt));
  }

  typedef map<const MCParticle*, const OfflineJet*> mapMCpJet;

  // Create the functor used to evaluate the sum of the particles
  SumJets<BaseJet> sumJets;

  // Associate partons to offlineJets only for H->bb/cc and tt->leptons+4jets decays
  int higgsDecayType = decayType.first;
  int ttDecayType = decayType.second;
  if ( (higgsDecayType == 0 || higgsDecayType == 1) && (ttDecayType == 11 || ttDecayType == 101 || ttDecayType == 1001) ) {

    // To hold the true Higgs b-jet pair (if any). Needed by the top triplet.
    auto_ptr<BaseJet> trueHiggsPair;

    AssociatorEt<MCParticle, OfflineJet> assocJetMCpartons(0.5);
    auto_ptr<mapMCpJet> mapPtr(assocJetMCpartons.Associate(ttHpartons, goodJets));
    for ( mapMCpJet::const_iterator mapIt = mapPtr->begin(); mapIt != mapPtr->end(); ++mapIt ) {
      // cout << "MC particle = " << mapIt->first->pid() << " of pt = " << mapIt->first->pt() << " associated to jet of et = " << mapIt->second->et() << endl;
      // Resolution histogram for jets associated to quarks
      jetVsMCpEt_->Fill(mapIt->first->pt(), mapIt->second->et());

      // -------------------------------------------------------------- //
      // -- THIS IS TEMPORARY, A MORE ACCURATE TAGGER SHOULD BE USED -- //
      // -------------------------------------------------------------- //
      // Consider as tagged those jets with highEff > 5.3.
      // high eff -> 50.30% b / 10.77% c / 0.92% uds /  0.98% g / 0.96% udsg // P.Schilling 23/10/07
      // Set the b-tag cut value
      float medium = 5.3;
      if ( mapIt->second->discriminatorHighEff()>medium ) {
        
        const OfflineJet * jet1 = mapIt->second;
        // Start from the next jet
        mapMCpJet::const_iterator subMapIt = mapIt;
        ++subMapIt;
        for ( ; subMapIt != mapPtr->end(); ++subMapIt ) {

          if ( subMapIt->second->discriminatorHighEff()>medium ) {

            const OfflineJet * jet2 = subMapIt->second;
            // Jet sum of the pair (ATTENTION to mass)
            BaseJet summedJet( sumJets( jet1, jet2 ) );

            const MCParticle * parton1 = mapIt->first;
            const MCParticle * parton2 = subMapIt->first;

            // DeltaR of the pair
            double dR = DeltaR(jet1->eta(), jet1->phi(),jet2->eta(), jet2->phi());
            double pairEt = summedJet.et();
            double pairEta = summedJet.eta();

            // Determine bin index
            int etId = int(pairEt/higgsEtBinSize_);
            int etaId = int(fabs(pairEta)/higgsEtaBinSize_);
            int dRId = int(dR/higgsDRbinSize_);
            if ( etId < 0 || etaId < 0 || dRId < 0 ) cout << "Error: index < 0, will crash..." << endl;
            if ( etId>int(higgsEtBinNum_-1) ) etId = higgsEtBinNum_-1;
            if ( etaId>int(higgsEtaBinNum_-1) ) etaId = higgsEtaBinNum_-1;
            if ( dRId>int(higgsDRbinNum_-1) ) dRId = higgsDRbinNum_-1;

            // If both b-tagged jets are associated bb from H fill HcombTrue.txt
            if ( abs(parton1->pid()) == 5 && abs(parton2->pid()) == 5 && parton1->mPid() == 25 && parton2->mPid() == 25 ) {
              // Evaluate the mass. Must do this by hand because of how E in BaseJet is evaluated.
              // (E1+E2)^2-P^2
              // Only if it is a true Higgs b-jet pair
              double mass = sqrt(pow(jet1->e()+jet2->e(),2) - pow(summedJet.e(),2));
              higgsMassTrue_->Fill(mass);
              trueHiggsPairEt_->Fill(pairEt);
              trueHiggsPairEta_->Fill(pairEta);
              trueHiggsPairDR_->Fill(dR);
              (trueH_[etId][etaId][dRId])++;
              // Save the true Higgs b-jet pair for the following top triplet
              trueHiggsPair = auto_ptr<BaseJet>(new BaseJet(summedJet));
            }
            // else fill HcombFalse.txt
            else {
              falseHiggsPairEt_->Fill(pairEt);
              falseHiggsPairEta_->Fill(pairEta);
              falseHiggsPairDR_->Fill(dR);
              (falseH_[etId][etaId][dRId])++;
            }
          }
        }
      }
    } // end of loop on b-pair combinations


    // Need to evaluate the DPhi between the hadronic top and the true Higgs.
    // Search for the top only if the trueHiggs was found
    if (trueHiggsPair.get() != 0) {
      for ( mapMCpJet::const_iterator mapIt = mapPtr->begin(); mapIt != mapPtr->end(); ++mapIt ) {
        const MCParticle * parton1 = mapIt->first;

        mapMCpJet::const_iterator subMapIt = mapIt;
        ++subMapIt;
        for ( ; subMapIt != mapPtr->end(); ++subMapIt ) {
          const MCParticle * parton2 = subMapIt->first;
          // Take the third jet, starting from the next to the second
          BaseJet summedJet( sumJets( mapIt->second, subMapIt->second ) );

          mapMCpJet::const_iterator subSubMapIt = subMapIt;
          ++subSubMapIt;
          for ( ; subSubMapIt != mapPtr->end(); ++subSubMapIt ) {
            const MCParticle * parton3 = subSubMapIt->first;
            const OfflineJet * jet3 = subSubMapIt->second;
            BaseJet tripletJet( sumJets( &summedJet, jet3 ) );

            double etHadronicTop = tripletJet.et();
            double etaHadronicTop = tripletJet.eta();
            double dPhiHiggsHadronicTop = DeltaPhi(trueHiggsPair->phi(), tripletJet.phi());
            // Determine bin index
            int etId = int(etHadronicTop/hadronicTopEtBinSize_);
            int etaId = int(fabs(etaHadronicTop)/hadronicTopEtaBinSize_);
            int dPhiHiggsHadronicTopId = int(dPhiHiggsHadronicTop/dPhiHiggsHadronicTopBinSize_);
            if ( etId < 0 || etaId < 0 || dPhiHiggsHadronicTopId < 0 ) cout << "Error: index < 0, will crash..." << endl;
            if ( etId>int(hadronicTopEtBinNum_-1) ) etId = hadronicTopEtBinNum_-1;
            if ( etaId>int(hadronicTopEtaBinNum_-1) ) etaId = hadronicTopEtaBinNum_-1;
            if ( dPhiHiggsHadronicTopId>int(dPhiHiggsHadronicTopBinNum_-1) ) dPhiHiggsHadronicTopId = dPhiHiggsHadronicTopId-1;

            // if all three jets are associated to partons from the same top
            // mPid = 6+24+24 = 54 or mPid = -6-24-24 = -54
            if ( abs(parton1->mPid() + parton2->mPid() + parton3->mPid()) == 54 ) {
              double hadronicTopMass = sqrt(pow((mapIt->second)->e()+(subMapIt->second)->e()+jet3->e(),2) - pow(tripletJet.e(),2));
              hadronicTopMassTrue_->Fill(hadronicTopMass);
              trueHadronicTopTripletEt_->Fill(etHadronicTop);
              trueHadronicTopTripletEta_->Fill(etaHadronicTop);
              trueHadronicTopTripletDphiHiggsHadronicTop_->Fill(dPhiHiggsHadronicTop);
              (trueHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId])++;
            }
            else {
              falseHadronicTopTripletEt_->Fill(etHadronicTop);
              falseHadronicTopTripletEta_->Fill(etaHadronicTop);
              falseHadronicTopTripletDphiHiggsHadronicTop_->Fill(dPhiHiggsHadronicTop);
              (falseHadronicTop_[etId][etaId][dPhiHiggsHadronicTopId])++;
            }
          }
        }
      }
    } // end if (true Higgs b-jet pair was found)
  } // end if( (H->cc || H->bb) && 4j )

  // Write the trueH_ and falseH_ to a txt file
  // ------------------------------------------
  ofstream higgsFile("HiggsPairProbability.txt");
  // The first line has informations on bin numbers and sizes
  higgsFile <<  "etBinNum = "  << higgsEtBinNum_  << " etBinSize = "  << higgsEtBinSize_
            << " etaBinNum = " << higgsEtaBinNum_ << " etaBinSize = " << higgsEtaBinSize_
            << " dRbinNum = "  << higgsDRbinNum_  << " dRbinSize = "  << higgsDRbinSize_ << endl;
  // The following lines have the counts for trueHiggs and falseHiggs pairs
  for(unsigned int i=0; i != higgsEtBinNum_; ++i) {
    for(unsigned int j=0; j != higgsEtaBinNum_; ++j) {
      for(unsigned int k=0; k != higgsDRbinNum_; ++k) {
        higgsFile <<  "trueHiggs["<<i<<"]["<<j<<"]["<<k<<"] = "  << trueH_[i][j][k]
                  << " falseHiggs["<<i<<"]["<<j<<"]["<<k<<"] = " << falseH_[i][j][k] << endl;
      }
    }
  }
  higgsFile.close();

  // Write the trueHadronicTop_ and falseHadronicTop_ to a txt file
  // --------------------------------------------------------------
  ofstream hadronicTopFile("HadronicTopTripletProbability.txt");
  // The first line has informations on bin numbers and sizes
  hadronicTopFile <<  "etBinNum = "                   << hadronicTopEtBinNum_        << " etBinSize = "                   << hadronicTopEtBinSize_
                  << " etaBinNum = "                  << hadronicTopEtaBinNum_       << " etaBinSize = "                  << hadronicTopEtaBinSize_
                  << " dPhiHiggsHadronicTopBinNum = " << dPhiHiggsHadronicTopBinNum_ << " dPhiHiggsHadronicTopBinSize = " << dPhiHiggsHadronicTopBinSize_ << endl;
  // The following lines have the counts for trueHadronicTop and falsHadronicTop pairs
  for(unsigned int i=0; i != hadronicTopEtBinNum_; ++i) {
    for(unsigned int j=0; j != hadronicTopEtaBinNum_; ++j) {
      for(unsigned int k=0; k != dPhiHiggsHadronicTopBinNum_; ++k) {
        hadronicTopFile <<  "trueHadronicTop["<<i<<"]["<<j<<"]["<<k<<"] = "  << trueHadronicTop_[i][j][k]
                        << " falseHadronicTop["<<i<<"]["<<j<<"]["<<k<<"] = " << falseHadronicTop_[i][j][k] << endl;
      }
    }
  }
  hadronicTopFile.close();
}

//       method called once each job just before starting event loop  
// -------------------------------------------------------------------------
void HadTopAndHiggsMassProbability::beginJob(const edm::EventSetup&) {
}


//       method called once each job just after ending the event loop 
// -------------------------------------------------------------------------
void HadTopAndHiggsMassProbability::endJob() {
  countTTHdecays_->writeDecays();
  delete countTTHdecays_;

  // delete the multidimensional arrays
  for(unsigned int i=0; i != higgsEtBinNum_; ++i) {
    for(unsigned int j=0; j != higgsEtaBinNum_; ++j) {
      delete[] trueH_[i][j];
      delete[] falseH_[i][j];
    }
    delete[] trueH_[i];
    delete[] falseH_[i];
  }

  for(unsigned int i=0; i != hadronicTopEtBinNum_; ++i) {
    for(unsigned int j=0; j != hadronicTopEtaBinNum_; ++j) {
      delete[] trueHadronicTop_[i][j];
      delete[] falseHadronicTop_[i][j];
    }
    delete[] trueHadronicTop_[i];
    delete[] falseHadronicTop_[i];
  }

  // Save histograms
  jetVsMCpEt_->Write();
  higgsMassTrue_->Write();

  trueHiggsPairEt_->Write();
  trueHiggsPairEta_->Write();
  trueHiggsPairDR_->Write();
  falseHiggsPairEt_->Write();
  falseHiggsPairEta_->Write();
  falseHiggsPairDR_->Write();

  trueHadronicTopTripletEt_->Write();
  trueHadronicTopTripletEta_->Write();
  trueHadronicTopTripletDphiHiggsHadronicTop_->Write();
  falseHadronicTopTripletEt_->Write();
  falseHadronicTopTripletEta_->Write();
  falseHadronicTopTripletDphiHiggsHadronicTop_->Write();

  outputFile_->Write();
}

// Define this as a plug-in
// ------------------------
// Also the line:
// <flags EDM_PLUGIN=1>
// should be added to the BuildFile before the export section
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HadTopAndHiggsMassProbability);

#endif // HADTOPANDHIGGSMASSPROBABILITY_CC
