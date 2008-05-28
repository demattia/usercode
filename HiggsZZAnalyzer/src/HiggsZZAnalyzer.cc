//
// Original Author:  Mia Tosi
//         Created:  Fri Feb 22 17:56:22 CET 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisExamples/HiggsZZAnalyzer/interface/HiggsZZAnalyzer.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"

#include "AnalysisExamples/AnalysisClasses/interface/DeltaPhi.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"
//#include "AnalysisExamples/AnalysisClasses/interface/DiMass.h"
#include "AnalysisExamples/AnalysisClasses/interface/DiParticleMass.h"

//#include "AnalysisExamples/AnalysisObjects/interface/MathParameters.h"
#include "AnalysisExamples/AnalysisObjects/interface/PythiaParticleIndex.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesMass.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesCharge.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>


HiggsZZAnalyzer::HiggsZZAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  offlineJetLabel_(            iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets"                  ) ),
  offlineMEtLabel_(            iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt"                   ) ),
  globalMuonLabel_(            iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons"                  ) ),
  simpleElectronLabel_(        iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons"              ) ),
  simpleTauLabel_(             iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus"                   ) ),
  combinedSVBJetTagsLabel_(    iConfig.getUntrackedParameter<edm::InputTag>( "combinedSVBJetTags"           ) ),
  leptonPtCut_(	               iConfig.getUntrackedParameter<double>(        "leptonPtCut"                  ) ),
  leptonIPsignificanceCut_(    iConfig.getUntrackedParameter<double>(        "leptonIPsignificanceCut"      ) ),
  muonIso03SumPtCut_(          iConfig.getUntrackedParameter<double>(        "muonIso03SumPtCut"            ) ),
  muonIso03emEtCut_(           iConfig.getUntrackedParameter<double>(        "muonIso03emEtCut"             ) ),
  centralJetEtCut_(            iConfig.getUntrackedParameter<double>(        "centralJetEtCut"              ) ),
  centralJetEtaCut_(           iConfig.getUntrackedParameter<double>(        "centralJetEtaCut"             ) ),
  centralJetLeptonDeltaRCut_(  iConfig.getUntrackedParameter<double>(        "centralJetLeptonDeltaRCut"    ) ),
  centralJetEMfracCut_(        iConfig.getUntrackedParameter<double>(        "centralJetEMfracCut"          ) ),
  forwardJetEtCut_(            iConfig.getUntrackedParameter<double>(        "forwardJetEtCut"              ) ),
  forwardJetEtaCut_(           iConfig.getUntrackedParameter<double>(        "forwardJetEtaCut"             ) ),
  forwardJetLeptonDeltaRCut_(  iConfig.getUntrackedParameter<double>(        "forwardJetLeptonDeltaRCut"    ) ),
  forwardDiJetMassCut_(        iConfig.getUntrackedParameter<double>(        "forwardDiJetMassCut"          ) ),
  pTbalanceCut_(               iConfig.getUntrackedParameter<double>(        "pTbalanceCut"                 ) ),
  ZdileptonMassminCut_(        iConfig.getUntrackedParameter<double>(        "ZdileptonMassminCut"          ) ),
  ZdileptonMassmaxCut_(        iConfig.getUntrackedParameter<double>(        "ZdileptonMassmaxCut"          ) ),
  ZdileptonDeltaRminCut_(      iConfig.getUntrackedParameter<double>(        "ZdileptonDeltaRminCut"        ) ),
  ZdileptonDeltaRmaxCut_(      iConfig.getUntrackedParameter<double>(        "ZdileptonDeltaRmaxCut"        ) ),
  ZdileptonShiftedEtaCut_(     iConfig.getUntrackedParameter<double>(        "ZdileptonShiftedEtaCut"       ) ),
  ZdijetMassminCut_(           iConfig.getUntrackedParameter<double>(        "ZdijetMassminCut"             ) ),
  ZdijetMassmaxCut_(           iConfig.getUntrackedParameter<double>(        "ZdijetMassmaxCut"             ) ),
  ZdijetDeltaRminCut_(         iConfig.getUntrackedParameter<double>(        "ZdijetDeltaRminCut"           ) ),
  ZdijetDeltaRmaxCut_(         iConfig.getUntrackedParameter<double>(        "ZdijetDeltaRmaxCut"           ) ),
  ZdijetShiftedEtaCut_(        iConfig.getUntrackedParameter<double>(        "ZdijetShiftedEtaCut"          ) ),
  ZZdeltaRminCut_(             iConfig.getUntrackedParameter<double>(        "ZZdeltaRminCut"               ) ),
  ZZdeltaRmaxCut_(             iConfig.getUntrackedParameter<double>(        "ZZdeltaRmaxCut"               ) ),
  MCParticleLabel_(            iConfig.getUntrackedParameter<std::string>(   "MCParticles"                  ) ),
  simVtxLabel_(                iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx"                       ) ),
  HMassCut_(                   iConfig.getUntrackedParameter<double>(        "HMassCut"                     ) )
{
  // Now do what ever initialization is needed
  // -----------------------------------------

  //
  // constants, enums and typedefs

  nZ_           = 2;
  nZleptons_    = 2;
  nZjets_       = 2;
  nforwardjets_ = 2;
  njets_        = 4;

  eventcounter_              = 0;
  signalTopologySizeCounter_ = 0;
  HtoZZeventcounter_ = 0;
  muonCounter_              = 0;
  goodMuonCounter_          = 0;
  electronCounter_          = 0;
  goodElectronCounter_      = 0;
  goodIc5JetCounter_        = 0;
  goodIc5centralJetCounter_ = 0;
  goodIc5forwardJetCounter_ = 0;
  Ic5forwardJetCounter_     = 0;

  nbin_         = 100;
  first_bin_et  =   0.;
  last_bin_et   = 200.;
  first_bin_eta =  -3.;
  last_bin_eta  =   3.;

}


HiggsZZAnalyzer::~HiggsZZAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // write out the file
  // ------------------
  OutputFile->Write();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HiggsZZAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace anaobj;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // Global muons
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;
  iEvent.getByLabel( globalMuonLabel_, globalMuons );
  // SimpleElectrons
  // ---------------
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  iEvent.getByLabel( simpleElectronLabel_, simpleElectrons);
  // Calorimeter jets
  // ----------------
  edm::Handle<OfflineJetCollection> offlineJets;
  iEvent.getByLabel( offlineJetLabel_, offlineJets );
  // Missing Et
  // ----------
  edm::Handle<OfflineMEt> MEt;
  iEvent.getByLabel( offlineMEtLabel_, MEt );
  // SimpleTaus
  // ----------
  edm::Handle < SimpleTauCollection > simpleTaus;
  iEvent.getByLabel( simpleTauLabel_, simpleTaus );
  // MCParticle
  // ----------
  edm::Handle < CandidateCollection > MCparticles;
  iEvent.getByLabel( MCParticleLabel_, MCparticles );

  unsigned globalMuonsSize     = globalMuons->size();
  unsigned simpleElectronsSize = simpleElectrons->size();
  unsigned offlineJetsSize     = offlineJets->size();
  unsigned MCparticlesSize     = MCparticles->size();
  unsigned simpleTausSize      = simpleTaus->size();
  if ( ( globalMuonsSize >= nZleptons_ || simpleElectronsSize >= nZleptons_ ) &&
       offlineJetsSize >= njets_ ) {
    // count number of events w/ at least 2 leptons and 4 jets
    // as we expect in H->ZZ->lljj via VBF
    ++signalTopologySizeCounter_;

    /////////////////////////////////// HEPG analysis ////////////////////////////
    std::vector< const Candidate* > Z1Candidates;
    std::vector< const Candidate* > Z2Candidates;
    std::vector< const Candidate* > HCandidates;
    std::vector< const Candidate* > bbCandidates;
    Particle * Z1particle = 0;
    Particle * Z2particle = 0;
    Particle * Hparticle  = 0;
    const Candidate* tmpMotherZ = 0;
    
    // Take CandidateCollection
    // ------------------------
    int NHcandidates  = 0;
    int NZ1candidates = 0;
    int NZ2candidates = 0;
    int Nbbcandidates = 0;
    // take Z decay particles
    // and save their ptr into a vector
    // --------------------------------
    CandidateCollection::const_iterator MCparticle = MCparticles->begin();
    for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
      int tmpParticleId = MCparticle->pdgId();
      if ( fabs(tmpParticleId) > pythiaH_ ) continue;
      int tmpMotherId   = 0;
      if( MCparticle->mother() != 0 ) tmpMotherId = MCparticle->mother()->pdgId();

      // H decay
      if ( tmpMotherId == pythiaH_ && tmpParticleId != tmpMotherId ) {
	HCandidates.push_back(&*MCparticle);
	NHcandidates++;
      }
      // Z decay
      if ( tmpMotherId == pythiaZ_ && tmpParticleId != tmpMotherId ) {
	if ( tmpMotherZ == 0 ) {
	  Z1Candidates.push_back(&*MCparticle);
	  tmpMotherZ = MCparticle->mother();
	  NZ1candidates++;
	} else {
	  if ( MCparticle->mother() == tmpMotherZ ) {
	    Z1Candidates.push_back(&*MCparticle);
	    NZ1candidates++;
	  } else {
	    Z2Candidates.push_back(&*MCparticle);
	    NZ2candidates++;
	  }
	}
      }
      if (tmpMotherId != pythiaH_ && tmpMotherId != pythiaZ_ ) {
	if ( fabs(tmpParticleId) == pythiab_ && tmpParticleId != tmpMotherId && tmpMotherId != 0 ) {
	  if ( Nbbcandidates == 2 ) continue;
	  bbCandidates.push_back(&*MCparticle);
	  Nbbcandidates++;
	}
      }
    }

    // manage info of particles from Z decay
    // -------------------------------------
    int Z1particlesId = 0;
    double Z1particle1Pt       = 0.;
    double Z1particle2Pt       = 0.;
    double Z1particle1Eta      = 0.;
    double Z1particle2Phi      = 0.;
    double Z1particlesDeltaPhi = 0.;
    double Z1particlesDeltaEta = 0.;
    double Z1particlesDeltaR   = 0.;
    double Z1mass = 0.;
    int Z2particlesId = 0;
    double Z2particle1Pt       = 0.;
    double Z2particle2Pt       = 0.;
    double Z2particle1Eta      = 0.;
    double Z2particle2Phi      = 0.;
    double Z2particlesDeltaPhi = 0.;
    double Z2particlesDeltaEta = 0.;
    double Z2particlesDeltaR   = 0.;
    double Z2mass = 0.;
    if ( Z1Candidates.size() == 2 ) {
      Z1particlesId       = fabs(int(Z1Candidates[0]->pdgId()));
      Z1particle1Pt       = Z1Candidates[0]->pt();
      Z1particle2Pt       = Z1Candidates[1]->pt();
      Z1particle1Eta      = Z1Candidates[0]->eta();
      Z1particle2Phi      = Z1Candidates[1]->phi();
      int Z1particle1Charge   = Z1Candidates[0]->charge();
      int Z1particle2Charge   = Z1Candidates[1]->charge();
      Z1particlesDeltaPhi = DeltaPhi(Z1Candidates[0],Z1Candidates[1]);
      Z1particlesDeltaEta = Z1Candidates[0]->eta()-Z1Candidates[1]->eta();
      Z1particlesDeltaR   = TMath::Sqrt(pow(Z1particlesDeltaEta,2)+
					pow(Z1particlesDeltaPhi,2));
      Z1particle = new Particle(Z1particle1Charge+Z1particle2Charge,
			        Z1Candidates[0]->p4()+
                               (Z1Candidates[1]->p4()),
				Z1Candidates[0]->vertex(),
				pythiaZ_,0,true);
      Z1mass = Z1particle->mass();
      if (Nbbcandidates == 2 ) {
	double bbMass = -99.;
	double m4     = -99.;
	bbMass = DiParticleMass<Candidate>(bbCandidates[0],bbCandidates[1]);
	Mbb_->Fill(bbMass); 
	Particle::LorentzVector particles4 = Particle::LorentzVector(bbCandidates[0]->p4()+
								     (bbCandidates[1]->p4())+
								     (Z1particle->p4()));
	m4 = particles4.mass();
	M4_->Fill(m4);
	MZvsMbb_->Fill(Z1mass,bbMass);
	if (bbMass  >= ZdijetMassminCut_    && bbMass <= ZdijetMassmaxCut_ &&
	    Z1mass  >= ZdileptonMassminCut_ && Z1mass <= ZdileptonMassmaxCut_ )
	  M4irrBkg_->Fill(m4);
      }
    }
    if ( Z2Candidates.size() == 2 ) { 
      Z2particlesId       = fabs(int(Z2Candidates[0]->pdgId()));
      Z2particle1Pt       = Z2Candidates[0]->pt();
      Z2particle2Pt       = Z2Candidates[1]->pt();
      Z2particle1Eta      = Z2Candidates[0]->eta();
      Z2particle2Phi      = Z2Candidates[1]->phi();
      int Z2particle1Charge   = Z2Candidates[0]->charge();
      int Z2particle2Charge   = Z2Candidates[1]->charge();
      Z2particlesDeltaPhi = DeltaPhi(Z2Candidates[0],Z2Candidates[1]);
      Z2particlesDeltaEta = Z2Candidates[0]->eta()-Z2Candidates[1]->eta();
      Z2particlesDeltaR   = TMath::Sqrt(pow(Z2particlesDeltaEta,2)+
					pow(Z2particlesDeltaPhi,2));
      Z2particle = new Particle(Z2particle1Charge+Z2particle2Charge,
				Z2Candidates[0]->p4()+
		  	       (Z2Candidates[1]->p4()),
				Z2Candidates[0]->vertex(),
				pythiaZ_,0,true);
      Z2mass = Z2particle->mass();
    }
	  
    std::vector< std::pair<int,int> >::const_iterator vecTypeIndexPair_itr = vecTypeIndexPair.begin();
    for ( ; vecTypeIndexPair_itr != vecTypeIndexPair.end(); ++vecTypeIndexPair_itr ) {
      int type  = vecTypeIndexPair_itr->first;
      int index = vecTypeIndexPair_itr->second;
      if ( Z1particlesId == type ) {
	vecZtoParticleseventcounter_[index]++;
	vecMZtoParticles_[index].Fill(Z1mass);
	vecPtZParticles_[index].Fill(Z1particle1Pt);
	vecPtZParticles_[index].Fill(Z1particle2Pt);
	vecPtZParticle1_[index].Fill(Z1particle1Pt);
	vecPtZParticle2_[index].Fill(Z1particle2Pt);
	vecDeltaPhiZParticles_[index].Fill(Z1particlesDeltaPhi);
	vecDeltaEtaZParticles_[index].Fill(Z1particlesDeltaEta);
	vecDeltaRZParticles_[index].Fill(Z1particlesDeltaR);
      }
      if ( Z2particlesId == type ) {
	vecZtoParticleseventcounter_[index]++;
	vecMZtoParticles_[index].Fill(Z2mass);
	vecPtZParticles_[index].Fill(Z2particle1Pt);
	vecPtZParticles_[index].Fill(Z2particle2Pt);
	vecPtZParticle1_[index].Fill(Z2particle1Pt);
	vecPtZParticle2_[index].Fill(Z2particle2Pt);
	vecDeltaPhiZParticles_[index].Fill(Z2particlesDeltaPhi);
	vecDeltaEtaZParticles_[index].Fill(Z2particlesDeltaEta);
	vecDeltaRZParticles_[index].Fill(Z2particlesDeltaR);
      }
    }
    if (Z1particlesId <= pythiat_ ){
      vecZtoParticleseventcounter_[3]++;
      vecDeltaPhiZParticles_[3].Fill(Z1particlesDeltaPhi);
      vecDeltaEtaZParticles_[3].Fill(Z1particlesDeltaEta);
      vecDeltaRZParticles_[3].Fill(Z1particlesDeltaR);
    }
    if (Z2particlesId <= pythiat_) {
      vecZtoParticleseventcounter_[3]++;
      vecDeltaPhiZParticles_[3].Fill(Z2particlesDeltaPhi);
      vecDeltaEtaZParticles_[3].Fill(Z2particlesDeltaEta);
      vecDeltaRZParticles_[3].Fill(Z2particlesDeltaR);
    }
    else if (Z1particlesId >= pythiae_ && Z1particlesId <= pythianutau_) 
      vecZtoParticleseventcounter_[4]++;
    if (Z2particlesId >= pythiae_ && Z2particlesId <= pythianutau_) 
      vecZtoParticleseventcounter_[4]++;
	  
    if ( Z1particle && Z2particle ) {
      double ZZDeltaPhi = -99.;
      double ZZDeltaEta = -99.;
      double ZZDeltaR   = -99.;
      double ZZdimass   = -99.;
      ZZDeltaPhi = DeltaPhi(Z1particle,Z2particle);
      ZZDeltaEta = (Z1particle->p4()).eta()-(Z2particle->p4()).eta();
      ZZDeltaR   = TMath::Sqrt(pow(ZZDeltaEta,2)+
			       pow(ZZDeltaPhi,2));
      Particle::LorentzVector ZZsystem = Particle::LorentzVector(Z1particle->p4()+
								 (Z2particle->p4()));
      ZZdimass = ZZsystem.mass();

      DeltaPhiZZ_->Fill(ZZDeltaPhi);
      DeltaEtaZZ_->Fill(ZZDeltaEta);
      DeltaRZZ_->Fill(ZZDeltaR);
      if (((fabs(Z1particlesId) == pythiab_ && (fabs(Z2particlesId) == pythiae_ || fabs(Z2particlesId) == pythiamu_)) ||
	   (fabs(Z2particlesId) == pythiab_ && (fabs(Z1particlesId) == pythiae_ || fabs(Z1particlesId) == pythiamu_)) )  ) 
	MZZ_->Fill(ZZdimass);
      MZ1vsMZ2_->Fill(Z1mass,Z2mass);
    }

    // manage info of particles from H decay
    // -------------------------------------
    int Hparticle1DaughterId = 0;
    int Hparticle2DaughterId = 0;
    int HparticlesId         = 0;
    double Hparticle1Pt       = 0.;
    double Hparticle2Pt       = 0.;
    double Hparticle1Eta      = 0.;
    double Hparticle1Phi      = 0.;
    double Hparticle2Eta      = 0.;
    double Hparticle2Phi      = 0.;
    double Hparticle1Mass     = 0.;
    double Hparticle2Mass     = 0.;
    double Hmass = 0.;
    if (HCandidates.size() != 0 ) {
      Hparticle1DaughterId = HCandidates[0]->daughter(0)->pdgId();
      Hparticle2DaughterId = HCandidates[1]->daughter(0)->pdgId();
      HparticlesId         = fabs(int(HCandidates[0]->pdgId()));
      Hparticle1Pt         = HCandidates[0]->pt();
      Hparticle2Pt         = HCandidates[1]->pt();
      Hparticle1Eta        = HCandidates[0]->eta();
      Hparticle1Phi        = HCandidates[0]->phi();
      Hparticle2Eta        = HCandidates[1]->eta();
      Hparticle2Phi        = HCandidates[1]->phi();
      Hparticle1Mass       = HCandidates[0]->mass();
      Hparticle2Mass       = HCandidates[1]->mass();
      Hparticle = new Particle(HCharge_,
			       Z1particle->p4()+
			      (Z2particle->p4()),
			       Z1particle->vertex(),
			       pythiaH_,0,true);
      Hmass = Hparticle->mass();

      // Fill histograms
      // ---------------
      if ( HparticlesId == pythiaZ_ ) {
	++HtoZZeventcounter_;
	PtHParticles_->Fill(Hparticle1Pt);
	PtHParticles_->Fill(Hparticle2Pt);
	PtHParticle1_->Fill(Hparticle1Pt);
	PtHParticle2_->Fill(Hparticle2Pt);
	if ( Z1mass <= 60. && Z2mass <= 60. ) std::cout << "PYTHIA BUG" << std::endl;
	else {
	  if (((fabs(Hparticle1DaughterId) <= pythiat_ && (fabs(Hparticle2DaughterId) == pythiae_ || fabs(Hparticle2DaughterId) == pythiamu_)) ||
	       (fabs(Hparticle2DaughterId) <= pythiat_ && (fabs(Hparticle1DaughterId) == pythiae_ || fabs(Hparticle1DaughterId) == pythiamu_)) )  
	      ) {
	    MHtoZZllqq_->Fill(Hmass);
	    if (fabs(Hparticle1DaughterId) == pythiab_ ) MHtoZZllbb_->Fill(Hmass);
	    if (fabs(Hparticle2DaughterId) == pythiab_ ) MHtoZZllbb_->Fill(Hmass);
	  }   	  
	}
      }
    }
  
    /////////////////////////////////// OFFLINE analysis ////////////////////////////
    std::vector<GlobalMuon>     goodMuonVec;
    std::vector<SimpleElectron> goodElectronVec;
    std::vector<OfflineJet>     goodIc5centralJetVec;
    std::vector<OfflineJet>     goodIc5forwardJetVec;
    // Global muons
    // ------------
    MuonNumber_->Fill(globalMuons->size());
    GlobalMuonCollection::const_iterator globalMuon_itr = globalMuons->begin();
    for ( ; globalMuon_itr != globalMuons->end(); ++globalMuon_itr ) {
      ++muonCounter_;
      double tmpMuonEta                = globalMuon_itr->eta();
      double tmpMuonPhi                = globalMuon_itr->phi();
      double tmpMuonPt                 = globalMuon_itr->pt();
      double tmpMuonPx                 = globalMuon_itr->px();
      double tmpMuonPy                 = globalMuon_itr->py();
      double tmpMuonPz                 = globalMuon_itr->pz();
      double tmpMuonEMcalo             = globalMuon_itr->caloEm();
      double tmpMuonHADcalo            = globalMuon_itr->caloHad();
      double tmpMuonHOcalo             = globalMuon_itr->caloHo();
      double tmpMuonIso03SumPt         = globalMuon_itr->iso03SumPt();
      double tmpMuonIso03emEt          = globalMuon_itr->iso03emEt();
      double tmpMuonIso03haEt          = globalMuon_itr->iso03hadEt();
      double tmpMuonIso03hoEt          = globalMuon_itr->iso03hoEt();
      int    tmpMuonIso03nJets         = globalMuon_itr->iso03nJets();
      int    tmpMuonIso03nTracks       = globalMuon_itr->iso03nTracks();
      double tmpMuonNormChi2           = globalMuon_itr->normChi2();
      double tmpMuonHitsNum            = globalMuon_itr->hitsNum();
      double tmpMuonImpactParamXY      = 0.;
      double tmpMuonErrorImpactParamXY = 0.;
      double tmpMuonImpactParamSig     = 99.;
      tmpMuonImpactParamXY             = globalMuon_itr->impactParamXY();
      tmpMuonErrorImpactParamXY        = globalMuon_itr->errorImpactParamXY();
      if(tmpMuonErrorImpactParamXY != 0.) tmpMuonImpactParamSig = tmpMuonImpactParamXY/tmpMuonErrorImpactParamXY;
      // they seem to be equal to 0!!
      //-----------------------------
      //      std::cout << "tmpMuonIso03SumPt:" << tmpMuonIso03SumPt << std::endl;
      //      std::cout << "tmpMuonIso03emEt:"  << tmpMuonIso03emEt  << std::endl;

      if ( tmpMuonPt > leptonPtCut_ 
	   && tmpMuonImpactParamSig < leptonIPsignificanceCut_ ) {
	goodMuonVec.push_back(*globalMuon_itr);
	++goodMuonCounter_;
	if ( tmpMuonIso03SumPt != 0. && tmpMuonIso03SumPt <= muonIso03emEtCut_ ) std::cout << "good SumPt" << std::endl;
	if ( tmpMuonIso03emEt  != 0. && tmpMuonIso03emEt  <= muonIso03emEtCut_ ) std::cout << "good emEt"  << std::endl;
      }

      if (tmpMuonEta        != 0.) MuonEta_->Fill(       tmpMuonEta       );
      if (tmpMuonPhi        != 0.) MuonPhi_->Fill(       tmpMuonPhi       );
      if (tmpMuonPt         != 0.) MuonPt_->Fill(        tmpMuonPt        );
      if (tmpMuonPx         != 0.) MuonPx_->Fill(        tmpMuonPx        );
      if (tmpMuonPy         != 0.) MuonPy_->Fill(        tmpMuonPy        );
      if (tmpMuonPz         != 0.) MuonPz_->Fill(        tmpMuonPz        );
      if (tmpMuonEMcalo     != 0.) MuonEMcalo_->Fill(    tmpMuonEMcalo    );
      if (tmpMuonHADcalo    != 0.) MuonHADcalo_->Fill(   tmpMuonHADcalo   );
      if (tmpMuonHOcalo     != 0.) MuonHOcalo_->Fill(    tmpMuonHOcalo    );
      if (tmpMuonIso03SumPt != 0.) MuonIso03SumPt_->Fill(tmpMuonIso03SumPt);
      if (tmpMuonIso03emEt  != 0.) MuonIso03emEt_->Fill( tmpMuonIso03emEt);
      if (tmpMuonNormChi2   != 0.) MuonNormChi2_->Fill(  tmpMuonNormChi2  );
      if (tmpMuonHitsNum    != 0.) MuonHitsNum_->Fill(   tmpMuonHitsNum   );    
      MuonImpactParamXY_->Fill(tmpMuonImpactParamXY);
      MuonErrorImpactParamXY_->Fill(tmpMuonErrorImpactParamXY);
      MuonImpactParamSig_->Fill(tmpMuonImpactParamSig);
    }
    goodMuonNumber_->Fill(goodMuonVec.size());

    // SimpleElectrons
    // ---------------
    ElectronNumber_->Fill(simpleElectrons->size());
    SimpleElectronCollection::const_iterator simpleElectron_itr = simpleElectrons->begin();
    for ( ; simpleElectron_itr != simpleElectrons->end(); ++simpleElectron_itr ) {
      ++electronCounter_;
      double tmpElectronEt                 = simpleElectron_itr->et();
      double tmpElectronEta                = simpleElectron_itr->eta();
      double tmpElectronPhi                = simpleElectron_itr->phi();
      double tmpElectronPt                 = simpleElectron_itr->pt();
      double tmpElectronPx                 = simpleElectron_itr->px();
      double tmpElectronPy                 = simpleElectron_itr->py();
      double tmpElectronPz                 = simpleElectron_itr->pz();
      double tmpElectronHadOverEm          = simpleElectron_itr->hadOverEm();
      double tmpElectronIsoVal             = simpleElectron_itr->isoVal();
      double tmpElectronImpactParamXY      = 0.; 
      double tmpElectronErrorImpactParamXY = 0.;
      double tmpElectronImpactParamSig     = 99.;
      tmpElectronImpactParamXY             = simpleElectron_itr->impactParamXY();
      tmpElectronErrorImpactParamXY        = simpleElectron_itr->errorImpactParamXY();
      if ( tmpElectronErrorImpactParamXY != 0. ) tmpElectronImpactParamSig = tmpElectronImpactParamXY / tmpElectronErrorImpactParamXY;
      if ( tmpElectronPt > leptonPtCut_ 
	   && tmpElectronImpactParamSig < leptonIPsignificanceCut_ ) {
	goodElectronVec.push_back(*simpleElectron_itr);
	++goodElectronCounter_;
      }

      if ( tmpElectronEt        != 0. ) ElectronEt_->Fill(       tmpElectronEt       );
      if ( tmpElectronEta       != 0. ) ElectronEta_->Fill(      tmpElectronEta      );
      if ( tmpElectronPhi       != 0. ) ElectronPhi_->Fill(      tmpElectronPhi      );
      if ( tmpElectronPt        != 0. ) ElectronPt_->Fill(       tmpElectronPt       );
      if ( tmpElectronPx        != 0. ) ElectronPx_->Fill(       tmpElectronPx       );
      if ( tmpElectronPy        != 0. ) ElectronPy_->Fill(       tmpElectronPy       );
      if ( tmpElectronPz        != 0. ) ElectronPz_->Fill(       tmpElectronPz       );
      if ( tmpElectronHadOverEm != 0. ) ElectronHadOverEm_->Fill(tmpElectronHadOverEm);
      if ( tmpElectronIsoVal    != 0. ) ElectronIsoVal_->Fill(   tmpElectronIsoVal   );
      ElectronImpactParamXY_->Fill(tmpElectronImpactParamXY);
      ElectronErrorImpactParamXY_->Fill(tmpElectronErrorImpactParamXY);
      ElectronImpactParamSig_->Fill(tmpElectronImpactParamSig);
    }
    goodElectronNumber_->Fill(goodElectronVec.size());

    // Missing Et
    // ----------
    double tmpSumEt  = MEt->sumEt();
    double tmpMEt    = MEt->et();
    double tmpMEtSig = MEt->mEtSig();
    double tmpMEtPhi = MEt->phi();
    MEtSumEt_->Fill(       tmpSumEt );
    MEt_->Fill(            tmpMEt   );
    MEtSignificance_->Fill(tmpMEtSig);
    MEtPhi_->Fill(         tmpMEtPhi);
  
    double DiLeptonMass      = -99.;
    double leptonDeltaR      = -99.;
    double lepton1ShiftedEta = 99999.;
    double lepton2ShiftedEta = 99999.;
    double centralDiJetMass      = -99.;
    double centralJetDeltaR      = -99.;
    double centralJet1ShiftedEta = 99999.;
    double centralJet2ShiftedEta = 99999.;
    double pTbalance   = 99999.;
    double lljjInvMass = -99.;
    std::vector< Particle > ZleptonsVec;
    std::vector< Particle > ZjetsVec;
    std::vector< Particle > VBFsystemVec;
    Particle * leptonicZrec = 0;
    Particle * hadronicZrec = 0;
    Particle * systemVBFrec = 0;
    Particle * recH = 0;
    if ( goodMuonVec.size() >= nZleptons_ || goodElectronVec.size() >= nZleptons_ ) {
      // Calorimeter jets
      // ----------------
      OfflinejetNumber_->Fill(offlineJets->size());
      std::vector<OfflineJet> Ic5forwardJetVec;
      //      if ( offlineJets->size() != 0 ) { 
	OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
	  double tmpOfflinejetUncorrEt     = offlinejet->uncorrEt();
	  double tmpOfflinejetEt           = offlinejet->et();
	  double tmpOfflinejetEtScale      = tmpOfflinejetEt/tmpOfflinejetUncorrEt;
	  double tmpOfflinejetEMenergyfrac = offlinejet->emEnergyFraction();
	  double tmpOfflinejetHighEffDiscr = double(offlinejet->discriminatorHighEff());
	  double tmpOfflinejetHighPurDiscr = double(offlinejet->discriminatorHighPur());
	  double tmpOfflinejetMass         = offlinejet->jetMass();
	  int    tmpOfflinejetTKnumS2      = offlinejet->tkNumS2();
	  double tmpOfflinejetTKsumPtS2    = offlinejet->tkSumPtS2();
	  double tmpOfflinejetTKtagmassS2  = offlinejet->tagTkMassS2();
	  double tmpOfflinejetPhi          = offlinejet->phi();
	  double tmpOfflinejetEta          = offlinejet->eta();
	  double tmp_jetlepton_deltaR = -99.;
	  if ( goodMuonVec.size() >= nZleptons_ ) {
	    std::vector<GlobalMuon>::const_iterator goodMuonVec_itr = goodMuonVec.begin();
	    for ( ; goodMuonVec_itr != goodMuonVec.end(); ++goodMuonVec_itr ) {
	      double muonEta = goodMuonVec_itr->eta();
	      double muonPhi = goodMuonVec_itr->phi();
	      double tmp_jetlepton_deltaEta = tmpOfflinejetEta-muonEta;
	      double tmp_jetlepton_deltaPhi = DeltaPhi(tmpOfflinejetPhi,muonPhi);
	      tmp_jetlepton_deltaR          = TMath::Sqrt(pow(tmp_jetlepton_deltaEta,2)+
							  pow(tmp_jetlepton_deltaPhi,2));
	      if ( tmp_jetlepton_deltaPhi < centralJetLeptonDeltaRCut_ ) continue;
	    }
	  } else if ( goodElectronVec.size() >= nZleptons_ ) {
	    std::vector<SimpleElectron>::const_iterator goodElectronVec_itr = goodElectronVec.begin();
	    for( ; goodElectronVec_itr != goodElectronVec.end(); ++goodElectronVec_itr ) {
	      double electronEta = goodElectronVec_itr->eta();
	      double electronPhi = goodElectronVec_itr->phi();
	      double tmp_jetlepton_deltaEta = tmpOfflinejetEta-electronEta;
	      double tmp_jetlepton_deltaPhi = DeltaPhi(tmpOfflinejetPhi,electronPhi);
	      tmp_jetlepton_deltaR          = TMath::Sqrt(pow(tmp_jetlepton_deltaEta,2)+
							  pow(tmp_jetlepton_deltaPhi,2));
	      if ( tmp_jetlepton_deltaPhi < centralJetLeptonDeltaRCut_ ) continue;
	    }
	  }
	  if (tmp_jetlepton_deltaR != -99.) OfflinejetLeptonDeltaR_->Fill(tmp_jetlepton_deltaR);

	  // jets selection
	  // --------------
	  if ( tmpOfflinejetEt         > centralJetEtCut_     
	       && tmp_jetlepton_deltaR > centralJetLeptonDeltaRCut_ ) {
	    ++goodIc5JetCounter_;
	    // jets from Z
	    // -----------
	    if ( fabs(tmpOfflinejetEta)       < centralJetEtaCut_ 
		 //		 && tmpOfflinejetEMenergyfrac < centralJetEMfracCut_ 
		 ) {
	      ++goodIc5centralJetCounter_;
	      goodIc5centralJetVec.push_back(*offlinejet);
	      centralJetLeptonDeltaR_->Fill(tmp_jetlepton_deltaR);	      
	    }
	    // forward jets
	    // ------------
	    if ( fabs(tmpOfflinejetEta) > forwardJetEtaCut_ ) {
	      ++Ic5forwardJetCounter_;
	      Ic5forwardJetVec.push_back(*offlinejet);
	      forwardJetLeptonDeltaR_->Fill(tmp_jetlepton_deltaR);
	    }
	  }
	  goodIc5JetNumber_->Fill(        goodIc5JetCounter_        );
	  goodIc5centralJetNumber_->Fill( goodIc5centralJetVec.size() );
	  Ic5forwardJetNumber_->Fill(     Ic5forwardJetCounter_     );
	  
	  OfflinejetUncorrEt_->Fill(                                           tmpOfflinejetUncorrEt     );
	  OfflinejetEt_->Fill(                                                 tmpOfflinejetEt           );
	  OfflinejetEtScale_->Fill(                                            tmpOfflinejetEtScale      );          
	  OfflinejetEtScaleVSUncorrEt_->Fill(                                  tmpOfflinejetUncorrEt,tmpOfflinejetEtScale);          
	  OfflinejetEtScaleVSEt_->Fill(                                        tmpOfflinejetEt,      tmpOfflinejetEtScale);                
	  OfflinejetEtScaleVSEta_->Fill(                                       tmpOfflinejetEta,     tmpOfflinejetEtScale);               
	  OfflinejetEtVSUncorrEt_->Fill(                                       tmpOfflinejetUncorrEt,tmpOfflinejetEt);          
	  OfflinejetEtVSEt_->Fill(                                             tmpOfflinejetEt,      tmpOfflinejetEt);                
	  OfflinejetEtVSEta_->Fill(                                            tmpOfflinejetEta,     tmpOfflinejetEt);               
	  OfflinejetUncorrEtVSUncorrEt_->Fill(                                 tmpOfflinejetUncorrEt,tmpOfflinejetUncorrEt);          
	  OfflinejetUncorrEtVSEt_->Fill(                                       tmpOfflinejetEt,      tmpOfflinejetUncorrEt);                
	  OfflinejetUncorrEtVSEta_->Fill(                                      tmpOfflinejetEta,     tmpOfflinejetUncorrEt);               
	  if (tmpOfflinejetEMenergyfrac != 0. ) OfflinejetEMenergyfrac_->Fill( tmpOfflinejetEMenergyfrac );
	  if (tmpOfflinejetHighEffDiscr >= 0. ) OfflinejetHighEffDiscr_->Fill( tmpOfflinejetHighEffDiscr );
	  if (tmpOfflinejetHighPurDiscr >= 0. ) OfflinejetHighPurDiscr_->Fill( tmpOfflinejetHighPurDiscr );
	  if (tmpOfflinejetMass         != 0. ) OfflinejetMass_->Fill(         tmpOfflinejetMass         );
	  if (tmpOfflinejetTKnumS2      != 0. ) OfflinejetTKnumS2_->Fill(      tmpOfflinejetTKnumS2      );
	  if (tmpOfflinejetTKsumPtS2    != 0. ) OfflinejetTKsumPtS2_->Fill(    tmpOfflinejetTKsumPtS2    );
	  if (tmpOfflinejetTKtagmassS2  != 0. ) OfflinejetTKtagmassS2_->Fill(  tmpOfflinejetTKtagmassS2  );
	  OfflinejetPhi_->Fill(                                                tmpOfflinejetPhi          );
	  OfflinejetEta_->Fill(                                                tmpOfflinejetEta          );
	}
	//      }

      // forward dijet system from VBF (hopefully)
      // -----------------------------------------
      if ( Ic5forwardJetVec.size() >= nforwardjets_ ) {
	double forwardjetEtaProduct = Ic5forwardJetVec[0].eta()*Ic5forwardJetVec[1].eta();
	if ( forwardjetEtaProduct < 0 ) {
	  ++goodIc5forwardJetCounter_;
	  goodIc5forwardJetVec.push_back(Ic5forwardJetVec[0]);
	  goodIc5forwardJetVec.push_back(Ic5forwardJetVec[1]);
	}
      }
      goodIc5forwardJetNumber_->Fill(goodIc5forwardJetVec.size());

      // higgs event selection
      // ---------------------
      if ( goodMuonVec.size() >= nZleptons_ || goodElectronVec.size() >= nZleptons_ ) {
	// leptonic Z
	// ----------
	int lepton1_charge = 0;
	int lepton2_charge = 0;
	int pdgId1 = 0;
	int pdgId2 = 0;
	double leptonicZmassUncertainty = 0.3*ZMass_;
	if ( goodMuonVec.size() >= nZleptons_ ) {
	  pdgId1 = pythiamu_;
	  pdgId2 = pythiamu_;
	  std::vector<GlobalMuon>::const_iterator goodMuonVec_itr1 = goodMuonVec.begin();
	  for ( ; goodMuonVec_itr1 != (goodMuonVec.end()-1); ++goodMuonVec_itr1) {
	    std::vector<GlobalMuon>::const_iterator goodMuonVec_itr2 = (goodMuonVec_itr1+1);
	    for ( ; goodMuonVec_itr2 != goodMuonVec.end(); ++goodMuonVec_itr2) {
	      lepton1_charge = goodMuonVec_itr1->charge();
	      lepton2_charge = goodMuonVec_itr2->charge();
	      if ( lepton1_charge*lepton2_charge < 0 ) {
		if (lepton1_charge < 0 ) pdgId1 = -pythiamu_;
		if (lepton2_charge < 0 ) pdgId2 = -pythiamu_;
		Particle lepton1 = Particle(lepton1_charge,goodMuonVec_itr1->p4(),goodMuonVec_itr1->vertex(),pdgId1,0,true);
		Particle lepton2 = Particle(lepton2_charge,goodMuonVec_itr2->p4(),goodMuonVec_itr2->vertex(),pdgId2,0,true);
		double dileptonmass = DiParticleMass<Particle>(&lepton1,&lepton2);
		if (fabs(dileptonmass-ZMass_) < leptonicZmassUncertainty) {
		  leptonicZmassUncertainty = fabs(dileptonmass-ZMass_);
		  ZleptonsVec.clear();
		  ZleptonsVec.push_back(lepton1);
		  ZleptonsVec.push_back(lepton2);
		}
	      }
	    }
	  }
	}
	if ( goodElectronVec.size() >= nZleptons_ ) {
	  pdgId1 = pythiae_;
	  pdgId2 = pythiae_;
	  std::vector<SimpleElectron>::const_iterator goodElectronVec_itr1 = goodElectronVec.begin();
	  for ( ; goodElectronVec_itr1 != (goodElectronVec.end()-1); ++goodElectronVec_itr1) {
	    std::vector<SimpleElectron>::const_iterator goodElectronVec_itr2 = (goodElectronVec_itr1+1);
	    for ( ; goodElectronVec_itr2 != goodElectronVec.end(); ++goodElectronVec_itr2 ) {
	      lepton1_charge = goodElectronVec_itr1->charge();
	      lepton2_charge = goodElectronVec_itr2->charge();
	      if ( lepton1_charge*lepton2_charge < 0 ) {
		if (lepton1_charge < 0 ) pdgId1 = -pythiae_;
		if (lepton2_charge < 0 ) pdgId2 = -pythiae_;
		Particle lepton1 = Particle(lepton1_charge,goodElectronVec_itr1->p4(),goodElectronVec_itr1->vertex(),pdgId1,0,true);
		Particle lepton2 = Particle(lepton2_charge,goodElectronVec_itr2->p4(),goodElectronVec_itr2->vertex(),pdgId2,0,true);
		double dileptonmass = DiParticleMass<Particle>(&lepton1,&lepton2);
		double leptondeltaEta = lepton1.eta()-lepton2.eta();
		double leptondeltaPhi = DeltaPhi(&lepton1,&lepton2);
		double leptondeltaR   = TMath::Sqrt(pow(leptondeltaEta,2)+
						    pow(leptondeltaPhi,2));
		dileptonMassVSdeltaR_->Fill(dileptonmass,leptondeltaR);
		if (fabs(dileptonmass-ZMass_) < leptonicZmassUncertainty) {
		  leptonicZmassUncertainty = fabs(dileptonmass-ZMass_);
		  ZleptonsVec.clear();
		  ZleptonsVec.push_back(lepton1);
		  ZleptonsVec.push_back(lepton2);
		}
	      }
	    }
	  }
	}
	if ( ZleptonsVec.size() == nZleptons_ ) {
	  leptonicZrec = new Particle(ZCharge_,ZleptonsVec[0].p4()+
				              (ZleptonsVec[1].p4()),ZleptonsVec[0].vertex(),pythiaZ_,0,true);
	  DiLeptonMass = leptonicZrec->mass();
	  double leptonDeltaEta = ZleptonsVec[0].eta()-ZleptonsVec[1].eta();
	  double leptonDeltaPhi = DeltaPhi(&ZleptonsVec[0],&ZleptonsVec[1]);
	  leptonDeltaR          = TMath::Sqrt(pow(leptonDeltaEta,2)+
					      pow(leptonDeltaPhi,2));
	  recZllMass_->Fill(DiLeptonMass);
	  recZllDeltaPhi_->Fill(leptonDeltaPhi);
	  recZllDeltaEta_->Fill(leptonDeltaEta);
	  recZllDeltaR_->Fill(leptonDeltaR);
	}
      }
    }
    if ( goodIc5centralJetVec.size() >= nZjets_ ) {
      // hadronic Z
      // ----------
      double hadronicZmassUncertainty = 0.3*ZMass_;
      std::vector<OfflineJet>::const_iterator goodIc5centralJetVec_itr1 = goodIc5centralJetVec.begin();
      for ( ; goodIc5centralJetVec_itr1 != (goodIc5centralJetVec.end()-1); ++goodIc5centralJetVec_itr1 ) {
	double jet1EtScale = -99.;
	if (goodIc5centralJetVec_itr1->uncorrEt() != 0.) jet1EtScale = goodIc5centralJetVec_itr1->et()/goodIc5centralJetVec_itr1->uncorrEt();
	std::vector<OfflineJet>::const_iterator goodIc5centralJetVec_itr2 = (goodIc5centralJetVec_itr1+1);
	for ( ; goodIc5centralJetVec_itr2 != goodIc5centralJetVec.end(); ++goodIc5centralJetVec_itr2 ) {
	  double jet2EtScale = -99.;
	  if (goodIc5centralJetVec_itr2->uncorrEt() != 0. ) jet2EtScale = goodIc5centralJetVec_itr2->et()/goodIc5centralJetVec_itr2->uncorrEt();
	  Particle bjet1 = Particle( bquarkCharge_, goodIc5centralJetVec_itr1->p4(),goodIc5centralJetVec_itr1->vertex(), pythiab_,0,false);
	  Particle bjet2 = Particle(-bquarkCharge_, goodIc5centralJetVec_itr2->p4(),goodIc5centralJetVec_itr2->vertex(),-pythiab_,0,false);
	  double dijetmass = DiParticleMass<Particle>(&bjet1,&bjet2);
	  Particle scaledbjet1 = Particle( bquarkCharge_, goodIc5centralJetVec_itr1->p4()*jet1EtScale,goodIc5centralJetVec_itr1->vertex(), pythiab_,0,false);
	  Particle scaledbjet2 = Particle(-bquarkCharge_, goodIc5centralJetVec_itr2->p4()*jet2EtScale,goodIc5centralJetVec_itr2->vertex(),-pythiab_,0,false);
	  double scaleddijetmass = DiParticleMass<Particle>(&scaledbjet1,&scaledbjet2);
	  double jetdeltaPhi = DeltaPhi(&scaledbjet1,&scaledbjet2);
	  double jetdeltaEta = scaledbjet1.eta()-scaledbjet2.eta();
	  double jetdeltaR   = TMath::Sqrt(pow(jetdeltaEta,2)+
					   pow(jetdeltaPhi,2));
	  dijetMassVSdeltaR_->Fill(scaleddijetmass,jetdeltaR);
	  dibjetMassVSdeltaR_->Fill(scaleddijetmass,jetdeltaR);
	  if (fabs(scaleddijetmass-ZMass_) <= hadronicZmassUncertainty) {
	    hadronicZmassUncertainty = fabs(scaleddijetmass-ZMass_);
	    double dijetMassComparison = (scaleddijetmass-dijetmass)/scaleddijetmass;
	    scaleddijetmass_->Fill(scaleddijetmass);
	    dijetMassComparison_->Fill(dijetMassComparison);
	    dijetMassVSscaleddibjetMass_->Fill(dijetmass,scaleddijetmass);
	    ZjetsVec.clear();
	    ZjetsVec.push_back(scaledbjet1);
	    ZjetsVec.push_back(scaledbjet2);
	  }
	}
      }
      if ( ZjetsVec.size() == nZjets_ ) {
	hadronicZrec = new Particle(ZCharge_,ZjetsVec[0].p4()+
				            (ZjetsVec[1].p4()),ZjetsVec[0].vertex(),pythiaZ_,0,true);
	centralDiJetMass = hadronicZrec->mass();
	double centralJetDeltaEta = ZjetsVec[0].eta()-ZjetsVec[1].eta();
	double centralJetDeltaPhi = DeltaPhi(&ZjetsVec[0],&ZjetsVec[1]);
	centralJetDeltaR          = TMath::Sqrt(pow(centralJetDeltaEta,2)+
						pow(centralJetDeltaPhi,2));
	recZjjMass_->Fill(centralDiJetMass);
	recZbbMass_->Fill(centralDiJetMass);
	recZjjDeltaR_->Fill(centralJetDeltaR);
	recZbbDeltaR_->Fill(centralJetDeltaR);
	recZjjDeltaPhi_->Fill(centralJetDeltaPhi);
	recZbbDeltaPhi_->Fill(centralJetDeltaPhi);
	recZjjDeltaEta_->Fill(centralJetDeltaEta);
	recZbbDeltaEta_->Fill(centralJetDeltaEta);
      }
    }
  
    double VBFdiJetMass = 100.;
    double forwardDiJetMassMinimum = -1.;
    if ( goodIc5forwardJetVec.size() >= nforwardjets_ ) {
      // VBF dijet system
      // ----------------
      std::vector<OfflineJet>::const_iterator goodIc5forwardJetVec_itr1 = goodIc5forwardJetVec.begin();
      for ( ; goodIc5forwardJetVec_itr1 != (goodIc5forwardJetVec.end()-1); ++goodIc5forwardJetVec_itr1 ) {
	Particle VBFjet1 = Particle( bquarkCharge_, goodIc5forwardJetVec_itr1->p4(),goodIc5forwardJetVec_itr1->vertex(), pythiab_,0,false);
	std::vector<OfflineJet>::const_iterator goodIc5forwardJetVec_itr2 = goodIc5forwardJetVec_itr1+1;
	for ( ; goodIc5forwardJetVec_itr2 != goodIc5forwardJetVec.end(); ++goodIc5forwardJetVec_itr2 ) {
	  Particle VBFjet2 = Particle( bquarkCharge_, goodIc5forwardJetVec_itr2->p4(),goodIc5forwardJetVec_itr2->vertex(), pythiab_,0,false);
	  double dijetmass = DiParticleMass<Particle>(&VBFjet1,&VBFjet2);
	  if ( dijetmass > forwardDiJetMassMinimum ) {
	    forwardDiJetMassMinimum = dijetmass;
	    VBFsystemVec.clear();
	    VBFsystemVec.push_back(VBFjet1);
	    VBFsystemVec.push_back(VBFjet2);
	  }
	}
      }
      if ( VBFsystemVec.size() == nforwardjets_ ) {
	systemVBFrec = new Particle(0,VBFsystemVec[0].p4()+
				     (VBFsystemVec[1].p4()),VBFsystemVec[0].vertex(),0,0,true);
	VBFdiJetMass = systemVBFrec->mass();
	double forwardDiJetDeltaPhi = VBFsystemVec[0].eta()-VBFsystemVec[1].eta();
	double forwardDiJetDeltaEta = DeltaPhi(VBFsystemVec[0].phi(),VBFsystemVec[1].phi());
	double forwardDiJetDeltaR = TMath::Sqrt(pow(forwardDiJetDeltaEta,2)+
						pow(forwardDiJetDeltaPhi,2));
	double forwardJet1Et  = VBFsystemVec[0].et();
	double forwardJet2Et  = VBFsystemVec[1].et();
	double forwardJet1Eta = VBFsystemVec[0].eta();
	double forwardJet2Eta = VBFsystemVec[1].eta();
	  
	forwardJetEt_->Fill(forwardJet1Et);
	forwardJetEt_->Fill(forwardJet2Et);
	forwardJetEta_->Fill(forwardJet1Eta);
	forwardJetEta_->Fill(forwardJet2Eta);
	forwardJet1Eta_->Fill(forwardJet1Eta);
	forwardJet2Eta_->Fill(forwardJet2Eta);
	forwardDiJetDeltaPhi_->Fill(forwardDiJetDeltaPhi);
	forwardDiJetDeltaEta_->Fill(forwardDiJetDeltaEta);
	forwardDiJetDeltaR_->Fill(forwardDiJetDeltaR);
	forwardDiJetMass_->Fill(VBFdiJetMass);  
      }
    }
      
    if ( ZleptonsVec.size() == nZleptons_ && ZjetsVec.size() == nZjets_ 
	 && leptonicZrec && hadronicZrec ) {
      // reconstructed Higgs
      // -------------------
      double recHZZdeltaPhi = DeltaPhi(leptonicZrec,hadronicZrec);
      double recHZZdeltaEta = leptonicZrec->eta()-hadronicZrec->eta();
      double recHZZdeltaR   = -99.;
      recH = new Particle(HCharge_,leptonicZrec->p4()+
		                  (hadronicZrec->p4()),leptonicZrec->vertex(),pythiaH_,true);
      lljjInvMass = recH->mass();
      recHZZdeltaR = TMath::Sqrt(pow(recHZZdeltaEta,2)+
				 pow(recHZZdeltaPhi,2));
      recHZZdeltaPhi_->Fill(recHZZdeltaPhi);
      recHZZdeltaEta_->Fill(recHZZdeltaEta);
      recHZZdeltaR_->Fill(recHZZdeltaR);
    
      if ( VBFsystemVec.size() == nforwardjets_ && systemVBFrec ) {
	lepton1ShiftedEta = fabs(ZleptonsVec[0].eta()-1/(2*(VBFsystemVec[0].eta()+VBFsystemVec[1].eta())));	
	lepton2ShiftedEta = fabs(ZleptonsVec[1].eta()-1/(2*(VBFsystemVec[0].eta()+VBFsystemVec[1].eta())));	
	recZll1ShiftedEta_->Fill(lepton1ShiftedEta);
	recZll2ShiftedEta_->Fill(lepton2ShiftedEta);
	recZllShiftedEta_->Fill(lepton1ShiftedEta);
	recZllShiftedEta_->Fill(lepton2ShiftedEta);

	centralJet1ShiftedEta = fabs(ZjetsVec[0].eta()-1/(2*(VBFsystemVec[0].eta()+VBFsystemVec[1].eta())));	
	centralJet2ShiftedEta = fabs(ZjetsVec[1].eta()-1/(2*(VBFsystemVec[0].eta()+VBFsystemVec[1].eta())));	
	recZjj1ShiftedEta_->Fill(centralJet1ShiftedEta);
	recZjj2ShiftedEta_->Fill(centralJet2ShiftedEta);
	recZbb1ShiftedEta_->Fill(centralJet1ShiftedEta);	  
	recZbb2ShiftedEta_->Fill(centralJet2ShiftedEta);	  
	recZjjShiftedEta_->Fill(centralJet1ShiftedEta);
	recZjjShiftedEta_->Fill(centralJet2ShiftedEta);
	recZbbShiftedEta_->Fill(centralJet1ShiftedEta);	  
	recZbbShiftedEta_->Fill(centralJet2ShiftedEta);	  
	
	Particle balancing = Particle(0,leptonicZrec->p4()+
				      (hadronicZrec->p4())+
				      (systemVBFrec->p4()),systemVBFrec->vertex(),0,true);
	pTbalance = balancing.pt();
	double pTbalancePhi = balancing.phi();
	double pTbalanceMEtDeltaPhi = DeltaPhi(pTbalancePhi,tmpMEtPhi);
	pTbalance_->Fill(pTbalance);
	pTbalancePhi_->Fill(pTbalancePhi);
	pTbalanceVSMEt_->Fill(tmpMEt,pTbalance);
	pTbalanceMEtDeltaPhi_->Fill(pTbalanceMEtDeltaPhi);
      }
    
      // selection
      // ---------
      recHlljjVSforwardDiJetMass_->Fill(VBFdiJetMass,lljjInvMass);
      recHlljjVSpTbalance_->Fill(pTbalance,lljjInvMass);       
      recHlljjVSZdileptonMass_->Fill(DiLeptonMass,lljjInvMass);   
      recHlljjVSleptonDeltaR_->Fill(leptonDeltaR,lljjInvMass);    
      recHlljjVSleptonShiftedEta_->Fill(lepton1ShiftedEta,lljjInvMass);
      recHlljjVSleptonShiftedEta_->Fill(lepton2ShiftedEta,lljjInvMass);
      recHlljjVSZdijetMass_->Fill(centralDiJetMass,lljjInvMass);      
      recHlljjVSZjetDeltaR_->Fill(centralJetDeltaR,lljjInvMass);      
      recHlljjVSZjetShiftedEta_->Fill(centralJet1ShiftedEta,lljjInvMass);  
      recHlljjVSZjetShiftedEta_->Fill(centralJet2ShiftedEta,lljjInvMass);  
      
      if ( ( VBFdiJetMass >= forwardDiJetMassCut_  
	     && pTbalance <= pTbalanceCut_ 
	     )
	   // leptonic Z cuts
	   && ( ( DiLeptonMass         >= ZdileptonMassminCut_    && DiLeptonMass      <= ZdileptonMassmaxCut_ )
		&& ( leptonDeltaR      >= ZdileptonDeltaRminCut_  && leptonDeltaR      <= ZdileptonDeltaRmaxCut_ )
		&& ( lepton1ShiftedEta <= ZdileptonShiftedEtaCut_ && lepton2ShiftedEta <= ZdileptonShiftedEtaCut_ ) 
		)
	   // hadronic Z
	   && ( ( centralDiJetMass         >= ZdijetMassminCut_    && centralDiJetMass      <= ZdijetMassmaxCut_ ) 
		&& ( centralJetDeltaR      >= ZdijetDeltaRminCut_  && centralJetDeltaR      <= ZdijetDeltaRmaxCut_ ) 
		&& ( centralJet1ShiftedEta <= ZdijetShiftedEtaCut_ && centralJet2ShiftedEta <= ZdijetShiftedEtaCut_ ) 
		) 
	   )
	recHlljj_->Fill(lljjInvMass);
      
      
      if ( VBFdiJetMass >= forwardDiJetMassCut_ ) recHlljj_forwardDiJetMass->Fill(lljjInvMass);	     
      if ( pTbalance        <= pTbalanceCut_        ) recHlljj_pTbalance->Fill(lljjInvMass);	     
      // leptonic Z cuts
      if ( DiLeptonMass>=ZdileptonMassminCut_         && DiLeptonMass<=ZdileptonMassmaxCut_           ) recHlljj_ZdileptonMass->Fill(lljjInvMass);	     
      if ( leptonDeltaR>=ZdileptonDeltaRminCut_       && leptonDeltaR<=ZdileptonDeltaRmaxCut_         ) recHlljj_leptonDeltaR->Fill(lljjInvMass);	     
      if ( lepton1ShiftedEta<=ZdileptonShiftedEtaCut_ && lepton2ShiftedEta<=ZdileptonShiftedEtaCut_   ) recHlljj_leptonShiftedEta->Fill(lljjInvMass);
      // hadronic Z
      if ( centralDiJetMass>=ZdijetMassminCut_         && centralDiJetMass<=ZdijetMassmaxCut_        ) recHlljj_ZdijetMass->Fill(lljjInvMass);	     
      if ( centralJetDeltaR>=ZdijetDeltaRminCut_       && centralJetDeltaR<=ZdijetDeltaRmaxCut_      ) recHlljj_ZjetDeltaR->Fill(lljjInvMass);	     
      if ( centralJet1ShiftedEta<=ZdijetShiftedEtaCut_ && centralJet2ShiftedEta<=ZdijetShiftedEtaCut_) recHlljj_ZjetShiftedEta->Fill(lljjInvMass);	     
      
      if ( VBFdiJetMass >= forwardDiJetMassCut_ ) {
	recHlljj_AFTERforwardDiJetMass->Fill(lljjInvMass);	     
	if ( pTbalance <= pTbalanceCut_ ) {
	  recHlljj_AFTERpTbalance->Fill(lljjInvMass);	     
	  // leptonic Z cuts
	  if ( DiLeptonMass>=ZdileptonMassminCut_ && DiLeptonMass<=ZdileptonMassmaxCut_ ) {
	    recHlljj_AFTERZdileptonMass->Fill(lljjInvMass);	     
	    if ( leptonDeltaR>=ZdileptonDeltaRminCut_ && leptonDeltaR<=ZdileptonDeltaRmaxCut_ ) {
	      recHlljj_AFTERleptonDeltaR->Fill(lljjInvMass);	     
	      if ( lepton1ShiftedEta<=ZdileptonShiftedEtaCut_ && lepton2ShiftedEta<=ZdileptonShiftedEtaCut_ ) {
		recHlljj_AFTERleptonShiftedEta->Fill(lljjInvMass);
		// hadronic Z
		if ( centralDiJetMass>=ZdijetMassminCut_ && centralDiJetMass<=ZdijetMassmaxCut_ ) {
		  recHlljj_AFTERZdijetMass->Fill(lljjInvMass);	     
		  if ( centralJetDeltaR>=ZdijetDeltaRminCut_ && centralJetDeltaR<=ZdijetDeltaRmaxCut_ ) {
		    recHlljj_AFTERZjetDeltaR->Fill(lljjInvMass);	     
		    if ( centralJet1ShiftedEta<=ZdijetShiftedEtaCut_ && centralJet2ShiftedEta<=ZdijetShiftedEtaCut_ ) {
		      recHlljj_AFTERZjetShiftedEta->Fill(lljjInvMass);	     
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if ( Z1particle || Z2particle ) {
      double hadronicZmass = -99.;
      double leptonicZmass = -99.;
      if ( Z1particlesId <= pythiat_ ) 
	hadronicZmass = Z1mass;
      else if ( Z1particlesId <= pythianutau_ ) {
	leptonicZmass = Z1mass;
      }
      if ( Z2particlesId <= pythiat_ ) 
	hadronicZmass = Z2mass;
      else if ( Z2particlesId <= pythianutau_ ) {
	leptonicZmass = Z2mass;
      }
      if ( leptonicZrec ) {
	if ( Z1particlesId == pythiae_ && fabs(ZleptonsVec[0].pdgId()) == pythiae_ ) {
	  recZeeMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	  recZllMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	}
	if ( Z2particlesId == pythiae_ && fabs(ZleptonsVec[0].pdgId()) == pythiae_ ) { 
	  recZeeMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	  recZllMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	}
	if ( Z1particlesId == pythiamu_ && fabs(ZleptonsVec[0].pdgId()) == pythiamu_ ) {
	  recZmmMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	  recZllMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	}
	if ( Z2particlesId == pythiamu_ && fabs(ZleptonsVec[0].pdgId()) == pythiamu_ ) {
	  recZmmMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	  recZllMassResolution_->Fill((leptonicZmass - DiLeptonMass)/leptonicZmass);
	}
      }
      if ( hadronicZrec ) {
	recZjjMassResolution_->Fill((hadronicZmass - centralDiJetMass)/hadronicZmass);
	if ( Z1particlesId == pythiab_ )
	  recZbbMassResolution_->Fill((hadronicZmass - centralDiJetMass)/hadronicZmass);
	if ( Z2particlesId == pythiab_ )
	  recZbbMassResolution_->Fill((hadronicZmass - centralDiJetMass)/hadronicZmass);
      }
      if ( Hparticle ) {
	if ( leptonicZrec && hadronicZrec && systemVBFrec ) {
	  if ( ( VBFdiJetMass >= forwardDiJetMassCut_  
		 && pTbalance <= pTbalanceCut_ 
		 )
	       // leptonic Z cuts
	       && ( ( DiLeptonMass         >= ZdileptonMassminCut_    && DiLeptonMass      <= ZdileptonMassmaxCut_ )
		    && ( leptonDeltaR      >= ZdileptonDeltaRminCut_  && leptonDeltaR      <= ZdileptonDeltaRmaxCut_ )
		    && ( lepton1ShiftedEta <= ZdileptonShiftedEtaCut_ && lepton2ShiftedEta <= ZdileptonShiftedEtaCut_ ) 
		    )
	       // hadronic Z
	       && ( ( centralDiJetMass         >= ZdijetMassminCut_    && centralDiJetMass      <= ZdijetMassmaxCut_ ) 
		    && ( centralJetDeltaR      >= ZdijetDeltaRminCut_  && centralJetDeltaR      <= ZdijetDeltaRmaxCut_ ) 
		    && ( centralJet1ShiftedEta <= ZdijetShiftedEtaCut_ && centralJet2ShiftedEta <= ZdijetShiftedEtaCut_ ) 
		    ) 
	       ) {
	    if (((fabs(Hparticle1DaughterId) <= pythiat_ && (fabs(Hparticle2DaughterId) == pythiae_ || fabs(Hparticle2DaughterId) == pythiamu_)) ||
		 (fabs(Hparticle2DaughterId) <= pythiat_ && (fabs(Hparticle1DaughterId) == pythiae_ || fabs(Hparticle1DaughterId) == pythiamu_)) ) 
		) { 
	      recHlljjMassResolution_->Fill((Hmass - lljjInvMass)/Hmass);
	      if (fabs(Hparticle1DaughterId) == pythiab_ ) recHllbbMassResolution_->Fill((Hmass - lljjInvMass)/Hmass);
	      if (fabs(Hparticle2DaughterId) == pythiab_ ) recHllbbMassResolution_->Fill((Hmass - lljjInvMass)/Hmass);
	    }
	  }
	}
      }
    }
    delete leptonicZrec;
    delete hadronicZrec;
    delete systemVBFrec;
    delete recH;
    delete Z1particle;
    delete Z2particle;
    delete Hparticle;

  } else {
    int ZtoTauTauCounter = 0;
    int hadronicZ1Counter = 0;
    int hadronicZ2Counter = 0;
    const Candidate* tmpMotherZ = 0;
    CandidateCollection::const_iterator MCparticle = MCparticles->begin();
    for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
      int tmpParticleId = MCparticle->pdgId();
      if ( fabs(tmpParticleId) > pythiaH_ ) continue;
      int tmpMotherId = 0;
      if( MCparticle->mother() != 0 ) tmpMotherId = MCparticle->mother()->pdgId();

      if ( tmpParticleId == pythiatau_ && fabs(tmpMotherId) == pythiaZ_ && MCparticle->mother()->mother()->pdgId() != tmpMotherId )
	ZtoTauTauCounter++;

      if ( tmpMotherId == pythiaZ_ && tmpParticleId != tmpMotherId ) {
	if ( tmpMotherZ == 0 ) {
	  tmpMotherZ = MCparticle->mother();
	  hadronicZ1Counter++;
	} else {
	  if ( MCparticle->mother() == tmpMotherZ ) {
	    hadronicZ1Counter++;
	  } else {
	    hadronicZ2Counter++;
	  }
	}
      }
    }
    bool hadronicZZ = false;
    if ( hadronicZ1Counter == 2 && hadronicZ2Counter == 2 ) hadronicZZ = true;

    double failedTopologyCounter = double(globalMuonsSize)*1.+double(simpleElectronsSize)*10.+double(offlineJetsSize)*100.;
    failedTopologyNumber_->Fill(failedTopologyCounter);
    //    std::cout << failedTopologyCounter; 
    //    if ( ZtoTauTauCounter >=1 ) std::cout << " while ZtoTauTau: " << ZtoTauTauCounter << std::endl;
    //    if ( hadronicZZ           ) std::cout << " both Z hadronically decay" << std::endl;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HiggsZZAnalyzer::beginJob(const edm::EventSetup&)
{

  vecTypeIndexPair.push_back(std::pair<int,int>(pythiab_, 0));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiae_, 1));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiamu_,2));

  std::vector< std::pair<int,int> >::const_iterator vecTypeIndexPair_itr = vecTypeIndexPair.begin();
  for ( ; vecTypeIndexPair_itr != vecTypeIndexPair.end(); ++vecTypeIndexPair_itr ) 
    vecZtoParticleseventcounter_.push_back(0);
  vecZtoParticleseventcounter_.push_back(0);
  vecZtoParticleseventcounter_.push_back(0);
  
  // File for output histograms
  // --------------------------
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","HiggsZZOutput");
  // The file must be opened first, 
  // so that becomes the default position for all the histograms
  OutputFile->cd();
  // White background for the canvases
  gROOT->SetStyle("Plain");
  
  OfflinejetNumber_             = new TH1D("OffilinejetNumber",           "number of jets",                                   nbin_,   0.,250.);
  OfflinejetUncorrEt_           = new TH1D("OfflinejetUncorrEt",          "offline jet uncor E_{T}",                          nbin_,   0.,500.);
  OfflinejetEt_                 = new TH1D("OfflinejetEt",                "offline jet E_{T}",                                nbin_,   0.,500.);
  OfflinejetEtScale_            = new TH1D("OfflinejetEtScale",           "offline jet E_{T} scale",                          nbin_,   0.,  5.);
  OfflinejetEtScaleVSUncorrEt_  = new TH2D("OfflinejetEtScaleVSUncorrEt", "offline jet E_{T} scale VS uncorr E_{T}",          nbin_,   0.,500.,nbin_,0.,  5.);
  OfflinejetEtScaleVSEt_        = new TH2D("OfflinejetEtScaleVSEt",       "offline jet E_{T} scale VS E_{T}",                 nbin_,   0.,500.,nbin_,0.,  5.);
  OfflinejetEtScaleVSEta_       = new TH2D("OfflinejetEtScaleVSEta",      "offline jet E_{T} scale VS #eta",                  nbin_, -10., 10.,nbin_,0.,  5.);
  OfflinejetEtVSUncorrEt_       = new TH2D("OfflinejetEtVSUncorrEt",      "offline jet E_{T} VS uncorr E_{T}",                nbin_,   0.,500.,nbin_,0.,500.);
  OfflinejetEtVSEt_             = new TH2D("OfflinejetEtVSEt",            "offline jet E_{T} VS E_{T}",                       nbin_,   0.,500.,nbin_,0.,500.);
  OfflinejetEtVSEta_            = new TH2D("OfflinejetEtVSEta",           "offline jet E_{T} VS #eta",                        nbin_, -10., 10.,nbin_,0.,500.);
  OfflinejetUncorrEtVSUncorrEt_ = new TH2D("OfflinejetUncorrEtVSUncorrEt","offline jet uncorr E_{T} VS uncorr E_{T}",         nbin_,   0.,500.,nbin_,0.,500.);
  OfflinejetUncorrEtVSEt_       = new TH2D("OfflinejetUncorrEtVSEt",      "offline jet uncorr E_{T} VS E_{T}",                nbin_,   0.,500.,nbin_,0.,500.);
  OfflinejetUncorrEtVSEta_      = new TH2D("OfflinejetUncorrEtVSEta",     "offline jet uncorr E_{T} VS #eta",                 nbin_, -10., 10.,nbin_,0.,500.);
  OfflinejetEMenergyfrac_       = new TH1D("OfflinejetEMenergyfrac",      "offline jet EM energy fraction",                   nbin_,   0.,         2.);
  OfflinejetHighEffDiscr_       = new TH1D("OfflinejetHighEffDiscr",      "offline jet High Eff discriminator (trk counting)",nbin_,   0.,        50.);
  OfflinejetHighPurDiscr_       = new TH1D("OfflinejetHighPurDiscr",      "offline jet High Pur discriminator (trk counting)",nbin_,   0.,        50.);
  OfflinejetMass_               = new TH1D("OfflinejetMass",              "offline jet mass",                                 nbin_,   0.,       100.);
  OfflinejetTKnumS2_            = new TH1D("OfflinejetTKnumS2",           "offline jet number of tracks S2",                  nbin_,   0.,        50.);
  OfflinejetTKsumPtS2_          = new TH1D("OfflinejetTKsumPtS2",         "offline jet #Sigma p_{T}^{tracks S2}",             nbin_,   0.,       100.);
  OfflinejetTKtagmassS2_        = new TH1D("OfflinejetTKtagmassS2",       "offline jet tagmass from tracks S2",               nbin_,   0.,        10.);
  OfflinejetPhi_                = new TH1D("OfflinejetPhi",               "offline jet #Phi",                                 nbin_, -TMath::Pi(),TMath::Pi() );
  OfflinejetEta_                = new TH1D("OfflinejetEta",               "offline jet #eta",                                 nbin_, -10.,         10.);
  OfflinejetLeptonDeltaR_ = new TH1D("OfflinejetLeptonDeltaR","#DeltaR between the jet and lepton",        nbin_,0.,10.);
  centralJetLeptonDeltaR_ = new TH1D("centralJetLeptonDeltaR","#DeltaR between the central jet and lepton",nbin_,0.,10.);

  forwardJetEt_           = new TH1D("forwardjetEt",          "offline forward jets E_{T}",           nbin_,  0., 5000.);
  forwardJetLeptonDeltaR_ = new TH1D("forwardJetLeptonDeltaR","#DeltaR beteen forward jet and lepton",nbin_,  0.,   10.);
  forwardJetEta_          = new TH1D("forwardJetEta",         "offline forward jets #eta",            nbin_,-10.,   10.);
  forwardJet1Eta_         = new TH1D("forwardJet1Eta",        "offline forward jet1 #eta",            nbin_,-10.,   10.);
  forwardJet2Eta_         = new TH1D("forwardJet2Eta",        "offline forward jet2 #eta",            nbin_,-10.,   10.);
  forwardDiJetMass_       = new TH1D("forwardDiJetMass",      "di-jet mass for forward jets",         nbin_,  0.,14000.);
  forwardDiJetDeltaPhi_   = new TH1D("forwardDiJetDeltaPhi",  "#Delta#Phi between the 2 forward jets",nbin_,  0.,TMath::Pi());
  forwardDiJetDeltaEta_   = new TH1D("forwardDiJetDeltaEta",  "#Delta#eta between the 2 forward jets",nbin_,-10.,   10.);
  forwardDiJetDeltaR_     = new TH1D("forwardDiJetDeltaR",    "#DeltaR between the 2 forward jets",   nbin_,  0.,   10.);
  
  pTbalance_            = new TH1D("pTbalance",   "balance p_{T}",nbin_,0.,1000.);
  pTbalancePhi_         = new TH1D("pTbalancePhi","balance #Phi", nbin_,0.,TMath::Pi());
  pTbalanceVSMEt_       = new TH2D("pTbalanceVSMEt","balance p_{T} VS missing E_{t}",nbin_,0.,1000.,nbin_,0.,1000.);
  pTbalanceMEtDeltaPhi_ = new TH1D("pTbalanceMEtDeltaPhi","#Delta#Phi between balance p_{T} and missing E_{T}",nbin_,0.,TMath::Pi());

  MuonNumber_     = new TH1D("MuonNumber",    "number of muons",        nbin_,   0.,         10.);
  MuonEta_        = new TH1D("MuonEta",       "muon #eta",              nbin_, -10.,         10.);
  MuonPhi_        = new TH1D("MuonPhi",       "muon #Phi",              nbin_, -TMath::Pi(),TMath::Pi() );
  MuonPt_         = new TH1D("MuonPt",        "muon p_{T}",             nbin_,   0.,       1000.);
  MuonPx_         = new TH1D("MuonPx",        "muon px",                nbin_,-250.,        250.);
  MuonPy_         = new TH1D("MuonPy",        "muon py",                nbin_,-250.,        250.);
  MuonPz_         = new TH1D("MuonPz",        "muon pz",                nbin_,-250.,        250.);
  MuonEMcalo_     = new TH1D("MuonEMcalo",    "muon EM calo",           nbin_,   0.,        100.);
  MuonHADcalo_    = new TH1D("MuonHADcalo",   "muon HAD calo",          nbin_,   0.,        100.);
  MuonHOcalo_     = new TH1D("MuonHOcalo",    "muon HO calo",           nbin_,   0.,        100.);
  MuonIso03SumPt_ = new TH1D("MuonIso03SumPt","muon Iso03 #Sigma p_{T}",nbin_,   0.,        100.);
  MuonNormChi2_   = new TH1D("MuonNormChi2",  "muon |#Chi^{2}|",        nbin_,   0.,         50.);
  MuonHitsNum_    = new TH1D("MuonHitsNum",   "muon number of hits",    nbin_,   0.,         50.);
  MuonImpactParamXY_      = new TH1D("MuonImpactParamXY",     "muon impact parameter",               nbin_,-1.,1.);
  MuonErrorImpactParamXY_ = new TH1D("MuonErrorImpactParamXY","muon error on impact parameter",      nbin_,0.,0.5);
  MuonImpactParamSig_     = new TH1D("MuonImpactParamSig",    "muon significance on impact parametr",nbin_,0.,10.);


  ElectronNumber_    = new TH1D("ElectronNumber",   "number of electrons", nbin_,   0.,         10.);
  ElectronEt_        = new TH1D("ElectronEt",       "electron E_{T}",      nbin_,   0.,       1000.);
  ElectronEta_       = new TH1D("ElectronEta",      "electron #eta",       nbin_, -10.,         10.);
  ElectronPhi_       = new TH1D("ElectronPhi",      "electron #Phi",       nbin_, -TMath::Pi(),TMath::Pi() );
  ElectronPt_        = new TH1D("ElectronPt",       "electron p_{T}",      nbin_,   0.,       1000.);
  ElectronPx_        = new TH1D("ElectronPx",       "electron px",         nbin_,-250.,        250.);
  ElectronPy_        = new TH1D("ElectronPy",       "electron py",         nbin_,-250.,        250.);
  ElectronPz_        = new TH1D("ElectronPz",       "electron pz",         nbin_,-250.,        250.);
  ElectronHadOverEm_ = new TH1D("ElectronHadOverEm","electron HAD/EM",     nbin_,   0.,          1.);
  ElectronIsoVal_    = new TH1D("ElectronIsoVal",   "electron Iso val",    nbin_,   0.,        100.);
  ElectronImpactParamXY_      = new TH1D("ElectronImpactParamXY",     "electron impact parameter",               nbin_,-1.,1.);
  ElectronErrorImpactParamXY_ = new TH1D("ElectronErrorImpactParamXY","electron error on impact parameter",      nbin_,0.,0.5);
  ElectronImpactParamSig_     = new TH1D("ElectronImpactParamSig",    "electron significance on impact parametr",nbin_,0.,10.);

  eventsNumber_             = new TH1D("eventsNumber",            "total number of events",                   1,0.,1.);
  signalTopologySizeNumber_ = new TH1D("signalTopologySizeNumber","number of events w/ 2 leptons and 4 jets", 1,0.,1.);
  failedTopologyNumber_     = new TH1D("failedTopologyNumber",    "number of events w/o 2 leptons and 4 jets",5000,0.,5000.);
  HtoZZEventsNumber_        = new TH1D("HtoZZeventsNumber",       "number of events w/ H->ZZ",                1,0.,1.);

  muonNumber_              = new TH1D("muonNumber",             "number of muons",                                                      nbin_,0.,100.);
  goodMuonNumber_          = new TH1D("goodMuonNumber",         "number of muons w/ p_{T}>20 GeV/c & \\frac{d_{0}}{#sigma_{d_{0}}}<3.5",    nbin_,0., 10.);
  electronNumber_          = new TH1D("electronNumber",         "number of electrons",                                                  nbin_,0.,100.);
  goodElectronNumber_      = new TH1D("goodElectronNumber",     "number of electrons w/ p_{T}>20 GeV/c & \\frac{d_{0}}{#sigma_{d_{0}}}<3.5",nbin_,0., 10.);
  goodIc5JetNumber_        = new TH1D("goodIc5JetNumber",       "number of jets w/ E_{T}>20. GeV/c & #DeltaR_{lj} > ",nbin_,0., 20.);
  goodIc5centralJetNumber_ = new TH1D("goodIc5centralJetNumber","number of good jets w/ |#eta|<2.0 & EMF<0.8",           nbin_,0., 10.);
  Ic5forwardJetNumber_     = new TH1D("Ic5forwardJetNumber",    "number of good jets w/ |#eta|>2.0",                   nbin_,0.,100.);
  goodIc5forwardJetNumber_ = new TH1D("goodIc5forwardJetNumber","number of forward jets/ #eta_{j1}*#eta_{j2}<0",       nbin_,0., 10.);

  vecZtoParticlesEventsNumber_.push_back( TH1D("ZtobbbarEventsNumber",    "number of events w/ Z->b\\bar{b}",      1,0.,1.));
  vecZtoParticlesEventsNumber_.push_back( TH1D("ZtoElectronsEventsNumber","number of events w/ Z->e^{+}e^{-}",     1,0.,1.));
  vecZtoParticlesEventsNumber_.push_back( TH1D("ZtoMuonsEventsNumber",    "number of events w/ Z->#mu^{+}#mu^{-}", 1,0.,1.));
  vecZtoParticlesEventsNumber_.push_back( TH1D("ZtoQuarksEventsNumber",   "number of events w/ Z->q\\bar{q}",      1,0.,1.));
  vecZtoParticlesEventsNumber_.push_back( TH1D("ZtoLeptonsEventsNumber",  "number of events w/ Z->l^{+}l^{-}",     1,0.,1.));

  MEt_             = new TH1D("MEt",            "missing E_{T}",                             nbin_,0.,1000.);
  MEtSumEt_        = new TH1D("MEtSumEt",       "#SigmaE_{T} from missing E_{T} information",nbin_,0.,1000.);
  MEtSignificance_ = new TH1D("MEtSignificance","missing E_{T} significance",                nbin_,0.,  50.);
  MEtPhi_          = new TH1D("MEtPhi",         "missing E_{T} #phi",                        nbin_,0.,TMath::Pi());

  dileptonMassVSdeltaR_ = new TH2D("dileptonmassVSdeltaR","di-lepton mass VS #DeltaR",nbin_,0.,200.,nbin_,0.,10.);
  dijetMassVSdeltaR_    = new TH2D("dijetmassVSdeltaR",   "di-jet mass VS #DeltaR",   nbin_,0.,200.,nbin_,0.,10.);
  dibjetMassVSdeltaR_   = new TH2D("dibjetmassVSdeltaR",  "di-bjet mass VS #DeltaR",  nbin_,0.,200.,nbin_,0.,10.);
  scaleddijetmass_             = new TH1D("scaleddijetmass","scaled di-jet mass",                           nbin_, 0.,200.);
  dijetMassComparison_         = new TH1D("dijetMassComparison","\\frac{m_{jj}^{scaled}-m_{jj}}{m_{jj}^{scaled}}",nbin_,-5.,  5.);
  dijetMassVSscaleddibjetMass_ = new TH2D("dijetMassVSscaleddibjetMass","scaled di-jet mass VS di-jet mass",nbin_, 0.,200.,nbin_,0.,200.);

  recZllMass_           = new TH1D("recZllMass",          "reconstructed dilepton invariant mass",             nbin_,  0.,200.);
  recZjjMass_           = new TH1D("recZjjMass",          "reconstructed dijet invariant mass",                nbin_,  0.,200.);
  recZbbMass_           = new TH1D("recZbbMass",          "reconstructed dib-jet invariant mass resolution",   nbin_,  0.,200.);
  recZllMassResolution_ = new TH1D("recZllMassResolution","reconstructed dilepton invariant mass resolution",  nbin_, -3.,  3.);
  recZeeMassResolution_ = new TH1D("recZeeMassResolution","reconstructed dielectron invariant mass resolution",nbin_, -3.,  3.);
  recZmmMassResolution_ = new TH1D("recZmmMassResolution","reconstructed dimuon invariant mass resolution",    nbin_, -3.,  3.);
  recZjjMassResolution_ = new TH1D("recZjjMassResolution","reconstructed dijet invariant mass resolution",     nbin_, -3.,  3.);
  recZbbMassResolution_ = new TH1D("recZbbMassResolution","reconstructed dib-jet invariant mass resolution",   nbin_, -3.,  3.);

  recZllDeltaEta_       = new TH1D("recZllDeltaEta",      "#Delta#eta between the reconstructed leptons from Z",nbin_,-10., 10.);
  recZjjDeltaEta_       = new TH1D("recZjjDeltaEta",      "#Delta#eta between the reconstructed jets from Z",   nbin_,-10., 10.);
  recZbbDeltaEta_       = new TH1D("recZbbDeltaEta",      "#Delta#eta between the reconstructed b-jets from Z", nbin_,-10., 10.);
  recZllDeltaPhi_       = new TH1D("recZllDeltaPhi",      "#Delta#Phi between the reconstructed leptons from Z",nbin_,  0., TMath::Pi());
  recZjjDeltaPhi_       = new TH1D("recZjjDeltaPhi",      "#Delta#Phi between the reconstructed jets from Z",   nbin_,  0., TMath::Pi());
  recZbbDeltaPhi_       = new TH1D("recZbbDeltaPhi",      "#Delta#Phi between the reconstructed b-jets from Z", nbin_,  0., TMath::Pi());
  recZllDeltaR_         = new TH1D("recZllDeltaR",        "#DeltaR between the reconstruced leptons from Z",    nbin_,  0., 10.);
  recZjjDeltaR_         = new TH1D("recZjjDeltaR",        "#DeltaR between the reconstruced jets from Z",       nbin_,  0., 10.);
  recZbbDeltaR_         = new TH1D("recZbbDeltaR",        "#DeltaR between the reconstruced b-jets from Z",     nbin_,  0., 10.);
  recZll1ShiftedEta_    = new TH1D("recZll1ShiftedEta",   "shifted #eta of the 1^{st} lepton from Z",           nbin_,  0., 10.);
  recZjj1ShiftedEta_    = new TH1D("recZjj1ShiftedEta",   "shifted #eta of the 1^{nd} jet from Z",              nbin_,  0., 10.);
  recZbb1ShiftedEta_    = new TH1D("recZbb1ShiftedEta",   "shifted #eta of the 1^{nd} b-jet from Z",            nbin_,  0., 10.);
  recZll2ShiftedEta_    = new TH1D("recZll2ShiftedEta",   "shifted #eta of the 2^{nd} lepton from Z",           nbin_,  0., 10.);
  recZjj2ShiftedEta_    = new TH1D("recZjj2ShiftedEta",   "shifted #eta of the 2^{nd} jet from Z",              nbin_,  0., 10.);
  recZbb2ShiftedEta_    = new TH1D("recZbb2ShiftedEta",   "shifted #eta of the 2^{nd} b-jet from Z",            nbin_,  0., 10.);
  recZllShiftedEta_     = new TH1D("recZllShiftedEta",    "shifted #eta of leptons from Z",                     nbin_,  0., 10.);
  recZjjShiftedEta_     = new TH1D("recZjjShiftedEta",    "shifted #eta of jets from Z",                        nbin_,  0., 10.);
  recZbbShiftedEta_     = new TH1D("recZbbShiftedEta",    "shifted #eta of b-jets from Z",                      nbin_,  0., 10.);

  recZjjLeptonDeltaR_ = new TH1D("recZjjLeptonDeltaR","#DeltaR between the jet from Z and the lepton",  nbin_,0.,10.);
  recZbbLeptonDeltaR_ = new TH1D("recZbbLeptonDeltaR","#DeltaR between the b-jet from Z and the lepton",nbin_,0.,10.);

  recHZZdeltaPhi_           = new TH1D("recHZZdeltaPhi",           "#Delta#Phi between the 2 reconstructed Z",                                        nbin_,  0.,TMath::Pi());
  recHZZdeltaEta_           = new TH1D("recHZZdeltaEta",           "#Delta#eta between the 2 reconstructed Z",                                        nbin_,-10.,  10.);
  recHZZdeltaR_             = new TH1D("recHZZdeltaR",             "#DeltaR between the 2 reconstructed Z",                                           nbin_,  0.,  10.);
  recHlljj_                 = new TH1D("recHlljj",                 "reconstructed 4-body invarint mass (lljj)",                                       nbin_,  0.,1000.);
  recHlljj_forwardDiJetMass = new TH1D("recHlljj_forwardDiJetMass","reconstructed 4-body invariant mass after forward/tagging m_{jj} cut",nbin_,  0.,1000.);
  recHlljj_pTbalance        = new TH1D("recHlljj_pTbalance",       "reconstructed 4-body invariant mass after p_{T}^{balance} cut",       nbin_,  0.,1000.);
  recHlljj_ZdileptonMass    = new TH1D("recHlljj_ZdileptonMass",   "reconstructed 4-body invariant mass after m_{ll} cut",                nbin_,  0.,1000.);
  recHlljj_leptonDeltaR     = new TH1D("recHlljj_leptonDeltaR",    "reconstructed 4-body invariant mass after #DeltaR_{ll} cut",          nbin_,  0.,1000.);
  recHlljj_leptonShiftedEta = new TH1D("recHlljj_leptonShiftedEta","reconstructed 4-body invariant mass after Shifted |#eta_{l}| cut",    nbin_,  0.,1000.);
  recHlljj_ZdijetMass       = new TH1D("recHlljj_ZdijetMass",      "reconstructed 4-body invariant mass after m_{jj} cut",                nbin_,  0.,1000.);
  recHlljj_ZjetDeltaR       = new TH1D("recHlljj_ZjetDeltaR",      "reconstructed 4-body invariant mass after #DeltaR_{jj} cut",          nbin_,  0.,1000.);
  recHlljj_ZjetShiftedEta   = new TH1D("recHlljj_ZjetShiftedEta",  "reconstructed 4-body invariant mass after Shifted |#eta_{j}| cut",    nbin_,  0.,1000.);
  recHlljjVSforwardDiJetMass_ = new TH2D("recHlljjVSforwardDiJetMass","reconstructed 4-body invariant mass VS forward/tagging m_{jj}",nbin_,0.,1000.,nbin_,0.,1000.);
  recHlljjVSpTbalance_        = new TH2D("recHlljjVSpTbalance",       "reconstructed 4-body invariant mass VS p_{T}^{balance}",       nbin_,0.,1000.,nbin_,0.,1000.);
  recHlljjVSZdileptonMass_    = new TH2D("recHlljjVSZdileptonMass",   "reconstructed 4-body invariant mass VS m_{ll}",                nbin_,0., 200.,nbin_,0.,1000.);
  recHlljjVSleptonDeltaR_     = new TH2D("recHlljjVSleptonDeltaR",    "reconstructed 4-body invariant mass VS #DeltaR_{ll}",          nbin_,0.,  10.,nbin_,0.,1000.);
  recHlljjVSleptonShiftedEta_ = new TH2D("recHlljjVSleptonShiftedEta","reconstructed 4-body invariant mass VS Shifted |#eta_{l}|",    nbin_,0.,  10.,nbin_,0.,1000.);
  recHlljjVSZdijetMass_       = new TH2D("recHlljjVSZdijetMass",      "reconstructed 4-body invariant mass VS m_{jj}",                nbin_,0., 200.,nbin_,0.,1000.);
  recHlljjVSZjetDeltaR_       = new TH2D("recHlljjVSZjetDeltaR",      "reconstructed 4-body invariant mass VS #DeltaR_{jj}",          nbin_,0.,  10.,nbin_,0.,1000.);
  recHlljjVSZjetShiftedEta_   = new TH2D("recHlljjVSZjetShiftedEta",  "reconstructed 4-body invariant mass VS Shifted |#eta_{j}|",    nbin_,0.,  10.,nbin_,0.,1000.);
  recHlljj_AFTERforwardDiJetMass = new TH1D("recHlljj_AFTERforwardDiJetMass","reconstructed 4-body invariant mass after forward/tagging m_{jj} cut",nbin_,0.,1000.);
  recHlljj_AFTERpTbalance        = new TH1D("recHlljj_AFTERpTbalance",       "reconstructed 4-body invariant mass after p_{T}^{balance} cut",       nbin_,0.,1000.);
  recHlljj_AFTERZdileptonMass    = new TH1D("recHlljj_AFTERZdileptonMass",   "reconstructed 4-body invariant mass after m_{ll} cut",                nbin_,0.,1000.);
  recHlljj_AFTERleptonDeltaR     = new TH1D("recHlljj_AFTERleptonDeltaR",    "reconstructed 4-body invariant mass after #DeltaR_{ll} cut",          nbin_,0.,1000.);
  recHlljj_AFTERleptonShiftedEta = new TH1D("recHlljj_AFTERleptonShiftedEta","reconstructed 4-body invariant mass after Shifted |#eta_{l}| cut",    nbin_,0.,1000.);
  recHlljj_AFTERZdijetMass       = new TH1D("recHlljj_AFTERZdijetMass",      "reconstructed 4-body invariant mass after m_{jj} cut",                nbin_,0.,1000.);
  recHlljj_AFTERZjetDeltaR       = new TH1D("recHlljj_AFTERZjetDeltaR",      "reconstructed 4-body invariant mass after #DeltaR_{jj} cut",          nbin_,0.,1000.);
  recHlljj_AFTERZjetShiftedEta   = new TH1D("recHlljj_AFTERZjetShiftedEta",  "reconstructed 4-body invariant mass after Shifted |#eta_{j}| cut",    nbin_,0.,1000.);
  recHlljjMassResolution_ = new TH1D("recHlljjMassResolution","reconstructed dilepton+dijet invariant mass resolution",  nbin_,-3.,3.);
  recHllbbMassResolution_ = new TH1D("recHllbbMassResolution","reconstructed dilepton+dib-jet invariant mass resolution",nbin_,-3.,3.);

  // histo from HEPG information
  // ---------------------------
  PtHParticles_ = new TH1D("PtHZZ", "p_{T} distribution of the ZZ associated to H",      nbin_,0.,1000.);
  PtHParticle1_ = new TH1D("PtHZ1", "p_{T} distribution of the 1^{st} Z associated to H",nbin_,0.,1000.);
  PtHParticle2_ = new TH1D("PtHZ2", "p_{T} distribution of the 2^{nd} Z associated to H",nbin_,0.,1000.);
  MHtoZZllqq_   = new TH1D("MHtoZZllqq","H mass from Zs (Z->ll & Z->qq) associated to H",nbin_,0.,1000.);
  MHtoZZllbb_   = new TH1D("MHtoZZllbb","H mass from Zs (Z->ll & Z->bb) associated to H",nbin_,0.,1000.);
  MZZ_          = new TH1D("MZZ",     "ZZ invariant mass",                                                   nbin_,0.,1000.);
  MZ1vsMZ2_     = new TH2D("MZ1vsMZ2","Z1 mass vs Z2 mass",                                                  nbin_,0.,200.,nbin_,0.,200.);
  M4_           = new TH1D("M4",      "Z+b\\bar{b} invariant mass",                                          nbin_,0.,1000.);
  M4irrBkg_     = new TH1D("M4irrBkg","Z+b\\bar{b} invariant mass for 80 GeV\\le m_{b\\bar{b}}\\le 100 GeV", nbin_,0.,1000.);
  Mbb_          = new TH1D("Mbb",     "b\\bar{b} invariant mass",                                            nbin_,0.,1000.);
  MZvsMbb_      = new TH2D("MZvsMbb", "Z mass vs b\\bar{b} invariant mass",                                  nbin_,0.,200.,nbin_,0.,250.);
  
  DeltaPhiZZ_ = new TH1D("DeltaPhiZZ", "#Delta#Phi between the 2 Z",nbin_,  0.,TMath::Pi());   
  DeltaEtaZZ_ = new TH1D("DeltaEtaZZ", "#Delta#eta between the 2 Z",nbin_,-10.,10.);   
  DeltaRZZ_   = new TH1D("DeltaRZZ",   "#DeltaR between the 2 Z",   nbin_,  0.,10.);   

  vecMZtoParticles_.push_back( TH1D("MZbbbar",    "Z mass from b-quarks associated to Z",  nbin_, 0., 200.));
  vecMZtoParticles_.push_back( TH1D("MZelectron", "Z mass from electrons associated to Z", nbin_, 0., 200.));
  vecMZtoParticles_.push_back( TH1D("MZmuons",    "Z mass from muons associated to Z",     nbin_, 0., 200.));

  vecPtZParticles_.push_back( TH1D("PtZbquarks",  "p_{T} distribution of the b-quarks associated to Z",           nbin_, 0., 1000.));
  vecPtZParticles_.push_back( TH1D("PtZelectrons","p_{T} distribution of the electrons associated to Z",          nbin_, 0., 1000.));
  vecPtZParticles_.push_back( TH1D("PtZmuons",    "p_{T} distribution of the muons associated to Z",              nbin_, 0., 1000.));
  vecPtZParticle1_.push_back( TH1D("PtZbquark1",  "p_{T} distribution of the first b-quark associated to Z",      nbin_, 0., 1000.));
  vecPtZParticle1_.push_back( TH1D("PtZelectron1","p_{T} distribution of the first electron associated to Z" ,    nbin_, 0., 1000.));
  vecPtZParticle1_.push_back( TH1D("PtZmuon1",    "p_{T} distribution of the first muon associated to Z",         nbin_, 0., 1000.));
  vecPtZParticle2_.push_back( TH1D("PtZbquark2",  "p_{T} distribution of the second b-quark associated to Z",     nbin_, 0., 1000.));
  vecPtZParticle2_.push_back( TH1D("PtZelectron2","p_{T} distribution of the second electron associated to Z",    nbin_, 0., 1000.));
  vecPtZParticle2_.push_back( TH1D("PtZmuon2",    "p_{T} distribution of the second muon associated to Z",        nbin_, 0., 1000.));
  vecDeltaPhiZParticles_.push_back( TH1D("DeltaPhiZbbbar",    "#Delta#Phi between the b-quarks associated to Z", nbin_, 0., TMath::Pi()));   
  vecDeltaPhiZParticles_.push_back( TH1D("DeltaPhiZelectrons","#Delta#Phi between the electrons associated to Z",nbin_, 0., TMath::Pi()));
  vecDeltaPhiZParticles_.push_back( TH1D("DeltaPhiZmuons",    "#Delta#Phi between the muons associated to Z",    nbin_, 0., TMath::Pi()));
  vecDeltaPhiZParticles_.push_back( TH1D("DeltaPhiZqqbar",    "#Delta#Phi between the quarks associated to Z",   nbin_, 0., TMath::Pi()));   
  vecDeltaEtaZParticles_.push_back( TH1D("DeltaEtaZbbbar",    "#Delta#eta between the b-quarks associated to Z", nbin_, -10., 10));   
  vecDeltaEtaZParticles_.push_back( TH1D("DeltaEtaZelectrons","#Delta#eta between the electrons associated to Z",nbin_, -10., 10));
  vecDeltaEtaZParticles_.push_back( TH1D("DeltaEtaZmuons",    "#Delta#eta between the muons associated to Z",    nbin_, -10., 10));
  vecDeltaEtaZParticles_.push_back( TH1D("DeltaEtaZqqbar",    "#Delta#eta between the quarks associated to Z",   nbin_, -10., 10));   
  vecDeltaRZParticles_.push_back( TH1D(  "DeltaRZbbbar",      "#DeltaR between the b-quarks associated to Z", nbin_, 0., 10.));   
  vecDeltaRZParticles_.push_back( TH1D(  "DeltaRZelectrons",  "#DeltaR between the electrons associated to Z",nbin_, 0., 10.));
  vecDeltaRZParticles_.push_back( TH1D(  "DeltaRZmuons",      "#DeltaR between the muons associated to Z",    nbin_, 0., 10.));
  vecDeltaRZParticles_.push_back( TH1D(  "DeltaRZqqbar",      "#DeltaR between the quarks associated to Z",   nbin_, 0., 10.));   

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsZZAnalyzer::endJob() {

  OutputFile->cd();
  // Fill histograms
  // ---------------
  eventsNumber_->SetBinContent(1,eventcounter_);
  eventsNumber_->Write();
  signalTopologySizeNumber_->SetBinContent(1,signalTopologySizeCounter_);
  signalTopologySizeNumber_->Write();
  failedTopologyNumber_->Write();
  HtoZZEventsNumber_->SetBinContent(1,HtoZZeventcounter_);
  HtoZZEventsNumber_->Write();

  muonNumber_->Write();
  goodMuonNumber_->Write();
  electronNumber_->Write();
  goodElectronNumber_->Write();
  goodIc5JetNumber_->Write();
  goodIc5centralJetNumber_->Write();
  Ic5forwardJetNumber_->Write();
  goodIc5forwardJetNumber_->Write();

  OfflinejetNumber_->Write();
  OfflinejetUncorrEt_->Write();
  OfflinejetEt_->Write();
  OfflinejetEtScale_->Write();
  OfflinejetEtScaleVSUncorrEt_->Write();
  OfflinejetEtScaleVSEt_->Write();
  OfflinejetEtScaleVSEta_->Write();
  OfflinejetEtVSUncorrEt_->Write();
  OfflinejetEtVSEt_->Write();
  OfflinejetEtVSEta_->Write();
  OfflinejetUncorrEtVSUncorrEt_->Write();
  OfflinejetUncorrEtVSEt_->Write();
  OfflinejetUncorrEtVSEta_->Write();
  OfflinejetEMenergyfrac_->Write();
  OfflinejetHighEffDiscr_->Write();
  OfflinejetHighPurDiscr_->Write();
  OfflinejetMass_->Write();
  OfflinejetTKnumS2_->Write();
  OfflinejetTKsumPtS2_->Write();
  OfflinejetTKtagmassS2_->Write();
  OfflinejetPhi_->Write();
  OfflinejetEta_->Write();
  OfflinejetLeptonDeltaR_->Write();
  centralJetLeptonDeltaR_->Write();

  pTbalance_->Write();
  pTbalancePhi_->Write();
  pTbalanceVSMEt_->Write();
  pTbalanceMEtDeltaPhi_->Write();

  MuonNumber_->Write();
  MuonEta_->Write();       
  MuonPhi_->Write();       
  MuonPt_->Write();        
  MuonPx_->Write();        
  MuonPy_->Write();        
  MuonPz_->Write();        
  MuonEMcalo_->Write();    
  MuonHADcalo_->Write();   
  MuonHOcalo_->Write();    
  MuonIso03SumPt_->Write();
  MuonNormChi2_->Write();  
  MuonHitsNum_->Write();   
  MuonImpactParamXY_->Write();
  MuonErrorImpactParamXY_->Write();
  MuonImpactParamSig_->Write();

  ElectronNumber_->Write();
  ElectronEt_->Write();	    
  ElectronEta_->Write();	    
  ElectronPhi_->Write();	    
  ElectronPt_->Write();	    
  ElectronPx_->Write();	    
  ElectronPy_->Write();	    
  ElectronPz_->Write();	    
  ElectronHadOverEm_->Write();
  ElectronIsoVal_->Write();   
  ElectronImpactParamXY_->Write();
  ElectronErrorImpactParamXY_->Write();
  ElectronImpactParamSig_->Write();

  MEt_->Write();
  MEtSumEt_->Write();
  MEtSignificance_->Write();
  MEtPhi_->Write();
  
  dileptonMassVSdeltaR_->Write();
  dijetMassVSdeltaR_->Write();   
  dibjetMassVSdeltaR_->Write();  
  scaleddijetmass_->Write();
  dijetMassComparison_->Write();
  dijetMassVSscaleddibjetMass_->Write();

  recZllMass_->Write();
  recZjjMass_->Write();
  recZbbMass_->Write();
  recZllMassResolution_->Write();
  recZeeMassResolution_->Write();
  recZmmMassResolution_->Write();
  recZjjMassResolution_->Write();
  recZbbMassResolution_->Write();
  recZllDeltaR_->Write();
  recZjjDeltaR_->Write();
  recZbbDeltaR_->Write();
  recZllDeltaPhi_->Write();
  recZjjDeltaPhi_->Write();
  recZbbDeltaPhi_->Write();
  recZllDeltaEta_->Write();
  recZjjDeltaEta_->Write();
  recZbbDeltaEta_->Write();
  recZll1ShiftedEta_->Write();
  recZjj1ShiftedEta_->Write();
  recZbb1ShiftedEta_->Write();
  recZll2ShiftedEta_->Write();
  recZjj2ShiftedEta_->Write();
  recZbb2ShiftedEta_->Write();
  recZllShiftedEta_->Write();
  recZjjShiftedEta_->Write();
  recZbbShiftedEta_->Write();

  recZjjLeptonDeltaR_->Write();
  recZbbLeptonDeltaR_->Write();

  forwardJetEt_->Write();
  forwardJetLeptonDeltaR_->Write();
  forwardJetEta_->Write();
  forwardJet1Eta_->Write();
  forwardJet2Eta_->Write();
  forwardDiJetMass_->Write();

  recHZZdeltaPhi_->Write();
  recHZZdeltaEta_->Write();
  recHZZdeltaR_->Write();
  recHlljj_->Write();	     
  recHlljj_forwardDiJetMass->Write();	     
  recHlljj_pTbalance->Write();	     
  recHlljj_ZdileptonMass->Write();	     
  recHlljj_leptonDeltaR->Write();	     
  recHlljj_leptonShiftedEta->Write();	     
  recHlljj_ZdijetMass->Write();	     
  recHlljj_ZjetDeltaR->Write();	     
  recHlljj_ZjetShiftedEta->Write();	     
  recHlljjVSforwardDiJetMass_->Write();	     
  recHlljjVSpTbalance_->Write();	     
  recHlljjVSZdileptonMass_->Write();	     
  recHlljjVSleptonDeltaR_->Write();	     
  recHlljjVSleptonShiftedEta_->Write();	     
  recHlljjVSZdijetMass_->Write();	     
  recHlljjVSZjetDeltaR_->Write();	     
  recHlljjVSZjetShiftedEta_->Write();	     
  recHlljj_AFTERforwardDiJetMass->Write();	     
  recHlljj_AFTERpTbalance->Write();	     
  recHlljj_AFTERZdileptonMass->Write();	     
  recHlljj_AFTERleptonDeltaR->Write();	     
  recHlljj_AFTERleptonShiftedEta->Write();	     
  recHlljj_AFTERZdijetMass->Write();	     
  recHlljj_AFTERZjetDeltaR->Write();	     
  recHlljj_AFTERZjetShiftedEta->Write();	     
  recHlljjMassResolution_->Write();	     
  recHllbbMassResolution_->Write();	     

  MHtoZZllqq_->Write();	     
  MHtoZZllbb_->Write();	     
  PtHParticles_->Write();	     
  PtHParticle1_->Write();	     
  PtHParticle2_->Write();	     
  MZ1vsMZ2_->Write();

  Mbb_->Write();
  M4_->Write();
  M4irrBkg_->Write();
  MZZ_->Write();
  MZvsMbb_->Write();

  DeltaPhiZZ_->Write();
  DeltaEtaZZ_->Write();
  DeltaRZZ_->Write();
  std::vector< std::pair<int,int> >::const_iterator vecTypeIndexPair_itr = vecTypeIndexPair.begin();
  for ( ; vecTypeIndexPair_itr != vecTypeIndexPair.end(); ++vecTypeIndexPair_itr ) {
    int index = vecTypeIndexPair_itr->second;
    vecZtoParticlesEventsNumber_[index].SetBinContent(1,(double)vecZtoParticleseventcounter_[index]);
    vecZtoParticlesEventsNumber_[index].Write();
    vecMZtoParticles_[index].Write();
    vecPtZParticles_[index].Write();
    vecPtZParticle1_[index].Write();
    vecPtZParticle2_[index].Write();
    vecDeltaPhiZParticles_[index].Write();
    vecDeltaEtaZParticles_[index].Write();
    vecDeltaRZParticles_[index].Write();
  }
  vecZtoParticlesEventsNumber_[3].SetBinContent(1,(double)vecZtoParticleseventcounter_[3]);
  vecZtoParticlesEventsNumber_[4].SetBinContent(1,(double)vecZtoParticleseventcounter_[4]);
  vecZtoParticlesEventsNumber_[3].Write();
  vecZtoParticlesEventsNumber_[4].Write();
  vecDeltaPhiZParticles_[3].Write();
  vecDeltaEtaZParticles_[3].Write();
  vecDeltaRZParticles_[3].Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsZZAnalyzer);
