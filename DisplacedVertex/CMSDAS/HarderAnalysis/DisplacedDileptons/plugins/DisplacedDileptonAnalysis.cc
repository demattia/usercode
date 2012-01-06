#include "HarderAnalysis/DisplacedDileptons/plugins/DisplacedDileptonAnalysis.h"
#include "HarderAnalysisTools/GenParticleInfo/interface/GenParticleInfo.h"
#include "HarderAnalysisTools/CutManager/interface/CutFlow.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/Common/interface/Vector.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TDCacheFile.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include <sstream>
#include <iomanip>

// PAT
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


DisplacedDileptonAnalysis::DisplacedDileptonAnalysis(const edm::ParameterSet& iConfig):
  checkHitPattern_(NULL),
  generatorTag_(iConfig.getParameter<edm::InputTag>("generatorSrc")),
  pileupTag_(edm::InputTag("addPileupInfo")),
  genEventInfoTag_(edm::InputTag("generator")),
  trigger_(iConfig.getParameter<edm::InputTag>("trigger")),
  triggerEvent_(iConfig.getParameter<edm::InputTag>("triggerEvent")),
  barrelSuperClusters_(iConfig.getParameter<edm::InputTag>("barrelSuperClusters")),
  endcapSuperClusters_(iConfig.getParameter<edm::InputTag>("endcapSuperClusters")),
  signalPDGId_(iConfig.getParameter<int>("signalPDGId")),
  leptonPDGId_(iConfig.getParameter<int>("leptonPDGId")),
  cutStudyMode_(iConfig.getParameter<bool>("cutStudyMode")),
  hltPaths_(iConfig.getParameter<std::vector<std::string> >("hltPaths")),
  runMin_(iConfig.getParameter<int>("minRunNumber")),
  runMax_(iConfig.getParameter<int>("maxRunNumber")),
  leptonPtCut_(iConfig.getParameter<double>("leptonPtCut")),
  leptonEtaCut_(iConfig.getParameter<double>("leptonEtaCut")),
  leptonIsolationCut_(iConfig.getParameter<double>("leptonIsolationCut")),
  leptonChargeCut_(iConfig.getParameter<bool>("leptonChargeCut")),
  vetoBackToBack_(iConfig.getParameter<double>("vetoBackToBack")),
  minDeltaRBetweenLeptons_(iConfig.getParameter<double>("minDeltaRBetweenLeptons")),
  vertexChi2Cut_(iConfig.getParameter<double>("vertexChi2Cut")),
  deltaPhiCut_(iConfig.getParameter<double>("deltaPhiCut")),
  decayLengthCut_(iConfig.getParameter<double>("decayLengthCut")),
  decayLengthSignificanceCut_(iConfig.getParameter<double>("decayLengthSignificanceCut")),
  minD0sig_(iConfig.getParameter<double>("minD0Significance")),
  maxHitsBeforeVertex_(iConfig.getParameter<int>("maxHitsBeforeVertex")),
  numTrigMatches_(iConfig.getParameter<int>("numTrigMatches")),
  maxNumStandAloneMuons_(iConfig.getParameter<int>("maxNumStandAloneMuons")),
  minNumCaloMatches_(iConfig.getParameter<int>("minNumCaloMatches")),
  numProcessedEvents_(0)
{
  if (leptonPDGId_==11) {
    thisLepton_=__lepTypeElectron;
    leptonName_="electron";
  } else if (leptonPDGId_==13) {
    thisLepton_=__lepTypeMuon;
    leptonName_="muon";
  } else if (leptonPDGId_==15) {
    thisLepton_=__lepTypeTau;
    leptonName_="tau";
  } else if (leptonPDGId_==-11) {
    thisLepton_=__lepTypeTrack;
    leptonName_="etrack";
  } else if (leptonPDGId_==-13) {
    thisLepton_=__lepTypeTrack;
    leptonName_="mutrack";
  } else if (leptonPDGId_==-15) {
    thisLepton_=__lepTypeTrack;
    leptonName_="tautrack";
  } else {
    throw cms::Exception("InvalidLeptonPDG") << "abs(leptonPDGId) must be 11, 13 or 15";
  }
  std::cout << "ANALYSIS CHANNEL: " << leptonName_ << std::endl;

  leptonCollTag_[__lepTypeElectron]=iConfig.getParameter<edm::InputTag>("electronSrc");
  leptonCollTag_[__lepTypeMuon]=iConfig.getParameter<edm::InputTag>("muonSrc");
  leptonCollTag_[__lepTypeTau]=iConfig.getParameter<edm::InputTag>("tauSrc");
  leptonCollTag_[__lepTypeTrack]=iConfig.getParameter<edm::InputTag>("trackSrc");

  // histogram managers  
  mcHistos_ = new HistMap("gen");
  weightHistos_ = new HistMap("weights");
  trigHistos_ = new HistMap("trig");
  lepHistos_ = new HistMap("leptons");
  dilepHistos_ = new HistMap("dileptons");

  // histogram boundaries
  histMaxRt_=100;
  histMaxM_=300;
  histMaxE_=1000;
  histMaxPt_=100;

  // number of bins
  histNBinsRt_=10;

  cutSummaryLeptonSignal_.clear();
  cutSummaryLeptonBackground_.clear();
  cutSummarySignal_.clear();
  cutSummaryPartial_.clear();
  cutSummaryBackground_.clear();

  // configure luminosity reweighting
  std::vector<double> lumiDistrMCdbl=iConfig.getParameter<std::vector<double> >("lumiDistrMC");
  std::vector<double> lumiDistrDatadbl=iConfig.getParameter<std::vector<double> >("lumiDistrData");
  std::vector<float> lumiDistrMC,lumiDistrData;
  for (unsigned i=0; i<std::min(lumiDistrMCdbl.size(),lumiDistrDatadbl.size()); i++) {
    lumiDistrMC.push_back(lumiDistrMCdbl[i]);
    lumiDistrData.push_back(lumiDistrDatadbl[i]);
  }
  LumiWeights_ = new edm::LumiReWeighting(lumiDistrMC,lumiDistrData);
  // for systematic uncertainty according to
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupSystematicErrors
  PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  PShiftUp_ = reweight::PoissonMeanShifter(0.6);
  PShiftUpUp_ = reweight::PoissonMeanShifter(1.2);

  // introduce short dCache read timeout to work around problems at RAL
  TDCacheFile::SetOpenTimeout(1800);
}


DisplacedDileptonAnalysis::~DisplacedDileptonAnalysis()
{
}


void DisplacedDileptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // event setup
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);
  if (checkHitPattern_==NULL) checkHitPattern_ = new CheckHitPattern(iSetup);

  bool isMC = findMCParticles(iEvent);

  // look at trigger information
  bool isTriggered=doTrigger(iEvent);

  // get beamspot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel ("offlineBeamSpot", beamSpotHandle);
  beamSpot_=(*beamSpotHandle);
  
  // get primary vertex
  edm::Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel ("offlinePrimaryVertices", primaryVertex);
  // there should always be a vertex in this collection. even if the primary vertex fit failed,
  // a fake vertex based on the beam spot should be present
  if (primaryVertex->size()==0) {
    std::cout << "no primary vertex found! skipping event" << std::endl;
    return;
  }
  reco::Vertex pv = primaryVertex->front();
  VertexState PVstate(GlobalPoint(pv.x(),pv.y(),pv.z()),GlobalError(pv.covariance()));
  primaryVertex_=&PVstate;
  // count number of proper primary vertices in event
  numPV_=0;
  for (unsigned i=0; i<primaryVertex->size(); i++) {
    if (!(*primaryVertex)[i].isFake()) ++numPV_;
  }

  // general track studies
  investigateTracks(iEvent);

  // evaluate generator-level event properties needed for reweighting,
  // such as PDF, pile-up vertices,
  // decay channels and lifetimes of long-lived particles etc
  GenEventProperties genProp(iEvent,signalPDGId_,generatorTag_,
			     pileupTag_,genEventInfoTag_);

  // calculate event weights
  weightMap_.clear();
  double lumiWeight=1;
  double lumiWeight_low=1;
  double lumiWeight_high=1;
  double lumiWeight_veryhigh=1;

  if (isMC) {
    // luminosity weight
    const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
    lumiWeight = LumiWeights_->weight( (*iEventB) );
    weightMap_["lumiWeight_central"]=lumiWeight;
    weightHistos_->fill("lumiWeight_central",lumiWeight,1,
			"lumiWeight_central",100,0,2);
    lumiWeight_low=lumiWeight
      *PShiftDown_.ShiftWeight( genProp.aveNumPileup3BX() );
    weightMap_["lumiWeight_low"]=lumiWeight_low;
    weightHistos_->fill("lumiWeight_low",lumiWeight_low,1,
			"lumiWeight_low",100,0,2);
    lumiWeight_high=lumiWeight
      *PShiftUp_.ShiftWeight( genProp.aveNumPileup3BX() );
    weightMap_["lumiWeight_high"]=lumiWeight_high;
    weightHistos_->fill("lumiWeight_high",lumiWeight_high,1,
			"lumiWeight_high",100,0,2);
    lumiWeight_veryhigh=lumiWeight
      *PShiftUpUp_.ShiftWeight( genProp.aveNumPileup3BX() );
    weightMap_["lumiWeight_veryhigh"]=lumiWeight_veryhigh;
    weightHistos_->fill("lumiWeight_veryhigh",lumiWeight_veryhigh,1,
			"lumiWeight_veryhigh",100,0,2);
  }
  if (isTriggered) {
    weightHistos_->fill("numPV",numPV_,1,
			"number of reconstructed primary vertices",
			25,0,25);
    weightHistos_->fill("numPV_central",numPV_,lumiWeight,
			"number of reconstructed primary vertices, reweighted",
			25,0,25);
    weightHistos_->fill("numPV_low",numPV_,lumiWeight_low,
			"number of reconstructed primary vertices, reweighted low",
			25,0,25);
    weightHistos_->fill("numPV_high",numPV_,lumiWeight_high,
			"number of reconstructed primary vertices, reweighted high",
			25,0,25);
    weightHistos_->fill("numPV_veryhigh",numPV_,lumiWeight_veryhigh,
			"number of reconstructed primary vertices, reweighted very high",
			25,0,25);
  }

  // dilepton reconstruction
  if (thisLepton_==__lepTypeElectron) {
    doDecayChannel<pat::Electron>(iEvent,__lepTypeElectron,genProp);
  } else if (thisLepton_==__lepTypeMuon) {
    doDecayChannel<pat::Muon>(iEvent,__lepTypeMuon,genProp);
  } else if (thisLepton_==__lepTypeTau) {
    doDecayChannel<pat::Tau>(iEvent,__lepTypeTau,genProp);
  } else {
    doDecayChannel<PseudoLepton>(iEvent,__lepTypeTrack,genProp);
  }

  ++numProcessedEvents_;
}


template<class LeptonCollection>
void DisplacedDileptonAnalysis::fillTrigMatches(const edm::Event& iEvent,leptonType lepType) {

  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEvent_, triggerEvent );
  if (triggerEvent.failedToGet()) {
    std::cout << "WARNING: cannot access triggerEvent for matching" << std::endl;
    return;
  }

  edm::Handle<LeptonCollection> leptons;
  iEvent.getByLabel(leptonCollTag_[lepType],leptons);

  leptonTriggerMatch_.clear();

  // loop over all lepton candidates
  for (unsigned ilep=0; ilep<leptons->size(); ilep++) {

    // skip leptons below p_t threshold
    if ((*leptons)[ilep].pt()<leptonPtCut_) {
      leptonTriggerMatch_.push_back(0);
      continue;
    }

    // find deltaR match among trigger objects
    double lep_eta=(*leptons)[ilep].eta();
    double lep_phi=(*leptons)[ilep].phi();
    double min_deltaR=9999.;
    const reco::Particle::LorentzVector* trig_mom=0;

    // loop over all triggers we are interested in
    for (unsigned itrig=0; itrig<hltPaths_.size(); itrig++) {
      if (!triggerEvent->path(hltPaths_[itrig])) continue;
      if (!triggerEvent->path(hltPaths_[itrig])->wasAccept()) continue;
      // get all trigger objects associated with this HLT path
      const pat::TriggerObjectRefVector triggerObjects(triggerEvent->pathObjects(hltPaths_[itrig],false));
      for (unsigned i=0; i<triggerObjects.size(); i++) {
	double trig_eta=triggerObjects[i]->eta();
	double trig_phi=triggerObjects[i]->phi();
	double dR=deltaR(lep_eta,lep_phi,trig_eta,trig_phi);
	if (dR<min_deltaR) {
	  min_deltaR=dR;
	  trig_mom=&(triggerObjects[i]->p4());
	}
      }
    }
    trigHistos_->fill(leptonName_+"_TrigMatchDeltaR",min_deltaR,1,
		      leptonName_+" deltaR distance to nearest trigger object",100,0,1.6);

    bool trigMatch=(min_deltaR<0.1);
    if (trigMatch) {
      leptonTriggerMatch_.push_back(trig_mom);
      if (trig_mom!=0) {
	trigHistos_->fill(leptonName_+"_ptratio",(*leptons)[ilep].pt()/trig_mom->pt(),1,
			  "ratio of track p_t over trigger p_t",100,0,2);
      }
    } else {
      leptonTriggerMatch_.push_back(0);
    }
    trigHistos_->fill(leptonName_+"_hasTrigMatch",trigMatch,1,
		      "does "+leptonName_+" have a trigger match",2,0,2);
  }
}


void DisplacedDileptonAnalysis::fillCaloMatches(const edm::Event& iEvent) {

  edm::Handle< std::vector<reco::SuperCluster> > barrelSuperClusters;
  iEvent.getByLabel( barrelSuperClusters_, barrelSuperClusters );
  if (barrelSuperClusters.failedToGet()) {
    std::cout << "WARNING: cannot access barrel superclusters for matching" << std::endl;
    return;
  }
  edm::Handle< std::vector<reco::SuperCluster> > endcapSuperClusters;
  iEvent.getByLabel( endcapSuperClusters_, endcapSuperClusters );
  if (endcapSuperClusters.failedToGet()) {
    std::cout << "WARNING: cannot access endcap superclusters for matching" << std::endl;
    return;
  }

  edm::Handle<PseudoLeptonCollection> leptons;
  iEvent.getByLabel(leptonCollTag_[__lepTypeTrack],leptons);

  leptonCaloMatch_.clear();

  // loop over all lepton candidates
  for (unsigned ilep=0; ilep<leptons->size(); ilep++) {

    // skip leptons below p_t threshold
    if ((*leptons)[ilep].pt()<leptonPtCut_) {
      leptonCaloMatch_.push_back(0);
      continue;
    }

    // find deltaR match among superclusters
    double lep_eta=(*leptons)[ilep].eta();
    double lep_phi=(*leptons)[ilep].phi();
    double min_deltaR=9999.;
    double calo_energy=0;

    const double deltaRcut=0.2;

    // loop over all barrel superclusters
    if (fabs(lep_eta)<1.48+1.5*deltaRcut) {
      for (unsigned iclus=0; iclus<barrelSuperClusters->size(); iclus++) {
	double clus_eta=(*barrelSuperClusters)[iclus].eta();
	double clus_phi=(*barrelSuperClusters)[iclus].phi();
	double dR=deltaR(lep_eta,lep_phi,clus_eta,clus_phi);
	if (dR<min_deltaR) {
	  min_deltaR=dR;
	  calo_energy=(*barrelSuperClusters)[iclus].energy();
	}
      }
    }
    // loop over all endcap superclusters
    if (fabs(lep_eta)>1.48-1.5*deltaRcut) {
      for (unsigned iclus=0; iclus<endcapSuperClusters->size(); iclus++) {
	double clus_eta=(*endcapSuperClusters)[iclus].eta();
	double clus_phi=(*endcapSuperClusters)[iclus].phi();
	double dR=deltaR(lep_eta,lep_phi,clus_eta,clus_phi);
	if (dR<min_deltaR) {
	  min_deltaR=dR;
	  calo_energy=(*endcapSuperClusters)[iclus].energy();
	}
      }
    }

    lepHistos_->fill(leptonName_+"_CaloMatchDeltaR",min_deltaR,1,
		     leptonName_+" deltaR distance to nearest supercluster",100,0,1.6);

    bool caloMatch=(min_deltaR<deltaRcut);
    if (caloMatch) {
      leptonCaloMatch_.push_back(calo_energy);
      if (calo_energy!=0) {
	lepHistos_->fill(leptonName_+"_Eratio",(*leptons)[ilep].p()/calo_energy,1,
			  "ratio of track momentum over supercluster energy",100,0,2);
      }
    } else {
      leptonCaloMatch_.push_back(0);
    }
    lepHistos_->fill(leptonName_+"_hasCaloMatch",caloMatch,1,
		      "does "+leptonName_+" have a supercluster match",2,0,2);
  }
}


bool DisplacedDileptonAnalysis::doTrigger(const edm::Event& iEvent) {

  // PAT trigger information
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEvent_, triggerEvent );
  trigHistos_->fill("failed",triggerEvent.failedToGet(),1,
		    "failed to retrieve TriggerEvent",2,0,2);
  if (triggerEvent.failedToGet()) {
    std::cout << "WARNING: no TriggerEvent found" << std::endl;
    return false;
  }

  edm::Handle< pat::TriggerPathCollection > triggerPaths;
  iEvent.getByLabel( trigger_, triggerPaths );
  edm::Handle< pat::TriggerFilterCollection > triggerFilters;
  iEvent.getByLabel( trigger_, triggerFilters );
  edm::Handle< pat::TriggerObjectCollection > triggerObjects;
  iEvent.getByLabel( trigger_, triggerObjects );

  trigHistos_->fill("wasRun",triggerEvent->wasRun(),1,"was trigger run",2,0,2);
  trigHistos_->fill("wasAccept",triggerEvent->wasAccept(),1,"was trigger accepted",2,0,2);
  trigHistos_->fill("wasError",triggerEvent->wasError(),1,"did trigger return error",2,0,2);
  trigHistos_->fill("eventsVsRun",iEvent.run(),1,"number of events found per run",
		    runMax_-runMin_+1,runMin_,runMax_);

  // analyze trigger objects in order to simulate triggers that weren't run in MC.
  // note that this will only work for un-prescaled triggers, i.e. while it will
  // hopefully be fine in MC, it will be very difficult to verify this in data.
  // (we could however do verification runs in special MC samples with newer trigger menus!)
  std::vector<std::string> pseudoTriggers;
  for (unsigned itrig=0; itrig<hltPaths_.size(); itrig++) {
    if (!triggerEvent->path(hltPaths_[itrig])) continue;
    if (!triggerEvent->path(hltPaths_[itrig])->wasAccept()) continue;
    // get all trigger objects associated with this HLT path
    const pat::TriggerObjectRefVector triggerObjects(triggerEvent->pathObjects(hltPaths_[itrig],false));
    trigHistos_->fill("numTriggerObjects_"+hltPaths_[itrig],triggerObjects.size(),1,
		      "number of trigger objects associated with "+hltPaths_[itrig],10,0,10);
    // count objects above various pt thresholds
    unsigned n30=0;
    unsigned n38=0;
    unsigned n43=0;
    for (unsigned i=0; i<triggerObjects.size(); i++) {
      double value=0;
      if (triggerObjects[i]->hasFilterId(trigger::TriggerPhoton)) {
	value=triggerObjects[i]->et();
      } else if (triggerObjects[i]->hasFilterId(trigger::TriggerMuon)) {
	value=triggerObjects[i]->pt();
      } else {
	// these are the L1 objects also associated with this trigger. ignore for now.
	continue;
      }
      if (value>30) ++n30;
      if (value>38) ++n38;
      if (value>43) ++n43;
    }

    if ((hltPaths_[itrig]=="HLT_L2DoubleMu23_NoVertex_v1") && n30>=2)
      pseudoTriggers.push_back("MC_L2DoubleMu30_NoVertex_v1");

    if ((hltPaths_[itrig]=="HLT_DoublePhoton33_v2") && n38>=2)
      pseudoTriggers.push_back("MC_DoublePhoton38_v2");
    if ((hltPaths_[itrig]=="HLT_DoublePhoton33_v2") && n43>=2)
      pseudoTriggers.push_back("MC_DoublePhoton43_v2");

  }

  // get list of trigger matches
  if (thisLepton_==__lepTypeElectron) {
    fillTrigMatches<pat::ElectronCollection>(iEvent,__lepTypeElectron);
  } else if (thisLepton_==__lepTypeMuon) {
    fillTrigMatches<pat::MuonCollection>(iEvent,__lepTypeMuon);
  } else if (thisLepton_==__lepTypeTau) {
    fillTrigMatches<pat::TauCollection>(iEvent,__lepTypeTau);
  } else {
    fillTrigMatches<PseudoLeptonCollection>(iEvent,__lepTypeTrack);
  }

  // now loop over all paths etc and fill histogram for each one that was found
  triggers_.clear();
  for (unsigned i=0; i<triggerEvent->acceptedPaths().size(); i++) {
    std::string trigname(triggerEvent->acceptedPaths()[i]->name());
    ++trigPathSummary_[trigname];
    for (unsigned itrg=0; itrg<hltPaths_.size(); itrg++) {
      if (trigname==hltPaths_[itrg]) {
  	triggers_.push_back(trigname);
  	trigHistos_->fill(trigname,iEvent.run(),1,trigname+" triggers per run",
  			  runMax_-runMin_+1,runMin_,runMax_);
      }
    }
  }
  // add all pseudotriggers
  for (unsigned i=0; i<pseudoTriggers.size(); i++) {
    triggers_.push_back(pseudoTriggers[i]);
  }
  // add a dummy trigger for each event (so we get candidates before any particular trigger)
  triggers_.push_back("anyTrigger");

  // list of accepted filters (not even sure I need this for anything)
  //for (pat::TriggerFilterRefVectorIterator it=triggerEvent->acceptedFilters().begin();
  //     it!=triggerEvent->acceptedFilters().end(); it++) {
  //  ++trigFilterSummary_[(*it)->label()];
  //}

  // return value: did any of our analysis triggers fire
  return (triggers_.size()>1);
}


void DisplacedDileptonAnalysis::investigateTracks(const edm::Event& iEvent) {

  edm::Handle< reco::TrackCollection > trackCollection;
  iEvent.getByLabel( "generalTracks", trackCollection );
  if (trackCollection.failedToGet()) {
    std::cout << "WARNING: generalTracks not found" << std::endl;
    return;
  }

  lepHistos_->fill("numTracks",trackCollection->size(),1,
		   "number of reconstructed general tracks in event",
		   100,0,100);
  if (trackCollection->size()>0) {
    lepHistos_->fill("trueLeptonRadWithTrack",smallestLeptonRadius_,1,
		     "true lepton radius with track (only for single particle events)",
		     100,0,histMaxRt_);
  }
  for (reco::TrackCollection::const_iterator trk=trackCollection->begin();
       trk!=trackCollection->end(); trk++) {
    lepHistos_->fill("allTrackD0",(*trk).d0(),1,"track 2d impact parameter",
		     100,0,histMaxRt_);
    lepHistos_->fill("trackD0Resolution",(*trk).d0(),smallestLeptonRadius_,1,
		     "reconstructed vs generated track impact parameter",
		     100,0,histMaxRt_,100,0,histMaxRt_);
    lepHistos_->fill("trackPtVsTrueD0",(*trk).pt(),smallestLeptonRadius_,1,
		     "reconstructed pt vs generated track impact parameter",
		     100,0,histMaxRt_,100,0,histMaxRt_);
    lepHistos_->fill("trackAlgoVsRadius",smallestLeptonRadius_,(*trk).algo(),1,
		     "track algorithm vs radius",100,0,50,
		     reco::TrackBase::algoSize,0,reco::TrackBase::algoSize);
    //NOTE: innerPosition() in TrackExtra -> not in AOD
    //lepHistos_->fill("trackInnerHitVsRadius",smallestLeptonRadius_,
    //		     (*trk).innerPosition().R(),1,
    //    	     "radius of innermost hit as function of impact parameter",
    //		     100,0,histMaxRt_,100,0,histMaxRt_);
    //lepHistos_->fill("trackInnerPhiVsRadius",smallestLeptonRadius_,
    //		     (*trk).innerPosition().Phi(),1,
    //		     "phi of innermost hit as function of impact parameter",
    //		     100,0,histMaxRt_,100,0,1.6);
    //lepHistos_->fill("trackInnerHitXY",(*trk).innerPosition().X(),
    //		     (*trk).innerPosition().Y(),1,
    //		     "xy position of innermost hit",
    //		     100,0,50,100,0,50);
    const reco::HitPattern& hp = (*trk).hitPattern();
    for (int i=0; i<hp.numberOfHits(); i++) {
      unsigned hit = hp.getHitPattern(i);
      if (hp.validHitFilter(hit) && hp.stripTOBHitFilter(hit)) {
      }
      lepHistos_->fill("trackHitPatternVsRadius",smallestLeptonRadius_,
		       hp.getLayer(hit),1,
		       "hit pattern in TOB",100,0,histMaxRt_,
		       10,0,10);
    }
  }
}


bool DisplacedDileptonAnalysis::findMCParticles(const edm::Event& iEvent) {

  edm::Handle<edm::View<reco::GenParticle> > mcParticles;
  iEvent.getByLabel(generatorTag_,mcParticles);
  if (mcParticles.failedToGet()) return false;

  unsigned nSignal=0;
  unsigned nGenLepton=0;
  smallestLeptonRadius_=9999;
  std::vector<double> signalLengths;
  std::vector<double> signalPhis;
  std::vector<double> signalEtas;
  std::vector<double> signalEnergies;
  std::vector<int> signalDecayChannels;
  for (unsigned imc=0; imc<mcParticles->size(); imc++) {

    // look for our signal particle
    if ((abs((*mcParticles)[imc].pdgId())==signalPDGId_) &&
	((*mcParticles)[imc].status()==3)) {
      ++nSignal;
      const reco::GenParticle part=(*mcParticles)[imc];
      GenParticleInfo genInfo(&part);
      int decayMode = decayChannel(part);
      mcHistos_->fill("gen_decaymode",decayMode,1,"decay channel",16,0,16);
      if (decayMode==abs(leptonPDGId_)) {
	mcHistos_->fill("gen_mass",part.mass(),1,"mass of signal particle, "
			+leptonName_+" channel",100,0,histMaxM_);
	mcHistos_->fill("gen_energy",part.energy(),1,"energy of signal particle, "
			+leptonName_+" channel",100,0,histMaxE_);
	mcHistos_->fill("gen_eta",part.eta(),1,"pseudo-rapidity of signal particle, "
			+leptonName_+" channel",100,-5,5);
	mcHistos_->fill("gen_phi",part.phi(),1,"phi of signal particle momentum, "
			+leptonName_+" channel",100,-5,5);
	mcHistos_->fill("gen_Ediff32",genInfo.fractionalDiffE(),1,
			"fractional energy difference signal particle status 3 vs status 2, "
			+leptonName_+" channel",100,-0.1,0.1);
	mcHistos_->fill("gen_numdaughters",genInfo.nDaughters(),1,
			"number of daughter particles, "+leptonName_+" channel",10,0,10);
	mcHistos_->fill("trueFlightDistance3D",genInfo.flightLength3D(),1,
			"3d flight length of signal particle, "+leptonName_+" channel",
			histNBinsRt_,0,histMaxRt_);
	mcHistos_->fill("trueFlightDistance2D",genInfo.flightLength2D(),1,
			"2d flight length of signal particle, "+leptonName_+" channel",
			histNBinsRt_,0,histMaxRt_);
	// properties of decay products
	if (genInfo.nDaughters()==2) {
	  mcHistos_->fill("gen_etaLep1Lep2",part.daughter(0)->eta(),part.daughter(1)->eta(),1,
			  "eta of signal lepton pair",100,-5,5,100,-5,5);
	  int nEB=0;
	  if (fabs(part.daughter(0)->eta())<1.442) ++nEB;
	  if (fabs(part.daughter(1)->eta())<1.442) ++nEB;
	  mcHistos_->fill("gen_nLepInBarrel",nEB,1,
			  "number of signal decay products in ECAL barrel region",3,0,3);
	}
      }
      // store info for correlation plots
      signalLengths.push_back(genInfo.flightLength2D());
      signalPhis.push_back(part.phi());
      signalEtas.push_back(part.eta());
      signalEnergies.push_back(part.energy());
      signalDecayChannels.push_back(decayMode);
    } else if ((abs((*mcParticles)[imc].pdgId())==abs(leptonPDGId_))
	       && ((*mcParticles)[imc].status()!=3)) {
      ++nGenLepton;
      double rad3d=(*mcParticles)[imc].vertex().R();
      double rad2d=(*mcParticles)[imc].vertex().Rho();
      if (rad2d<smallestLeptonRadius_) smallestLeptonRadius_=rad2d;
      mcHistos_->fill("leptonProdVtxRadius2D",rad2d,1,
		      "2d radius of production vertex of true "+leptonName_+"s",
		      100,0,histMaxRt_);
      mcHistos_->fill("leptonProdVtxRadius3D",rad3d,1,
		      "3d radius of production vertex of true "+leptonName_+"s",
		      100,0,histMaxRt_);
      if (signalOrigin(&(*mcParticles)[imc])) {
	mcHistos_->fill("leptonProdVtxRadius2D_signal",rad2d,1,
			"2d radius of production vertex of true "+leptonName_+"s from signal",
			100,0,histMaxRt_);
	mcHistos_->fill("leptonProdVtxRadius3D_signal",rad3d,1,
			"3d radius of production vertex of true "+leptonName_+"s from signal",
			100,0,histMaxRt_);
      }
    }
  }
  mcHistos_->fill("gen_num_signal",nSignal,1,
		  "number of signal particles in each event",5,0,5);
  mcHistos_->fill("gen_num_leptons",nGenLepton,1,
		  "number of true "+leptonName_+"s in each event",10,0,10);
  if (nSignal>1) {
    mcHistos_->fill("gen_correl_flightLength2D",signalLengths[0],signalLengths[1],1,
		    "flight length of two signal particles in event",
		    histNBinsRt_,0,histMaxRt_,histNBinsRt_,0,histMaxRt_);
    mcHistos_->fill("gen_correl_phi",signalPhis[0],signalPhis[1],1,
		    "phi of two signal particles in event",
		    100,-5,5,100,-5,5);
    mcHistos_->fill("gen_correl_eta",signalEtas[0],signalEtas[1],1,
		    "eta of two signal particles in event",
		    100,-5,5,100,-5,5);
    mcHistos_->fill("gen_correl_energy",signalEnergies[0],signalEnergies[1],1,
		    "energy of two signal particles in event",
		    100,0,histMaxE_,100,0,histMaxE_);
    mcHistos_->fill("gen_correl_decayChannel",signalDecayChannels[0],signalDecayChannels[1],1,
		    "decay channels of two signal particles in event",
		    5,11,16,5,11,16);
  }
  return true;
}


const reco::GenParticle* DisplacedDileptonAnalysis::getGenParticleMatch(const edm::Event& iEvent, const reco::Particle::LorentzVector& mom) {

  edm::Handle<edm::View<reco::GenParticle> > mcParticles;
  iEvent.getByLabel(generatorTag_,mcParticles);
  if (mcParticles.failedToGet()) {
    return NULL;
  }

  double minDeltaR=999;
  unsigned genIndex=0;
  for (unsigned imc=0; imc<mcParticles->size(); imc++) {
    reco::GenParticle mcPart=(*mcParticles)[imc];
    if (mcPart.status()==3) continue;
    if (mcPart.charge()==0) continue;
    if (mcPart.numberOfDaughters()>0) continue;
    double dR=deltaR(mom.eta(),mom.phi(),mcPart.eta(),mcPart.phi());
    if (dR<minDeltaR) {
      minDeltaR=dR;
      genIndex=imc;
    }
  }
  if (minDeltaR<0.1) {
    return &((*mcParticles)[genIndex]);
  } else {
    return NULL;
  }
}


double DisplacedDileptonAnalysis::trackerIsolation(const edm::Event& iEvent,
						   const reco::Particle::LorentzVector& plept,
						   const reco::Particle::LorentzVector& pveto) {

  // get PseudoLepton collection. these are high quality tracks
  // and standalone muons. we will only use the tracks for isolation, not muons.
  edm::Handle< PseudoLeptonCollection > trackCollection;
  iEvent.getByLabel( leptonCollTag_[__lepTypeTrack], trackCollection );
  if (trackCollection.failedToGet()) {
    std::cout << "WARNING: PseudoLepton collection not found" << std::endl;
    return -1;
  }

  // loop over all tracks and get deltaR with respect to both leptons
  double sumPt=0;
  for (PseudoLeptonCollection::const_iterator trk=trackCollection->begin();
       trk!=trackCollection->end(); trk++) {
    if (trk->isStandAloneMuon()) continue;
    double deltaRlept=deltaR(plept.eta(),plept.phi(),trk->eta(),trk->phi());
    double deltaRveto=deltaR(pveto.eta(),pveto.phi(),trk->eta(),trk->phi());
    if (deltaRlept<0.3 && deltaRlept>0.03 && deltaRveto>0.03) {
      // inside isolation code, but not close enough to any of the two tracks to suspect
      // that this might actually be one of the two tracks
      sumPt+=trk->pt();
    }
  }
  return sumPt;
}
  

reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const PseudoLepton &track) {
  reco::TransientTrack ttrack=trackBuilder_->build(track);
  ttrack.setBeamSpot(beamSpot_);
  return ttrack;
}


reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Electron &electron) {
  if (!electron.gsfTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(electron.gsfTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  } else {
    return 0;
  }
}


reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Muon &muon) {
  if (!muon.innerTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(muon.innerTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  } else if (!muon.outerTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(muon.outerTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  } else {
    return 0;
  }
}


reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Tau &tau) {
  if (!tau.leadTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(tau.leadTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  } else {
    return 0;
  }
}


// lepton IDs are specific to lepton type, but we provide a default template here
template<class Lepton>
bool DisplacedDileptonAnalysis::leptonID(const Lepton& lepton) {
  return true;
}

// override template for the electron channel, where we do use a specific ID requirement
bool DisplacedDileptonAnalysis::leptonID(const pat::Electron& electron) {
  return int(electron.electronID("eidLoose"))&1;
}


double DisplacedDileptonAnalysis::leptonTiming(const PseudoLepton &track) {
  return -999;
}


double DisplacedDileptonAnalysis::leptonTiming(const pat::Electron &electron) {
  //DetId seedDetID = electron->superCluster()->seed()->seed();
  return -999;
}


double DisplacedDileptonAnalysis::leptonTiming(const pat::Muon &muon) {
  if (muon.isTimeValid()) {
    // use time under the assumption the muon came from inside of the detector
    // and had beta==1.
    return muon.time().timeAtIpInOut;
  } else {
    return -999;
  }
}


double DisplacedDileptonAnalysis::leptonTiming(const pat::Tau &tau) {
  // not implemented
  return -999;
}


template<class Lepton>
bool DisplacedDileptonAnalysis::leptonCuts(const Lepton& lepton, CutFlow& cuts) {

  // basic kinematic and ID cuts
  if (leptonChargeCut_) cuts.applyCut("chargedTrackRequired",fabs(lepton.charge())>0);
  cuts.applyCut("leptonID",leptonID(lepton));
  cuts.applyLifetimeCut("leptonT0",leptonTiming(lepton),-99999,99999);
  return cuts.passAll();
}


template<class Lepton>
void DisplacedDileptonAnalysis::doDecayChannel(const edm::Event& iEvent,
					       const leptonType lepType,
					       GenEventProperties& genProp) {

  edm::Handle<edm::View<Lepton> > particles;
  iEvent.getByLabel(leptonCollTag_[lepType],particles);

  lepHistos_->fill("number",particles->size(),1,
		   "number of "+leptonName_+" candidates in event",20,0,20);

  // generator level information relevant for evnet weights
  int numDecays=0;
  for (unsigned i=0; i<genProp.numDecays(); i++) {
    if (genProp.decayMode(i)==abs(leptonPDGId_)) ++numDecays;
  }
  double ctau1=-1.0;
  if (genProp.numDecays()>0) ctau1=genProp.ctau(0);
  double ctau2=-1.0;
  if (genProp.numDecays()>1) ctau2=genProp.ctau(1);

  // make sure we have cut flow analysis for all triggers in this event
  for (unsigned itrig=0; itrig<triggers_.size(); itrig++) {
    if (!cutSummaryLeptonSignal_[triggers_[itrig]])
      cutSummaryLeptonSignal_[triggers_[itrig]]=new CutSummary("leptons_signal_"+triggers_[itrig]);
    if (!cutSummaryLeptonBackground_[triggers_[itrig]])
      cutSummaryLeptonBackground_[triggers_[itrig]]=new CutSummary("leptons_background_"+triggers_[itrig]);
    if (!cutSummarySignal_[triggers_[itrig]])
      cutSummarySignal_[triggers_[itrig]]=new CutSummary("dileptons_signal_"+triggers_[itrig]);
    if (!cutSummaryPartial_[triggers_[itrig]])
      cutSummaryPartial_[triggers_[itrig]]=new CutSummary("dileptons_partial_"+triggers_[itrig]);
    if (!cutSummaryBackground_[triggers_[itrig]])
      cutSummaryBackground_[triggers_[itrig]]=new CutSummary("dileptons_background_"+triggers_[itrig]);
  }

  // fill calo supercluster matches for etrack analysis
  if (leptonPDGId_==-11) {
    fillCaloMatches(iEvent);
  }

  // loop over leptons once and apply cuts.
  // we do this in a separate loop before the dilepton finding in order
  // to facilitate cut efficiency measurements, and also in order to avoid
  // applying the same cuts again and again while looping over all possible pairs
  std::vector<const reco::Candidate*> ancestor;
  std::vector<const reco::GenParticle*> genParticle;
  std::vector<bool> passCuts;
  std::vector<bool> passLifetimeCuts;
  std::vector<double> leptonD0;
  std::vector<double> leptonD0significance;
  unsigned nLeptonsPassCuts=0;
  for (unsigned ipart1=0; ipart1<particles->size(); ipart1++) {
    const Lepton& part1 = (*particles)[ipart1];

    // generator level information
    const reco::GenParticle* genPart=(getGenParticleMatch(iEvent,part1.p4()));
    genParticle.push_back(genPart);
    const reco::Candidate* origin=NULL;
    if (genPart) origin = signalOrigin(genPart);
    ancestor.push_back(origin);

    // selection cuts
    CutFlow leptonCutFlow;
    bool cutResult=leptonCuts(part1,leptonCutFlow);
    passCuts.push_back(leptonCutFlow.passAllIgnoreLifetime());
    passLifetimeCuts.push_back(leptonCutFlow.passAll());

    // impact parameter
    double d0=-9999;
    double d0significance=-9999;
    try {
      d0=leptonTrack(part1).stateAtBeamLine().transverseImpactParameter().value();
      d0significance=leptonTrack(part1).stateAtBeamLine().transverseImpactParameter().significance();
  } catch (cms::Exception& e) {};

    leptonD0.push_back(d0);
    leptonD0significance.push_back(d0significance);
    if (cutResult) ++nLeptonsPassCuts;
    for (unsigned itrig=0; itrig<triggers_.size(); itrig++) {
      if (origin) {
	cutSummaryLeptonSignal_[triggers_[itrig]]->addEntry(leptonCutFlow,
							    weightMap_,
							    numDecays,
							    ctau1,
							    ctau2,
							    leptonD0[ipart1]);
      } else {
	cutSummaryLeptonBackground_[triggers_[itrig]]->addEntry(leptonCutFlow,
								weightMap_,
								numDecays,
								ctau1,
								ctau2,
								leptonD0[ipart1]);
      }
    }

    // plots
    leptonPlots<Lepton>(&part1,genPart,lepType,leptonCutFlow.passAll(),(origin!=NULL));
  }

  // plot number of leptons passing all cuts in events with specific triggers
  for (unsigned itrig=0; itrig<triggers_.size(); itrig++) {
    lepHistos_->fill("nLeptonsPassAllCuts_"+triggers_[itrig],nLeptonsPassCuts,1,
		     "number of leptons passing all cuts in events with trigger "+triggers_[itrig],
		     10,0,10);
  }

  // loop over particle pairs of sufficiently high p_t
  unsigned ncand=0;
  for (unsigned ipart1=0; ipart1<particles->size(); ipart1++) {
    const Lepton part1 = (*particles)[ipart1];
    if (part1.pt()<leptonPtCut_) continue;
    for (unsigned ipart2=ipart1+1; ipart2<particles->size(); ipart2++) {
      const Lepton part2 = (*particles)[ipart2];
      if (part2.pt()<leptonPtCut_) continue;

      CutFlow cuts;

      // first set of cuts (before vertex reconstruction attempt)
      const Lepton* highPtLepton=&part1;
      const Lepton* lowPtLepton=&part2;
      unsigned iHighPt=ipart1;
      unsigned iLowPt=ipart2;
      if (part1.pt()<part2.pt()) {
	highPtLepton=&part2;
	lowPtLepton=&part1;
	iHighPt=ipart2;
	iLowPt=ipart1;
      }
      cuts.applyCut("leptonQualityL",passCuts[iLowPt]);
      cuts.applyCut("leptonQualityH",passCuts[iHighPt]);
      cuts.applyLifetimeCut("leptonAbsD0significanceL",
			    fabs(leptonD0significance[iLowPt]),minD0sig_,99999);
      cuts.applyLifetimeCut("leptonAbsD0significanceH",
			    fabs(leptonD0significance[iHighPt]),minD0sig_,99999);
      cuts.applyLifetimeCut("leptonD0L",leptonD0[iLowPt],-99999,99999);
      cuts.applyLifetimeCut("leptonD0H",leptonD0[iHighPt],-99999,99999);
      cuts.applyLifetimeCut("leptonD0significanceL",
			    leptonD0significance[iLowPt],-99999,99999);
      cuts.applyLifetimeCut("leptonD0significanceH",
			    leptonD0significance[iHighPt],-99999,99999);
      cuts.applyCut("leptonPtL",lowPtLepton->pt(),leptonPtCut_,999999.);
      cuts.applyCut("leptonPtH",highPtLepton->pt(),leptonPtCut_,999999.);
      cuts.applyCut("leptonEtaL",lowPtLepton->eta(),-leptonEtaCut_,leptonEtaCut_);
      cuts.applyCut("leptonEtaH",highPtLepton->eta(),-leptonEtaCut_,leptonEtaCut_);
      if (leptonChargeCut_) cuts.applyCut("oppositeCharge",(part1.charge()*part2.charge()<0));
      double isol1=trackerIsolation(iEvent,lowPtLepton->p4(),highPtLepton->p4());
      double isol2=trackerIsolation(iEvent,highPtLepton->p4(),lowPtLepton->p4());
      cuts.applyCut("trackerIsolationL",isol1,-999,leptonIsolationCut_);
      cuts.applyCut("trackerIsolationH",isol2,-999,leptonIsolationCut_);

      // require matching to trigger objects
      int nTrigMatches=0;
      if (leptonTriggerMatch_[ipart1]!=0) ++nTrigMatches;
      if (leptonTriggerMatch_[ipart2]!=0) ++nTrigMatches;
      cuts.applyCut("numTrigMatches",nTrigMatches,numTrigMatches_,2);

      // make sure we actually match to *different* trigger objects!
      // (this will kill very low mass signal particles, but it will
      // also remove candidates from split tracks and very close-by tracks
      cuts.applyCut("differentTrigObjects",nTrigMatches<2 ||
		    leptonTriggerMatch_[ipart1]!=leptonTriggerMatch_[ipart2]);

      // cosmic veto: reject tracks that are exactly back-to-back
      // (code adapted from S.E. Stoynev)
      reco::TrackBase::Vector mom1 = leptonTrack(part1).track().momentum();
      reco::TrackBase::Vector mom2 = leptonTrack(part2).track().momentum();
      double cosine = mom1.Dot(mom2)/mom1.R()/mom2.R();
      cuts.applyCut("vetoBackToBack",cosine,vetoBackToBack_,9999.);

      // require minimum deltaR between the leptons. This is especially
      // useful in the muon channel, where the trigger efficiency seems to
      // be difficult to understand in that region
      double dR=deltaR<reco::TrackBase::Vector,reco::TrackBase::Vector>(mom1,mom2);
      cuts.applyCut("deltaRBetweenLeptons",dR,minDeltaRBetweenLeptons_,9999);

      // skip the rest if not all cuts were passed so far and we want to save some CPU time
      if (!cutStudyMode_ && !cuts.passAll()) continue;

     // attempt vertex fit with our lepton pair
      std::vector<reco::TransientTrack> fitTracks;
      fitTracks.push_back(leptonTrack(part1));
      fitTracks.push_back(leptonTrack(part2));
      bool validTracks = (fitTracks[0].isValid() && fitTracks[1].isValid()); 
      cuts.applyCut("validTracks",validTracks);
      if (!cutStudyMode_ && !cuts.passAll()) continue;

      double decayLength=-999;
      double decayLengthSignificance=-999;
      double vertexChi2=9999;
      double dileptonMass = (part1.p4()+part2.p4()).mass();
      double corrDileptonMass=-999;
      double scaleCorrMass=-999;
      double caloCorrMass=-999;
      double triggerCorrDileptonMass=-999;
      double dPhi=-999;
      double dPhicorr=-999;
      double dPhitriggerCorr=-999;
      double cosThetaStar=-999;
      double eta_cand=-999;
      double phi_cand=-999;
      bool validVertex=false;
      int maxHitsBeforeVertex=-999;
      int maxHitsMissedAfterVertex=-999;
      if (validTracks) {
	TransientVertex leptonVertex;
	try {
	  leptonVertex = vertexFitter_->vertex(fitTracks);
	  validVertex=leptonVertex.isValid();
	} catch(cms::Exception& e) {};

	if (validVertex) {

	  vertexChi2=leptonVertex.normalisedChiSquared();

	  // calculate p_t corrected invariant mass
	  ROOT::Math::PxPyPzMVector totalMom(0,0,0,0);
	  ROOT::Math::PxPyPzMVector triggerCorrectedMom(0,0,0,0);
	  ROOT::Math::PxPyPzMVector caloCorrectedMom(0,0,0,0);
	  const std::vector<reco::TransientTrack> tracks = leptonVertex.refittedTracks();
	  double leptonMass=0.000511;
	  if (lepType==__lepTypeMuon) leptonMass=0.106;
	  if (lepType==__lepTypeTau) leptonMass=1.777;
	  for (unsigned i=0; i<tracks.size(); i++) {
	    totalMom+=ROOT::Math::PxPyPzMVector(tracks[i].track().px(),tracks[i].track().py(),
						tracks[i].track().pz(),leptonMass);
	    if (i==0 && leptonTriggerMatch_[ipart1]!=0) {
	      triggerCorrectedMom+=*(leptonTriggerMatch_[ipart1]);
	    } else if (i==1 && leptonTriggerMatch_[ipart2]!=0) {
	      triggerCorrectedMom+=*(leptonTriggerMatch_[ipart2]);
	    } else {
	      triggerCorrectedMom+=reco::Particle::LorentzVector(
				      	 tracks[i].track().px(),
				      	 tracks[i].track().py(),
				      	 tracks[i].track().pz(),
					 sqrt(tracks[i].track().p()*tracks[i].track().p()+leptonMass*leptonMass));
	    }
	    if (leptonCaloMatch_.size()>ipart2) {
	      if (i==0 && leptonCaloMatch_[ipart1]!=0) {
		caloCorrectedMom+=(leptonCaloMatch_[ipart1]/part1.p())*part1.p4();
	      } else if (i==1 && leptonCaloMatch_[ipart2]!=0) {
		caloCorrectedMom+=(leptonCaloMatch_[ipart2]/part2.p())*part2.p4();
	      } else {
		caloCorrectedMom+=reco::Particle::LorentzVector(
								tracks[i].track().px(),
								tracks[i].track().py(),
								tracks[i].track().pz(),
								sqrt(tracks[i].track().p()*tracks[i].track().p()+leptonMass*leptonMass));
	      }
	    }
	  }
	  corrDileptonMass = totalMom.M();
	  triggerCorrDileptonMass = triggerCorrectedMom.M();
	  caloCorrMass = caloCorrectedMom.M();

	  // decay length and significance
	  VertexDistanceXY vertTool;
	  Measurement1D vertexDist = vertTool.distance(leptonVertex.vertexState(),*primaryVertex_);
	  decayLength = vertexDist.value();
	  decayLengthSignificance = vertexDist.significance();

	  // vertex flight direction vs dilepton momentum direction
	  GlobalVector vertexDir = leptonVertex.position() - primaryVertex_->position();  
	  math::XYZTLorentzVector diLeptonMom = part1.p4()+part2.p4();
	  dPhi= fabs(deltaPhi<GlobalVector, math::XYZTLorentzVector>(vertexDir, diLeptonMom));
	  dPhicorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, totalMom));
	  dPhitriggerCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, triggerCorrectedMom));

	  // general kinematic information about candidate
	  eta_cand=diLeptonMom.eta();
	  phi_cand=diLeptonMom.phi();

	  // scale momenta to force dR=0 (phi component only).
	  // we scale up because we assume
	  // the momentum is underestimated due to energy loss
	  double momScale=-(mom1.y()*vertexDir.x()-mom1.x()*vertexDir.y())
	    /(mom2.y()*vertexDir.x()-mom2.x()*vertexDir.y());
	  dilepHistos_->fill("momScale",momScale,1,"momentum scaling factor for dR==0",100,-5,5);
	  if (momScale>1) {
	    scaleCorrMass=(part1.p4()+momScale*part2.p4()).mass();
	  } else {
	    scaleCorrMass=(part1.p4()/momScale+part2.p4()).mass();
	  }
	    

	  // angle between positive lepton momentum in dilepton rest frame and dilepton momentum
	  math::XYZTLorentzVector positiveMom;
	  if (part1.charge()>0) positiveMom = part1.p4(); else positiveMom = part2.p4();
	  ROOT::Math::Boost boost;
	  boost.SetComponents(-diLeptonMom/diLeptonMom.E());
	  math::XYZTLorentzVector cmsMom = boost(positiveMom);
          cosThetaStar = TMath::Cos(angle<math::XYZTLorentzVector,math::XYZTLorentzVector>
				    (diLeptonMom,cmsMom));

	  // check hit patterns. each track should have either all hits at radii larger than
	  // the vertex (when going outwards) or all hits at radii smaller than the vertex,
	  // unless it happens to go inside for a bit and then back out, in which case it would
          // have more than one hit on at least one later. (is that even possible with the CMS
	  // track reconstruction?)
	  // first simple attempt: require tracks to be outside of vertex only
	  for (unsigned i=0; i<fitTracks.size(); i++) {
	    std::pair<unsigned int, unsigned int> hitPat
	      = checkHitPattern_->analyze(fitTracks[i].track(),leptonVertex);
	    // first number: number of hits before vertex
	    // second number: number of missed layers after vertex
	    if (int(hitPat.first)>maxHitsBeforeVertex) maxHitsBeforeVertex=hitPat.first;
	    if (int(hitPat.second)>maxHitsMissedAfterVertex) maxHitsMissedAfterVertex=hitPat.second;
	  }

	  // maybe we can use some more sophisticated criterion: reject track if it is compatible
	  // with originating from IP and there is actually at least one hit between vertex and
	  // IP. Then things are simpler. I could restrict the hit check to pixel hits only.
	  // (if the track were coming from the IP, it would be very unlikely not to have pixel
	  // hits). try e.g. fitTracks[0].trajectoryStateClosestToPoint(primaryVertex_)
	}
      }

      cuts.applyCut("validVertex",validVertex);
      cuts.applyCut("vertexChi2",vertexChi2,-999,vertexChi2Cut_);
      cuts.applyLifetimeCut("dPhicorr",dPhicorr,-99,deltaPhiCut_);
      cuts.applyCut("maxHitsBeforeVertex",maxHitsBeforeVertex,-999,maxHitsBeforeVertex_);
      cuts.applyLifetimeCut("decayLength2D",decayLength,
			    decayLengthCut_,999999);
      cuts.applyLifetimeCut("decayLengthSignificance2D",
			    decayLengthSignificance,
			    decayLengthSignificanceCut_,999999);
      
      // pseudo-cuts that shouldn't do anything other than create diagnostic plots
      cuts.applyCut("cosThetaStar",cosThetaStar,-2,2);
      cuts.applyCut("eta_cand",eta_cand,-999,999);
      cuts.applyCut("phi_cand",phi_cand,-999,999);
      cuts.applyCut("mass_corr",corrDileptonMass,-9999,9999);
      cuts.applyCut("mass_triggercorr",triggerCorrDileptonMass,-9999,9999);
      cuts.applyCut("mass_scalecorr",scaleCorrMass,-9999,9999);
      cuts.applyCut("mass_calocorr",caloCorrMass,-9999,9999);
      cuts.applyCut("dPhi",dPhi,-999,999);
      cuts.applyCut("dPhitriggerCorr",dPhitriggerCorr,-999,999);
      cuts.applyCut("maxHitsMissedAfterVertex",maxHitsMissedAfterVertex,-999,999);
      cuts.applyCut("numPrimaryVertices",numPV_,0,999);

      int numStandAloneMuons=0;
      if (part1.isStandAloneMuon()) ++numStandAloneMuons;
      if (part2.isStandAloneMuon()) ++numStandAloneMuons;
      cuts.applyCut("numStandAloneMuons",numStandAloneMuons,-999,maxNumStandAloneMuons_);
      if (leptonCaloMatch_.size()>ipart2) {
	int nCaloMatches=0;
	if (leptonCaloMatch_[ipart1]!=0) ++nCaloMatches;
	if (leptonCaloMatch_[ipart2]!=0 && leptonCaloMatch_[ipart2]!=leptonCaloMatch_[ipart1]) ++nCaloMatches;
	cuts.applyCut("numCaloMatches",nCaloMatches,minNumCaloMatches_,999);
      }

      // do generator level analysis of our dilepton candidate to come up with a candidate category
      categoryType cat;
      int pdgIdCommonOrigin=0;
      int pdgSignalDecayLepton=0;
      if (ancestor[ipart1]!=NULL && ancestor[ipart1]==ancestor[ipart2]) {
	// fully reconstructed signal!
	pdgSignalDecayLepton=decayChannel(*(ancestor[ipart1]));
	if (pdgSignalDecayLepton==abs(leptonPDGId_)) {
	  // this is the proper decay channel too
	  cat=__catSignal;
	} else {
          // leakage from other signal decay channel (e.g. tautau reconstructed as ee)
	  cat=__catSignalWrongChannel;
	}
      } else if (ancestor[ipart1]!=NULL && ancestor[ipart2]!=NULL) {
	// mixture of contributions from different signal particles
	cat=__catSignalMixed;
      } else if (ancestor[ipart1]!=NULL) {
	// one lepton from signal, one from other source
	if (genParticle[ipart2]) {
	  cat=__catOneSigOneBkgWithGen;
	} else {
	  cat=__catOneSigOneBadReco;
	}
      } else if (ancestor[ipart2]!=NULL) {
	// one lepton from signal, one from other source
	if (genParticle[ipart1]) {
	  cat=__catOneSigOneBkgWithGen;
	} else {
	  cat=__catOneSigOneBadReco;
	}
      } else {
	// neither lepton originates from signal
	// check whether there is any other common origin (e.g. a Z)
	const reco::Candidate* origin=commonOrigin(genParticle[ipart1],genParticle[ipart2]);
	if (origin) {
	  // common (non-signal) origin found
	  cat=__catGoodButWrongType;
	  pdgIdCommonOrigin=abs(origin->pdgId());
	} else {
	  cat=__catBackground;
	}
      }

      // choose mass variable to show along with cut variable plots
      double plot_mass=dileptonMass; // good for electron and muon channels
      if (lepType==__lepTypeTrack) {
	if (leptonName_=="etrack") {
	  plot_mass=caloCorrMass;
	} else if (leptonName_=="mutrack") {
	  plot_mass=corrDileptonMass;
	} else {
	  throw cms::Exception("InvalidTrackChannel") << "unknown track channel for plot_mass";
	}
      }

      // cut efficiency bookkeeping
      for (unsigned itrig=0; itrig<triggers_.size(); itrig++) {
	if (cat==__catSignal) {
	  cutSummarySignal_[triggers_[itrig]]->addEntry(cuts,weightMap_,
							numDecays,
							ctau1,
							ctau2,
							leptonD0[ipart1],
							leptonD0[ipart2],
							plot_mass);
	} else if (cat==__catBackground || cat==__catGoodButWrongType) {
	  cutSummaryBackground_[triggers_[itrig]]->addEntry(cuts,weightMap_,
							    numDecays,
							    ctau1,
							    ctau2,
							    leptonD0[ipart1],
							    leptonD0[ipart2],
							    plot_mass);
	} else {
	  cutSummaryPartial_[triggers_[itrig]]->addEntry(cuts,weightMap_,
							 numDecays,
							 ctau1,
							 ctau2,
							 leptonD0[ipart1],
							 leptonD0[ipart2],
							 plot_mass);
	}
      }


      // skip all the rest (filling histograms!) for candidates failing the cuts
      if (!cuts.passAll()) continue;


      ++ncand;

      // print candidate information
      std::cout << "found candidate, run " << iEvent.id().run() <<", lumi block "
		<< iEvent.luminosityBlock() << ", event " << iEvent.id().event()
		<< ", mass" << corrDileptonMass << " (" << leptonName_ << ")"
		<< std::endl;
      std::cout << "   found candidate has constituent momenta " <<
	part1.momentum() << ", " << part2.momentum() << std::endl;

      dilepHistos_->fill("mass_all",plot_mass,1,
			 "reconstructed di"+leptonName_+" mass",100,0,histMaxM_);
      if (leptonTriggerMatch_[ipart1]!=0 || leptonTriggerMatch_[ipart2]!=0) {
	dilepHistos_->fill("mass_all_matched",plot_mass,1,
			   "reconstructed di"+leptonName_+
			   " mass, >0 trigger matched leptons",100,0,histMaxM_);
      }

      int numTrueLeptons=0;
      if (genParticle[ipart1]) {
	if (abs(genParticle[ipart1]->pdgId())==abs(leptonPDGId_)) ++numTrueLeptons;
      }
      if (genParticle[ipart2]) {
	if (abs(genParticle[ipart2]->pdgId())==abs(leptonPDGId_)) ++numTrueLeptons;
      }
      dilepHistos_->fill("number_of_actual_leptons",numTrueLeptons,1,
			 "number of true leptons in di"+leptonName_+" candidate",
			 3,0,3);

      // do some MC analysis of our dilepton candidate
      if (cat==__catSignal) {
	// fully reconstructed signal!
	dilepHistos_->fill("mass_signal",dileptonMass,1,
			   "reconstructed di"+leptonName_+" mass",100,0,histMaxM_);
	GenParticleInfo genInfo(ancestor[ipart1]);
	dilepHistos_->fill("trueDecayLength2D",genInfo.flightLength2D(),1,
			   "true decay length, properly reconstructed di"+leptonName_+"s",
			   histNBinsRt_,0,histMaxRt_);
	if (numDecays==1) {
	  dilepHistos_->fill("trueDecayLength2D_1",genInfo.flightLength2D(),1,
			     "true decay length, properly reconstructed di"+leptonName_+"s",
			     histNBinsRt_,0,histMaxRt_);
	} else if (numDecays==2) {
	  dilepHistos_->fill("trueDecayLength2D_2",genInfo.flightLength2D(),1,
			     "true decay length, properly reconstructed di"+leptonName_+"s",
			     histNBinsRt_,0,histMaxRt_);
	}
	dilepHistos_->fill("di"+leptonName_+"_r_resolution",
			   genInfo.flightLength3D(),decayLength,1,
			   "true vs reconstructed 3d decay length, properly reconstructed di"
			   +leptonName_+"s",histNBinsRt_,0,histMaxRt_,histNBinsRt_,0,histMaxRt_);
	dilepHistos_->fill("numHitsBeforeVertex_signal",double(maxHitsBeforeVertex),1,
			   "number of hits inside of presumed vertex, signal only",10,0,10);
	dilepHistos_->fill("numHitsMissedAfterVertex_signal",double(maxHitsMissedAfterVertex),1,
			   "number of missed layers after presumed vertex, signal only",10,0,10);

	// study of tracker isolation as function of pile-up
	dilepHistos_->fill("isolation_vs_pileup_signal",
			   isol1,genProp.aveNumPileup3BX(),1,
			   "tracker isolation vs number of pile-up events",
			   100,0.,100.,25,0,25);
	dilepHistos_->fill("isolation_vs_pileup_signal",
			   isol2,genProp.aveNumPileup3BX(),1,
			   "tracker isolation vs number of pile-up events",
			   100,0.,100.,25,0,25);
			   

	// check trigger matching of properly reconstructed signal events
	if (leptonTriggerMatch_[ipart1]!=0 || leptonTriggerMatch_[ipart2]!=0) {
	  dilepHistos_->fill("mass_signal_matched",plot_mass,1,"reconstructed di"
			     +leptonName_+" mass, >0 trigger matched leptons",
			     100,0,histMaxM_);
	  dilepHistos_->fill("trueDecayLength2D_matched",genInfo.flightLength2D(),1,
			     "true decay length, properly reconstructed di"
			     +leptonName_+"s, >0 trigger matched leptons",
			     histNBinsRt_,0,histMaxRt_);
	}

      } else if (cat==__catSignalWrongChannel) {
	// wrong decay mode, but otherwise properly reconstructed signal
	dilepHistos_->fill("mass_signal_leakage",plot_mass,1,"reconstructed di"
			   +leptonName_+" mass",100,0,histMaxM_);
	dilepHistos_->fill("decayChannel_leakage",pdgSignalDecayLepton,1,"decay channels leaking "
			   "into "+leptonName_+" selection",5,11,16);
      } else if (cat==__catSignalMixed) {
	// mixture of contributions from different signal particles
	dilepHistos_->fill("mass_signal_mixed",plot_mass,1,"reconstructed di"
			   +leptonName_+" mass",100,0,histMaxM_);
      } else if (cat==__catOneSigOneBkgWithGen || cat==__catOneSigOneBadReco) {
	// one signal particle, one non-signal particle
	dilepHistos_->fill("mass_onesignal_oneother",plot_mass,1,
			   "reconstructed di"+leptonName_+" mass",100,0,histMaxM_);
      } else if (cat==__catGoodButWrongType) {
	// common origin found, but not a signal particle
	std::stringstream name,title;
	name << "mass_" << pdgIdCommonOrigin << "_same";
	title << "di"+leptonName_+" mass, same mother particle of PDG ID " << pdgIdCommonOrigin;
	dilepHistos_->fill(name.str(),plot_mass,1,title.str(),100,0,histMaxM_);
	dilepHistos_->fill("mass_other_resonance",plot_mass,1,
			   "mass distribution, same non-signal mother particle",
			   100,0,histMaxM_);
      } else {
	// no common origin whatsoever
	std::stringstream name,title;
	name << "mass_background";
	title << "di"+leptonName_+" mass, no common origin";
	dilepHistos_->fill(name.str(),plot_mass,1,title.str(),100,0,histMaxM_);
      }
    }
  }

  dilepHistos_->fill("dilepton_num",ncand,1,
		     "number of di"+leptonName_+" candidates in event",10,0,10);
  
}


template<class Lepton> 
void DisplacedDileptonAnalysis::leptonPlots(const Lepton* cand,
					    const reco::GenParticle* genPart,
					    const leptonType lepType,
					    const bool passesAllCuts,
					    const bool isSignal) {

  // fill a few particle-level histograms
  if (genPart) {
    lepHistos_->fill("pt_resolution",genPart->pt(),cand->pt(),1,
		     "reconstructed vs generated "+leptonName_+" p_{t}",
		     100,0,100,100,0,100);
    lepHistos_->fill("charge_identification",
		     genPart->charge(),cand->charge(),1,
		     "reconstructed vs generated "+leptonName_+" charge",
		     3,-1.1,1.1,3,-1.1,1.1);
    lepHistos_->fill("trueProdVtxRad",genPart->vertex().R(),1,
		     "3d radius of MC production vertex of "+leptonName_,
		     100,0,histMaxRt_);
    lepHistos_->fill("pt_vs_radius",genPart->vertex().R(),cand->pt(),1,
		     "reconstructed p_t vs generated 3d radius, "+leptonName_+"s",
		     100,0,histMaxRt_,100,0,100);
    if (leptonName_=="muon") {
      if (!cand->isGlobalMuon() && !cand->isTrackerMuon()) {
	lepHistos_->fill("pt_vs_radius_standalone",
			 genPart->vertex().R(),cand->pt(),1,
			 "reconstructed p_t vs generated 3d radius, "
			 +leptonName_+"s",100,0,histMaxRt_,100,0,100);
      } else {
	lepHistos_->fill("pt_vs_radius_allothers",
			 genPart->vertex().R(),cand->pt(),1,
			 "reconstructed p_t vs generated 3d radius, "
			 +leptonName_+"s",100,0,histMaxRt_,100,0,100);
      }
    }
  }

  if (passesAllCuts) {
    lepHistos_->fill("leptonsVsRun",1,1,"number of leptons passing all cuts, per run",
		     runMax_-runMin_+1,runMin_,runMax_);
  }

  int ngen=0;
  if (genPart) ngen=1;
  lepHistos_->fill("num_genparticles",ngen,1,
		   "number of generator particles associated to "
		   +leptonName_+" candidate",5,0,5);
}


const reco::Candidate* DisplacedDileptonAnalysis::signalOrigin(const reco::Candidate* part) {

  const reco::Candidate* cand = part;
  const reco::Candidate* found2 = NULL;
  const reco::Candidate* found3 = NULL;
  while (cand!=0 && !found3) {
    if (abs(cand->pdgId())==signalPDGId_) {
      if (cand->status()==3) found3=cand; else found2=cand;
    }
    cand=cand->mother();
  }
  if (found3) return found3; else return found2;
}


int DisplacedDileptonAnalysis::decayChannel(const reco::Candidate& part) {

  // first check direct descendants
  for (unsigned i=0; i<part.numberOfDaughters(); i++) {
    const reco::Candidate* daughter = part.daughter(i);
    int pid = abs(daughter->pdgId());
    if (pid==11 || pid==13 || pid==15) return pid;
  }
  return 0;
}


const reco::Candidate* DisplacedDileptonAnalysis::commonOrigin(const reco::Candidate* part1,
							       const reco::Candidate* part2) {

  if (!part1 || !part2) return 0;

  const reco::Candidate* cand1=part1;
  const reco::Candidate* cand2=part2;
  const reco::Candidate* common=NULL;
  while (cand1!=NULL && common==NULL) {
    cand2=part2;
    while (cand2!=NULL && common==NULL) {
      if (cand1==cand2) common=cand1;
      cand2=cand2->mother();
    }
    cand1=cand1->mother();
  }
  if (common) {
    // if the first common particle is a gluon or quark in the hard interaction,
    // then we consider this no common origin
    int pdgId = abs(common->pdgId());
    if (common->status()==3 && (pdgId<6 || pdgId==21)) common=NULL;
  }
  return common;
}


void DisplacedDileptonAnalysis::beginJob()
{
  vertexFitter_ = new KalmanVertexFitter(true);
}


void DisplacedDileptonAnalysis::endJob() 
{

  // cut flow histograms
  std::string name=leptonName_;

  // track reconstruction efficiency
  lepHistos_->makeEfficiencyHistogram("trueLeptonRadWithTrack",
				      mcHistos_->getHist1D("leptonProdVtxRadius3D"),
				      "track_reco_efficiency_decaylength3D");
  // lepton reconstruction efficiency
  lepHistos_->makeEfficiencyHistogram("trueProdVtxRad",
				      mcHistos_->getHist1D("leptonProdVtxRadius3D"),
				      "lepton_reco_efficiency_decaylength3D");

  // dilepton reconstruction efficiency histogram
  dilepHistos_->makeEfficiencyHistogram("trueDecayLength2D",
					mcHistos_->getHist1D("trueDecayLength2D"),
					"dilepton_reco_efficiency_decayLength2D");
  // efficiency to have at least one lepton trigger matched
  dilepHistos_->makeEfficiencyHistogram("trueDecayLength2D_matched",
					"trueDecayLength2D",
					"dilepton_match_efficiency_decayLength2D");
  // trigger matching and reconstruction efficiency combined
  dilepHistos_->makeEfficiencyHistogram("trueDecayLength2D_matched",
					mcHistos_->getHist1D("trueDecayLength2D"),
					"dilepton_trigreco_efficiency_decayLength2D");

  // print trigger summary information. use a stupid simple sort algorithm
  std::cout << "===========================================================" << std::endl;
  std::cout << "=== TRIGGER SUMMARY =======================================" << std::endl;
  std::cout << "===========================================================" << std::endl;
  unsigned highestRate=numProcessedEvents_;
  while(highestRate>0) {
    unsigned nextLowest=0;
    for (std::map<std::string,unsigned int>::iterator it=trigPathSummary_.begin();
	 it!=trigPathSummary_.end(); it++) {
      if (it->second==highestRate) {
	std::cout << "PATH " << std::setw(6) << it->second << "  " << it->first << std::endl;
      } else if (it->second<highestRate && it->second>nextLowest) {
	nextLowest=it->second;
      }
    }
    highestRate=nextLowest;
  }
  std::cout << "===========================================================" << std::endl;
  highestRate=numProcessedEvents_;
  while(highestRate>0) {
    unsigned nextLowest=0;
    for (std::map<std::string,unsigned int>::iterator it=trigFilterSummary_.begin();
	 it!=trigFilterSummary_.end(); it++) {
      if (it->second==highestRate) {
	std::cout << "FILTER " << std::setw(6) << it->second << "  " << it->first << std::endl;
      } else if (it->second<highestRate && it->second>nextLowest) {
	nextLowest=it->second;
      }
    }
    highestRate=nextLowest;
  }
  std::cout << "===========================================================" << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DisplacedDileptonAnalysis);
