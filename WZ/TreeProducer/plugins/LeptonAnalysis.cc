#include "TreeProducer/TreeProducer/plugins/LeptonAnalysis.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
// #include "DataFormats/Common/interface/Vector.h"
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
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include <sstream>
#include <iomanip>

// PAT
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

void LeptonAnalysis::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & setup)
{
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  numProcessedEvents_ += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
  lumi.getByLabel("nEventsPostHLTFilter", nEventsFilteredCounter);
  numEventsPassingTrigger_ += nEventsFilteredCounter->value;
}

LeptonAnalysis::LeptonAnalysis(const edm::ParameterSet& iConfig):
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
  leptonChargeCut_(iConfig.getParameter<bool>("leptonChargeCut")),
  useMCTruth_(iConfig.getParameter<bool>("UseMCTruth")),
  numProcessedEvents_(0),
  numEventsPassingTrigger_(0)
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

  // introduce short dCache read timeout to work around problems at RAL
  TDCacheFile::SetOpenTimeout(1800);

  // Create the tree
  edm::Service<TFileService> fs_;
  outputTree_ = fs_->make<TTree>("outputTree","outputTree");
  // Single lepton variables
  outputTree_->Branch("candidates","Candidates",&candidates_,32000,0);
  outputTree_->Branch("triggers","std::vector<std::string>",&triggers_,32000,0);

  // To save the total number of processed events
  hTotalProcessedEvents_ = new TH1F("totalProcessedEvents", "Total processed events", 1, 0, 1);
  hEventsPassingTrigger_ = new TH1F("eventsPassingTrigger", "Events passing trigger selection", 1, 0, 1);
}


LeptonAnalysis::~LeptonAnalysis()
{
}


void LeptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // ++numProcessedEvents_;

  leptons_.clear();

  // std::cout << "number of processed events = " << numProcessedEvents_ << std::endl;

  // event setup
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);
  if (checkHitPattern_==NULL) checkHitPattern_ = new CheckHitPattern(iSetup);

  // look at trigger information, this fills the vars_.triggers vector of strings
  bool passesTrigger = doTrigger(iEvent);

  if( passesTrigger ) {
    // ++numEventsPassingTrigger_;

    // get beamspot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel ("offlineBeamSpot", beamSpotHandle);
    beamSpot_=(*beamSpotHandle);

    edm::Handle<std::vector<pat::MET> > metHandle;
    iEvent.getByLabel ("patMETs", metHandle);
    MET_ = metHandle->begin()->pt();
    METPhi_ = metHandle->begin()->phi();

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

    // evaluate generator-level event properties needed for reweighting,
    // such as PDF, pile-up vertices,
    // decay channels and lifetimes of long-lived particles etc
    if( useMCTruth_ ) {
      GenEventProperties genProp(iEvent,signalPDGId_,generatorTag_,
                                 pileupTag_,genEventInfoTag_);
      nvtx_m1_ = genProp.numPileup(-1);
      nvtx_0_ = genProp.numPileup(0);
      nvtx_p1_ = genProp.numPileup(1);
    }

    // dilepton reconstruction
    // if (thisLepton_==__lepTypeElectron) {
    //   doDecayChannel<pat::Electron>(iEvent,__lepTypeElectron);
    // } else if (thisLepton_==__lepTypeMuon) {
    //   doDecayChannel<pat::Muon>(iEvent,__lepTypeMuon);
    // } else if (thisLepton_==__lepTypeTau) {
    //   doDecayChannel<pat::Tau>(iEvent,__lepTypeTau);
    // } else {
      // doDecayChannel<PseudoLepton>(iEvent,__lepTypeTrack,genProp);
      doDecayChannel<PseudoLepton>(iEvent,__lepTypeTrack);
      // }
  }
}

template<class Lepton>
const reco::Particle::LorentzVector * LeptonAnalysis::findTrigMatch(const Lepton & lepton,
								    const edm::Handle< pat::TriggerEvent > & triggerEvent)
{
  // find deltaR match among trigger objects
  double lep_eta = lepton.eta();
  double lep_phi = lepton.phi();
  double min_deltaR=9999.;
  const reco::Particle::LorentzVector* trig_mom=0;

  // loop over all triggers we are interested in
  std::cout << "hltPaths_.size() = " << hltPaths_.size() << std::endl;
  for (unsigned itrig=0; itrig<hltPaths_.size(); itrig++) {
    // std::cout << "Trigger path number " << itrig << std::endl;
    if (!triggerEvent->path(hltPaths_[itrig])) continue;
    if (!triggerEvent->path(hltPaths_[itrig])->wasAccept()) continue;
    // get all trigger objects associated with this HLT path
    const pat::TriggerObjectRefVector triggerObjects(triggerEvent->pathObjects(hltPaths_[itrig],false));
    std::cout << "triggerObjects.size() = " << triggerObjects.size() << std::endl;
    for (unsigned i=0; i<triggerObjects.size(); i++) {
      // std::cout << "Trigger object number " << i << std::endl;
      double trig_eta=triggerObjects[i]->eta();
      double trig_phi=triggerObjects[i]->phi();
      double dR=deltaR(lep_eta,lep_phi,trig_eta,trig_phi);
      if (dR<min_deltaR) {
        min_deltaR=dR;
        trig_mom=&(triggerObjects[i]->p4());
      }
    }
  }

  bool trigMatch=(min_deltaR<0.1);
  if (trigMatch) {
    return trig_mom;
  }
  return 0;
}

// Given a PseudoLepton and the barrel and endcap superclusters returns a pointer to a maching supercluster or a null pointer if no match was found.
template<class Lepton>
const reco::SuperCluster * LeptonAnalysis::findCaloMatch(const Lepton & lepton,
							 const edm::Handle< std::vector<reco::SuperCluster> > & barrelSuperClusters,
							 const edm::Handle< std::vector<reco::SuperCluster> > & endcapSuperClusters)
{
  // find deltaR match among superclusters
  double lep_eta=lepton.eta();
  double lep_phi=lepton.phi();
  double min_deltaR=9999.;
  // double calo_energy=0;
  const reco::SuperCluster * matchingSC = 0;

  const double deltaRcut=0.2;

  // loop over all barrel superclusters
  std::cout << barrelSuperClusters->size() << std::endl;
  if (fabs(lep_eta)<1.48+1.5*deltaRcut) {
    for (unsigned iclus=0; iclus<barrelSuperClusters->size(); iclus++) {
      double clus_eta=(*barrelSuperClusters)[iclus].eta();
      double clus_phi=(*barrelSuperClusters)[iclus].phi();
      double dR=deltaR(lep_eta,lep_phi,clus_eta,clus_phi);
      if (dR<min_deltaR) {
        min_deltaR=dR;
        // calo_energy=(*barrelSuperClusters)[iclus].energy();
        matchingSC = &(*barrelSuperClusters)[iclus];
      }
    }
  }
  // loop over all endcap superclusters
  std::cout << endcapSuperClusters->size() << std::endl;
  if (fabs(lep_eta)>1.48-1.5*deltaRcut) {
    for (unsigned iclus=0; iclus<endcapSuperClusters->size(); iclus++) {
      double clus_eta=(*endcapSuperClusters)[iclus].eta();
      double clus_phi=(*endcapSuperClusters)[iclus].phi();
      double dR=deltaR(lep_eta,lep_phi,clus_eta,clus_phi);
      if (dR<min_deltaR) {
        min_deltaR=dR;
        // calo_energy=(*endcapSuperClusters)[iclus].energy();
        matchingSC = &(*endcapSuperClusters)[iclus];
      }
    }
  }

  bool caloMatch = (min_deltaR<deltaRcut);
  if( caloMatch ) return matchingSC;
  return 0;
}

/// Check if the triggers selected via cfg are firing for this event.
bool LeptonAnalysis::doTrigger(const edm::Event& iEvent)
{
  // PAT trigger information
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEvent_, triggerEvent );
  if (triggerEvent.failedToGet()) {
    std::cout << "WARNING: no TriggerEvent found" << std::endl;
    return false;
  }

  bool takeAllPaths = (hltPaths_.size() == 1 && hltPaths_[0] == "*");
  // now loop over all paths etc and fill histogram for each one that was found
  triggers_.clear();
  for (unsigned i=0; i<triggerEvent->acceptedPaths().size(); i++) {
    std::string trigname(triggerEvent->acceptedPaths()[i]->name());
    ++trigPathSummary_[trigname];
    for (unsigned itrg=0; itrg<hltPaths_.size(); itrg++) {
      if( takeAllPaths || trigname==hltPaths_[itrg] ) {
        // FIXME: might want to select a reduced number of triggers to save space in the tree.
        triggers_.push_back(trigname);
      }
    }
  }
  return (triggers_.size()>0);
}

/// Return the closest stable genParticle in deltaR
const reco::GenParticle* LeptonAnalysis::getGenParticleMatch(const edm::Handle<edm::View<reco::GenParticle> > & mcParticles,
							     const reco::Particle::LorentzVector& mom)
{
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

/// Compute the sumPt in a cone. This is absolute isolation.
double LeptonAnalysis::trackerIsolation(const edm::Event& iEvent,
					const reco::Track & plept,
					const reco::Track & pveto)
{
  // FIXME: move to relative isolation
  // FIXME: move the Handle outside and pass it as an argument

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

/// Convert to a transient track and assign the beamspot.
reco::TransientTrack LeptonAnalysis::leptonTrack(const PseudoLepton &track)
{
  reco::TransientTrack ttrack=trackBuilder_->build(track);
  ttrack.setBeamSpot(beamSpot_);
  return ttrack;
}

/// Convert to a transient track and assign the beamspot.
reco::TransientTrack LeptonAnalysis::leptonTrack(const pat::Electron &electron)
{
  if (!electron.gsfTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(electron.gsfTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  }
  // FIXME
  return reco::TransientTrack();
}

/// Convert to a transient track and assign the beamspot. Take innerTrack if available or outerTrack otherwise.
reco::TransientTrack LeptonAnalysis::leptonTrack(const pat::Muon &muon)
{
  if (!muon.innerTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(muon.innerTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  } else if (!muon.outerTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(muon.outerTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  }
  // FIXME
  return reco::TransientTrack();
}

/// Convert to a transient track and assign the beamspot.
reco::TransientTrack LeptonAnalysis::leptonTrack(const pat::Tau &tau)
{
  if (!tau.leadTrack().isNull()) {
    reco::TransientTrack ttrack=trackBuilder_->build(tau.leadTrack());
    ttrack.setBeamSpot(beamSpot_);
    return ttrack;
  }
  // FIXME
  return reco::TransientTrack();
}

/// lepton IDs are specific to lepton type, but we provide a default template here
template<class Lepton>
bool LeptonAnalysis::leptonID(const Lepton& lepton)
{
  return true;
}

/// override template for the electron channel, where we do use a specific ID requirement
bool LeptonAnalysis::leptonID(const pat::Electron& electron)
{
  return int(electron.electronID("eidLoose"))&1;
}

/// Default method, returns -999
double LeptonAnalysis::leptonTiming(const PseudoLepton &track)
{
  return -999;
}

/// Not implemented, returns -999
double LeptonAnalysis::leptonTiming(const pat::Electron &electron)
{
  //DetId seedDetID = electron->superCluster()->seed()->seed();
  return -999;
}

/// Muon timing. Assumes muon coming from inside the detector and has beta == 1
double LeptonAnalysis::leptonTiming(const pat::Muon &muon)
{
  if (muon.isTimeValid()) {
    // use time under the assumption the muon came from inside of the detector
    // and had beta==1.
    return muon.time().timeAtIpInOut;
  } else {
    return -999;
  }
}

/// Not implemented, returns -999
double LeptonAnalysis::leptonTiming(const pat::Tau &tau)
{
  // not implemented
  return -999;
}


template<class Lepton>
void LeptonAnalysis::doDecayChannel(const edm::Event& iEvent,
				    const leptonType lepType)
{
  edm::Handle<edm::View<Lepton> > particles;
  iEvent.getByLabel(leptonCollTag_[lepType],particles);

  // Load superCluster collections if this is for electrons
  edm::Handle< std::vector<reco::SuperCluster> > barrelSuperClusters;
  edm::Handle< std::vector<reco::SuperCluster> > endcapSuperClusters;
  if( leptonName_ == "etrack" ) {
    iEvent.getByLabel( barrelSuperClusters_, barrelSuperClusters );
    if (barrelSuperClusters.failedToGet()) {
      std::cout << "WARNING: cannot access barrel superclusters for matching" << std::endl;
    }
    iEvent.getByLabel( endcapSuperClusters_, endcapSuperClusters );
    if (endcapSuperClusters.failedToGet()) {
      std::cout << "WARNING: cannot access endcap superclusters for matching" << std::endl;
    }
  }

  // Load gen particles or set the flag to false if they are not available
  bool MCTruthFound = true;
  edm::Handle<edm::View<reco::GenParticle> > mcParticles;
  if( useMCTruth_ ) {
    iEvent.getByLabel(generatorTag_,mcParticles);
    if(mcParticles.failedToGet()) {
      std::cout << "WARNING: cannot access MC truth information." << std::endl;
      MCTruthFound = false;
    }
  }

  // Load trigger objects
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEvent_, triggerEvent );
  if (triggerEvent.failedToGet()) {
    std::cout << "WARNING: cannot access triggerEvent for matching" << std::endl;
    //    return;
  }

  // First loop on leptons to fill lepton-based information (trigger matches, calo-matches, muon-matches, etc...)
  // ------------------------------------------------------------------------------------------------------------
  for( typename edm::View<Lepton>::const_iterator it = particles->begin(); it != particles->end(); ++it ) {
    // std::cout << "it->pt() = " << it->pt() << std::endl;
    if (it->pt() < leptonPtCut_) continue;
    // Build and save the transient track. This contains information on the beamspot and will be used later for extrapolations and vertex fitting.
    LeptonContainer lepton(leptonTrack(*it));
    lepton.p4 = &(it->p4());
    lepton.isStandAloneMuon = it->isStandAloneMuon();
    lepton.isGlobalMuon = it->isGlobalMuon();
    lepton.isTrackerMuon = it->isTrackerMuon();
    lepton.trackIso = it->trackIso();
    lepton.ecalIso = it->ecalIso();
    lepton.hcalIso = it->hcalIso();
    lepton.normChi2 = it->normChi2();
    lepton.numberOfValidTrackerHits = it->numberOfValidTrackerHits();
    lepton.numberOfValidPixelHits = it->numberOfValidPixelHits();
    lepton.numberOfValidMuonHits = it->numberOfValidMuonHits();
    lepton.numberOfMatchedStations = it->numberOfMatchedStations();
    // Calo matching if it is an electron track analysis
    if ( leptonName_ == "etrack" ) {

      lepton.matchedSC = findCaloMatch(*it, barrelSuperClusters, endcapSuperClusters);

      if ( lepton.matchedSC == 0 ) std::cout << "Did not find a match" << std::endl;
    }
    // Gen matching if this is MC
    if( useMCTruth_ && MCTruthFound ) {
      lepton.genPart = getGenParticleMatch(mcParticles,it->p4());
      if( lepton.genPart ) lepton.motherPart = signalOrigin(lepton.genPart);
    }
    // Trigger matching
    lepton.triggerMatch = findTrigMatch(*it, triggerEvent);

    lepton.vx = it->vx();
    lepton.vy = it->vy();
    lepton.vz = it->vz();

    leptons_.push_back(lepton);
  }

  // Initialize all variables
  initializeVars();


  // Event level variables
  candidates_.run = iEvent.run();
  candidates_.event = iEvent.id().event();
  candidates_.numPV = numPV_;
  if( useMCTruth_ ) {
    candidates_.nvtx_m1 = nvtx_m1_;
    candidates_.nvtx_0 = nvtx_0_;
    candidates_.nvtx_p1 = nvtx_p1_;
  }
  candidates_.MET = MET_;
  candidates_.METPhi = METPhi_;

  // Find candidates
  // ---------------
  // loop over the saved leptons
  for (unsigned ipart1=0; ipart1<leptons_.size(); ++ipart1) {
    const LeptonContainer part1 = leptons_[ipart1];

    // Save all the leptons
    candidates_.leptons_.push_back(TreeLepton());
    TreeLepton & lepton = candidates_.leptons_.back();
    reco::Track tempLepton = leptons_[ipart1].transientTrack.track();
    lepton.pt = tempLepton.pt();
    lepton.eta = tempLepton.eta();
    lepton.phi = tempLepton.phi();
    lepton.charge = tempLepton.charge();
    lepton.index = ipart1;
    lepton.vx = leptons_[ipart1].vx;
    lepton.vy = leptons_[ipart1].vy;
    lepton.vz = leptons_[ipart1].vz;
    lepton.trackIso = leptons_[ipart1].trackIso;
    lepton.ecalIso = leptons_[ipart1].ecalIso;
    lepton.hcalIso = leptons_[ipart1].hcalIso;
    lepton.normChi2 = leptons_[ipart1].normChi2;
    lepton.numberOfValidTrackerHits = leptons_[ipart1].numberOfValidTrackerHits;
    lepton.numberOfValidPixelHits = leptons_[ipart1].numberOfValidPixelHits;
    lepton.numberOfValidMuonHits = leptons_[ipart1].numberOfValidMuonHits;
    lepton.numberOfMatchedStations = leptons_[ipart1].numberOfMatchedStations;
    try {
      lepton.d0 = leptons_[ipart1].transientTrack.stateAtBeamLine().transverseImpactParameter().value();
      lepton.absD0Significance = fabs(leptons_[ipart1].transientTrack.stateAtBeamLine().transverseImpactParameter().significance());
    } catch (cms::Exception& e) {};
    lepton.iso = trackerIsolation(iEvent,leptons_[ipart1].transientTrack.track(),leptons_[ipart1].transientTrack.track());
    lepton.triggerMatch = 0;
    if( leptons_[ipart1].triggerMatch != 0 ) lepton.triggerMatch = 1;
    lepton.isStandAlone = leptons_[ipart1].isStandAloneMuon;
    lepton.isGlobalMuon = leptons_[ipart1].isGlobalMuon;
    lepton.isTrackerMuon = leptons_[ipart1].isTrackerMuon;
    lepton.hasCaloMatch = 0;
    if( leptons_[ipart1].matchedSC != 0 ) lepton.hasCaloMatch = 1;

    for (unsigned ipart2=ipart1+1; ipart2<leptons_.size(); ++ipart2) {
      const LeptonContainer part2 = leptons_[ipart2];

      candidates_.candidates_.push_back(TreeCandidate());
      TreeCandidate & candidate = candidates_.candidates_.back();

      unsigned iHighPt=ipart1;
      unsigned iLowPt=ipart2;
      if (part1.transientTrack.track().pt()<part2.transientTrack.track().pt()) {
        iHighPt=ipart2;
        iLowPt=ipart1;
      }
      reco::Track lowPtLepton = leptons_[iLowPt].transientTrack.track();
      reco::Track highPtLepton = leptons_[iHighPt].transientTrack.track();
      candidate.leptonPtL = lowPtLepton.pt();
      candidate.leptonPtH = highPtLepton.pt();
      candidate.leptonEtaL = lowPtLepton.eta();
      candidate.leptonEtaH = highPtLepton.eta();
      candidate.leptonChargeL = lowPtLepton.charge();
      candidate.leptonChargeH = highPtLepton.charge();
      candidate.leptonIndexL = iLowPt;
      candidate.leptonIndexH = iHighPt;

      try {
        candidate.leptonD0L = leptons_[iLowPt].transientTrack.stateAtBeamLine().transverseImpactParameter().value();
        candidate.leptonAbsD0SignificanceL = fabs(leptons_[iLowPt].transientTrack.stateAtBeamLine().transverseImpactParameter().significance());
        candidate.leptonD0H = leptons_[iHighPt].transientTrack.stateAtBeamLine().transverseImpactParameter().value();
        candidate.leptonAbsD0SignificanceH = fabs(leptons_[iHighPt].transientTrack.stateAtBeamLine().transverseImpactParameter().significance());
      } catch (cms::Exception& e) {};

      candidate.leptonIsoL = trackerIsolation(iEvent,leptons_[iLowPt].transientTrack.track(),leptons_[iLowPt].transientTrack.track());
      candidate.leptonIsoH = trackerIsolation(iEvent,leptons_[iHighPt].transientTrack.track(),leptons_[iHighPt].transientTrack.track());

      // FIXME: this can be improved
      candidate.triggerMatchL = 0;
      candidate.triggerMatchH = 0;
      if( leptons_[iLowPt].triggerMatch != leptons_[iHighPt].triggerMatch ) {
        // if( leptonTriggerMatch_[iLowPt] != 0 ) candidate.triggerMatchL = 1;
        // if( leptonTriggerMatch_[iHighPt] != 0 ) candidate.triggerMatchH = 2;
        if( leptons_[iLowPt].triggerMatch != 0 ) candidate.triggerMatchL = 1;
        if( leptons_[iHighPt].triggerMatch != 0 ) candidate.triggerMatchH = 2;
      }

      reco::TrackBase::Vector mom1 = lowPtLepton.momentum();
      reco::TrackBase::Vector mom2 = highPtLepton.momentum();
      candidate.cosine = mom1.Dot(mom2)/mom1.R()/mom2.R();
      candidate.deltaR = deltaR<reco::TrackBase::Vector,reco::TrackBase::Vector>(mom1,mom2);

      // ---------- //
      // Vertex fit //
      // ---------- //

      // attempt vertex fit with our lepton pair
      std::vector<reco::TransientTrack> fitTracks;
      fitTracks.push_back(leptons_[iLowPt].transientTrack);
      fitTracks.push_back(leptons_[iHighPt].transientTrack);
      bool validTracks = (fitTracks[0].isValid() && fitTracks[1].isValid()); 
      if( !validTracks ) continue;

      candidate.mass = (*(part1.p4)+*(part2.p4)).mass();

      if (validTracks) {
        TransientVertex leptonVertex;
        try {
          leptonVertex = vertexFitter_->vertex(fitTracks);
          candidate.validVertex=leptonVertex.isValid();
        } catch(cms::Exception& e) {};

        if (leptonVertex.isValid()) {
          candidate.vertexChi2 = leptonVertex.normalisedChiSquared();
          // calculate p_t corrected invariant mass
          ROOT::Math::PxPyPzMVector totalMom(0,0,0,0);
          ROOT::Math::PxPyPzMVector triggerCorrectedMom(0,0,0,0);
          ROOT::Math::PxPyPzMVector caloCorrectedMom(0,0,0,0);
          // const std::vector<reco::TransientTrack> tracks = leptonVertex.refittedTracks();
          reco::TransientTrack refittedTrackLowPt(leptonVertex.refittedTrack(leptons_[iLowPt].transientTrack));
          reco::TransientTrack refittedTrackHighPt(leptonVertex.refittedTrack(leptons_[iHighPt].transientTrack));
          double leptonMass=0.000511;
          if (lepType==__lepTypeMuon) leptonMass=0.106;
          if (lepType==__lepTypeTau) leptonMass=1.777;
          // Note that this is done for simplicity. In case of more than two tracks the main loop must be changed and
          // this can also become a loop on the refitted tracks.

          ROOT::Math::PxPyPzMVector refittedTrackMomLowPt(refittedTrackLowPt.track().px(),refittedTrackLowPt.track().py(),
                                                          refittedTrackLowPt.track().pz(),leptonMass);
          ROOT::Math::PxPyPzMVector refittedTrackMomHighPt(refittedTrackHighPt.track().px(),refittedTrackHighPt.track().py(),
                                                           refittedTrackHighPt.track().pz(),leptonMass);
          // Refitted momentum
          totalMom = refittedTrackMomLowPt + refittedTrackMomHighPt;
          // Trigger corrected momentum
          if (leptons_[iLowPt].triggerMatch != 0) triggerCorrectedMom += *(leptons_[iLowPt].triggerMatch);
          else triggerCorrectedMom += refittedTrackMomLowPt;
          if (leptons_[iHighPt].triggerMatch != 0) triggerCorrectedMom += *(leptons_[iHighPt].triggerMatch);
          else triggerCorrectedMom += refittedTrackMomHighPt;
          // Calo corrected momentum
          if( leptons_[iLowPt].matchedSC != 0 ) caloCorrectedMom +=
              (leptons_[iLowPt].matchedSC->energy()/leptons_[iLowPt].transientTrack.track().p())*(*(leptons_[iLowPt].p4));
          else caloCorrectedMom += refittedTrackMomLowPt;
          if( leptons_[iHighPt].matchedSC != 0 ) caloCorrectedMom +=
              (leptons_[iHighPt].matchedSC->energy()/leptons_[iHighPt].transientTrack.track().p())*(*(leptons_[iHighPt].p4));
          else caloCorrectedMom += refittedTrackMomHighPt;

          candidate.corrDileptonMass = totalMom.M();
          candidate.triggerCorrDileptonMass = triggerCorrectedMom.M();
          candidate.caloCorrMass = caloCorrectedMom.M();

          // decay length and significance
          // FIXME: this could be not very appropriate for standAloneMuons with big displacements (>~ 1 m).
          VertexDistanceXY vertTool;
          Measurement1D vertexDist = vertTool.distance(leptonVertex.vertexState(),*primaryVertex_);

          candidate.decayLength = vertexDist.value();
          candidate.decayLengthSignificance = vertexDist.significance();

          // vertex flight direction vs dilepton momentum direction
          GlobalVector vertexDir = leptonVertex.position() - primaryVertex_->position();
          math::XYZTLorentzVector diLeptonMom = *(part1.p4)+*(part2.p4);
          candidate.dPhi= fabs(deltaPhi<GlobalVector, math::XYZTLorentzVector>(vertexDir, diLeptonMom));
          candidate.dPhiCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, totalMom));
          candidate.dPhiTriggerCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, triggerCorrectedMom));
          candidate.dPhiCaloCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, caloCorrectedMom));

	  candidate.vx = leptonVertex.position().x();
	  candidate.vy = leptonVertex.position().y();
	  candidate.vz = leptonVertex.position().z();

          // general kinematic information about candidate
          candidate.pt=diLeptonMom.pt();
          candidate.eta=diLeptonMom.eta();
          candidate.phi=diLeptonMom.phi();
          candidate.ptCorr = totalMom.pt();
          candidate.etaCorr = totalMom.eta();
          candidate.phiCorr = totalMom.phi();
          candidate.ptTriggerCorr = triggerCorrectedMom.pt();
          candidate.etaTriggerCorr = triggerCorrectedMom.eta();
          candidate.phiTriggerCorr = triggerCorrectedMom.phi();
          candidate.ptCaloCorr = caloCorrectedMom.pt();
          candidate.etaCaloCorr = caloCorrectedMom.eta();
          candidate.phiCaloCorr = caloCorrectedMom.phi();

          // scale momenta to force dR=0 (phi component only).
          // we scale up because we assume
          // the momentum is underestimated due to energy loss
          double momScale=-(mom1.y()*vertexDir.x()-mom1.x()*vertexDir.y())/(mom2.y()*vertexDir.x()-mom2.x()*vertexDir.y());
          if (momScale>1) {
            candidate.scaleCorrMass=((*(part1.p4))+momScale*(*(part2.p4))).mass();
          } else {
            candidate.scaleCorrMass=((*(part1.p4))/momScale+(*(part2.p4))).mass();
          }

          // angle between positive lepton momentum in dilepton rest frame and dilepton momentum
          math::XYZTLorentzVector positiveMom;
          if (part1.transientTrack.track().charge()>0) positiveMom = *(part1.p4); else positiveMom = *(part2.p4);
          ROOT::Math::Boost boost;
          boost.SetComponents(-diLeptonMom/diLeptonMom.E());
          math::XYZTLorentzVector cmsMom = boost(positiveMom);
          candidate.cosThetaStar = TMath::Cos(angle<math::XYZTLorentzVector,math::XYZTLorentzVector>
                                          (diLeptonMom,cmsMom));

          double leptonPt[2];
          int hitsBeforeVertex[2];
          int missedLayersAfterVertex[2];
          if( fitTracks.size() != 2 ) {
            std::cout << "fitTracks size = " << fitTracks.size() << std::endl;
          }

          // Check hit pattern
          for (unsigned i=0; i<fitTracks.size(); i++) {
            std::pair<unsigned int, unsigned int> hitPat
                = checkHitPattern_->analyze(fitTracks[i].track(),leptonVertex);
            // first number: number of hits before vertex
            // second number: number of missed layers after vertex
            leptonPt[i] = fitTracks[i].track().pt();
            hitsBeforeVertex[i] = int(hitPat.first);
            missedLayersAfterVertex[i] = int(hitPat.second);
          }
          int iPtL = 0;
          int iPtH = 1;
          if( leptonPt[1] < leptonPt[0] ) {
            iPtL = 1;
            iPtH = 0;
          }
          candidate.hitsBeforeVertexL = hitsBeforeVertex[iPtL];
          candidate.hitsBeforeVertexH = hitsBeforeVertex[iPtH];
          candidate.missedLayersAfterVertexL = missedLayersAfterVertex[iPtL];
          candidate.missedLayersAfterVertexH = missedLayersAfterVertex[iPtH];
        }
      }

      candidate.isStandAloneL = leptons_[iLowPt].isStandAloneMuon;
      candidate.isStandAloneH = leptons_[iHighPt].isStandAloneMuon;
      candidate.isGlobalMuonL = leptons_[iLowPt].isGlobalMuon;
      candidate.isGlobalMuonH = leptons_[iHighPt].isGlobalMuon;
      candidate.isTrackerMuonL = leptons_[iLowPt].isTrackerMuon;
      candidate.isTrackerMuonH = leptons_[iHighPt].isTrackerMuon;
      if( leptons_[iLowPt].genPart != 0 ) candidate.pdgIdL = leptons_[iLowPt].genPart->pdgId();
      if( leptons_[iHighPt].genPart != 0 ) candidate.pdgIdH = leptons_[iHighPt].genPart->pdgId();
      if( leptons_[iLowPt].motherPart != 0 ) candidate.originPdgIdL = leptons_[iLowPt].motherPart->pdgId();
      if( leptons_[iHighPt].motherPart != 0 ) candidate.originPdgIdH = leptons_[iHighPt].motherPart->pdgId();
      if( leptons_[iLowPt].motherPart != 0 && (leptons_[iLowPt].motherPart == leptons_[iHighPt].motherPart) ) {
        GenEventProperties::DecayLengthAndType dlt(GenEventProperties::getDecayLengthAndType(*(leptons_[iLowPt].genPart)));
        candidate.genDecayLength2D = dlt.decayLength2D;
        candidate.genDecayLength3D = dlt.decayLength3D;
        candidate.genctau = dlt.ctau;
      }
      // FIXME: Improve calo matching part
      if( leptons_[iLowPt].matchedSC != 0 ) candidate.hasCaloMatchL = 1;
      if( leptons_[iHighPt].matchedSC != 0 ) candidate.hasCaloMatchH = 1;

    }
  }

  // Fill the tree for each candidate
  // std::cout << "Filling tree" << std::endl;
  // std::cout << "Filled outputTree" << std::endl;
  outputTree_->Fill();
}


const reco::Candidate* LeptonAnalysis::signalOrigin(const reco::Candidate* part)
{
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

int LeptonAnalysis::decayChannel(const reco::Candidate& part)
{
  // first check direct descendants
  for (unsigned i=0; i<part.numberOfDaughters(); i++) {
    const reco::Candidate* daughter = part.daughter(i);
    int pid = abs(daughter->pdgId());
    if (pid==11 || pid==13 || pid==15) return pid;
  }
  return 0;
}

void LeptonAnalysis::beginJob()
{
  vertexFitter_ = new KalmanVertexFitter(true);
}


void LeptonAnalysis::endJob() 
{
  // cut flow histograms
  std::string name=leptonName_;

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

  hTotalProcessedEvents_->SetBinContent(1, numProcessedEvents_);
  hEventsPassingTrigger_->SetBinContent(1, numEventsPassingTrigger_);
  fs_->cd();
  hTotalProcessedEvents_->Write();
  hEventsPassingTrigger_->Write();
}

void LeptonAnalysis::initializeVars()
{
  candidates_.clear();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LeptonAnalysis);
