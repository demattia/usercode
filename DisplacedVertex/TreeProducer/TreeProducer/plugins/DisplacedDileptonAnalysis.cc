#include "TreeProducer/TreeProducer/plugins/DisplacedDileptonAnalysis.h"

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

void DisplacedDileptonAnalysis::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & setup) {
  // Total number of events is the sum of the events in each of these luminosity blocks
  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;
  lumi.getByLabel("nEventsTotal", nEventsTotalCounter);
  numProcessedEvents_ += nEventsTotalCounter->value;

  edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
  lumi.getByLabel("nEventsPostHLTFilter", nEventsFilteredCounter);
  numEventsPassingTrigger_ += nEventsFilteredCounter->value;
}

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
  outputTree_->Branch("leptonPtL",&vars_.leptonPtL,"leptonPtL/F");
  outputTree_->Branch("leptonPtH",&vars_.leptonPtH,"leptonPtH/F");
  outputTree_->Branch("leptonEtaL",&vars_.leptonEtaL,"leptonEtaL/F");
  outputTree_->Branch("leptonEtaH",&vars_.leptonEtaH,"leptonEtaH/F");
  outputTree_->Branch("leptonChargeL",&vars_.leptonChargeL,"leptonChargeL/I");
  outputTree_->Branch("leptonChargeH",&vars_.leptonChargeH,"leptonChargeH/I");
  outputTree_->Branch("leptonD0L",&vars_.leptonD0L,"leptonD0L/F");
  outputTree_->Branch("leptonD0H",&vars_.leptonD0H,"leptonD0H/F");
  outputTree_->Branch("leptonAbsD0SignificanceL",&vars_.leptonAbsD0SignificanceL,"leptonAbsD0SignificanceL/F");
  outputTree_->Branch("leptonAbsD0SignificanceH",&vars_.leptonAbsD0SignificanceH,"leptonAbsD0SignificanceH/F");
  outputTree_->Branch("leptonIsoL",&vars_.leptonIsoL,"leptonIsoL/F");
  outputTree_->Branch("leptonIsoH",&vars_.leptonIsoH,"leptonIsoH/F");
  outputTree_->Branch("triggerMatchL",&vars_.triggerMatchL,"triggerMatchL/I");
  outputTree_->Branch("triggerMatchH",&vars_.triggerMatchH,"triggerMatchH/I");
  // Di-lepton variables
  outputTree_->Branch("cosine",&vars_.cosine,"cosine/F");
  outputTree_->Branch("deltaR",&vars_.deltaR,"deltaR/F");
  outputTree_->Branch("cosThetaStar",&vars_.cosThetaStar,"cosThetaStar/F");
  outputTree_->Branch("mass",&vars_.mass,"mass/F");
  outputTree_->Branch("corrDileptonMass",&vars_.corrDileptonMass,"corrDileptonMass/F");
  outputTree_->Branch("triggerCorrDileptonMass",&vars_.triggerCorrDileptonMass,"triggerCorrDileptonMass/F");
  outputTree_->Branch("caloCorrMass",&vars_.caloCorrMass,"caloCorrMass/F");
  outputTree_->Branch("scaleCorrMass",&vars_.scaleCorrMass,"scaleCorrMass/F");
  outputTree_->Branch("eta",&vars_.eta,"eta/F");
  outputTree_->Branch("phi",&vars_.phi,"phi/F");
  outputTree_->Branch("etaCorr",&vars_.etaCorr,"etaCorr/F");
  outputTree_->Branch("phiCorr",&vars_.phiCorr,"phiCorr/F");
  outputTree_->Branch("etaTriggerCorr",&vars_.etaTriggerCorr,"etaTriggerCorr/F");
  outputTree_->Branch("phiTriggerCorr",&vars_.phiTriggerCorr,"phiTriggerCorr/F");
  outputTree_->Branch("etaCaloCorr",&vars_.etaCaloCorr,"etaCaloCorr/F");
  outputTree_->Branch("phiCaloCorr",&vars_.phiCaloCorr,"phiCaloCorr/F");
  outputTree_->Branch("decayLength",&vars_.decayLength,"decayLength/F");
  outputTree_->Branch("decayLengthSignificance",&vars_.decayLengthSignificance,"decayLengthSignificance/F");
  outputTree_->Branch("dPhi",&vars_.dPhi,"dPhi/F");
  outputTree_->Branch("dPhiCorr",&vars_.dPhiCorr,"dPhiCorr/F");
  outputTree_->Branch("dPhiTriggerCorr",&vars_.dPhiTriggerCorr,"dPhiTriggerCorr/F");
  outputTree_->Branch("dPhiCaloCorr",&vars_.dPhiCaloCorr,"dPhiCaloCorr/F");
  // Single lepton variables from refitted vertex
  outputTree_->Branch("hitsBeforeVertexL",&vars_.hitsBeforeVertexL,"hitsBeforeVertexL/F");
  outputTree_->Branch("hitsBeforeVertexH",&vars_.hitsBeforeVertexH,"hitsBeforeVertexH/F");
  outputTree_->Branch("missedLayersAfterVertexL",&vars_.missedLayersAfterVertexL,"missedLayersAfterVertexL/F");
  outputTree_->Branch("missedLayersAfterVertexH",&vars_.missedLayersAfterVertexH,"missedLayersAfterVertexH/F");
  outputTree_->Branch("isStandAloneL",&vars_.isStandAloneL,"isStandAloneL/I");
  outputTree_->Branch("isStandAloneH",&vars_.isStandAloneH,"isStandAloneH/I");
  outputTree_->Branch("hasCaloMatchL",&vars_.hasCaloMatchL,"hasCaloMatchL/I");
  outputTree_->Branch("hasCaloMatchH",&vars_.hasCaloMatchH,"hasCaloMatchH/I");
  outputTree_->Branch("validVertex",&vars_.validVertex,"validVertex/I");
  outputTree_->Branch("vertexChi2",&vars_.vertexChi2,"vertexChi2/F");
  // Generator level information
  outputTree_->Branch("originPdgIdL",&vars_.originPdgIdL,"originPdgIdL/I");
  outputTree_->Branch("originPdgIdH",&vars_.originPdgIdH,"originPdgIdH/I");
  outputTree_->Branch("pdgIdL",&vars_.pdgIdL,"pdgIdL/I");
  outputTree_->Branch("pdgIdH",&vars_.pdgIdH,"pdgIdH/I");
  outputTree_->Branch("genctau",&vars_.genctau, "genctau/F");
  outputTree_->Branch("genDecayLength2D",&vars_.genDecayLength2D, "genDecayLength2D/F");
  outputTree_->Branch("genDecayLength3D",&vars_.genDecayLength3D, "genDecayLength3D/F");

  // Event information
  outputTree_->Branch("numPV",&vars_.numPV,"numPV/I");
  outputTree_->Branch("nvtx_m1",&vars_.nvtx_m1,"nvtx_m1/I");
  outputTree_->Branch("nvtx_0",&vars_.nvtx_0,"nvtx_0/I");
  outputTree_->Branch("nvtx_p1",&vars_.nvtx_p1,"nvtx_p1/I");
  outputTree_->Branch("run",&vars_.run,"run/I");
  outputTree_->Branch("event",&vars_.event,"event/I");
  outputTree_->Branch("triggers","std::vector<std::string>",&vars_.triggers);

  // To save the total number of processed events
  hTotalProcessedEvents_ = new TH1F("totalProcessedEvents", "Total processed events", 1, 0, 1);
  hEventsPassingTrigger_ = new TH1F("eventsPassingTrigger", "Events passing trigger selection", 1, 0, 1);
}


DisplacedDileptonAnalysis::~DisplacedDileptonAnalysis()
{
}


void DisplacedDileptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    if (thisLepton_==__lepTypeElectron) {
      doDecayChannel<pat::Electron>(iEvent,__lepTypeElectron);
    } else if (thisLepton_==__lepTypeMuon) {
      doDecayChannel<pat::Muon>(iEvent,__lepTypeMuon);
    } else if (thisLepton_==__lepTypeTau) {
      doDecayChannel<pat::Tau>(iEvent,__lepTypeTau);
    } else {
      // doDecayChannel<PseudoLepton>(iEvent,__lepTypeTrack,genProp);
      doDecayChannel<PseudoLepton>(iEvent,__lepTypeTrack);
    }
  }
}

template<class Lepton>
const reco::Particle::LorentzVector * DisplacedDileptonAnalysis::findTrigMatch(const Lepton & lepton,
                                                                               const edm::Handle< pat::TriggerEvent > & triggerEvent)
{
  // find deltaR match among trigger objects
  double lep_eta = lepton.eta();
  double lep_phi = lepton.phi();
  double min_deltaR=9999.;
  const reco::Particle::LorentzVector* trig_mom=0;

  // loop over all triggers we are interested in
  for (unsigned itrig=0; itrig<hltPaths_.size(); itrig++) {
    // std::cout << "Trigger path number " << itrig << std::endl;
    if (!triggerEvent->path(hltPaths_[itrig])) continue;
    if (!triggerEvent->path(hltPaths_[itrig])->wasAccept()) continue;
    // get all trigger objects associated with this HLT path
    const pat::TriggerObjectRefVector triggerObjects(triggerEvent->pathObjects(hltPaths_[itrig],false));
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
const reco::SuperCluster * DisplacedDileptonAnalysis::findCaloMatch(const Lepton & lepton,
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
bool DisplacedDileptonAnalysis::doTrigger(const edm::Event& iEvent)
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
  vars_.triggers.clear();
  for (unsigned i=0; i<triggerEvent->acceptedPaths().size(); i++) {
    std::string trigname(triggerEvent->acceptedPaths()[i]->name());
    ++trigPathSummary_[trigname];
    for (unsigned itrg=0; itrg<hltPaths_.size(); itrg++) {
      if( takeAllPaths || trigname==hltPaths_[itrg] ) {
        // FIXME: might want to select a reduced number of triggers to save space in the tree.
        vars_.triggers.push_back(trigname);
      }
    }
  }
  return (vars_.triggers.size()>0);
}

/// Return the closest stable genParticle in deltaR
const reco::GenParticle* DisplacedDileptonAnalysis::getGenParticleMatch(const edm::Handle<edm::View<reco::GenParticle> > & mcParticles,
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
double DisplacedDileptonAnalysis::trackerIsolation(const edm::Event& iEvent,
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
reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const PseudoLepton &track)
{
  reco::TransientTrack ttrack=trackBuilder_->build(track);
  ttrack.setBeamSpot(beamSpot_);
  return ttrack;
}

/// Convert to a transient track and assign the beamspot.
reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Electron &electron)
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
reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Muon &muon)
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
reco::TransientTrack DisplacedDileptonAnalysis::leptonTrack(const pat::Tau &tau)
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
bool DisplacedDileptonAnalysis::leptonID(const Lepton& lepton)
{
  return true;
}

/// override template for the electron channel, where we do use a specific ID requirement
bool DisplacedDileptonAnalysis::leptonID(const pat::Electron& electron)
{
  return int(electron.electronID("eidLoose"))&1;
}

/// Default method, returns -999
double DisplacedDileptonAnalysis::leptonTiming(const PseudoLepton &track)
{
  return -999;
}

/// Not implemented, returns -999
double DisplacedDileptonAnalysis::leptonTiming(const pat::Electron &electron)
{
  //DetId seedDetID = electron->superCluster()->seed()->seed();
  return -999;
}

/// Muon timing. Assumes muon coming from inside the detector and has beta == 1
double DisplacedDileptonAnalysis::leptonTiming(const pat::Muon &muon)
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
double DisplacedDileptonAnalysis::leptonTiming(const pat::Tau &tau) {
  // not implemented
  return -999;
}


template<class Lepton>
void DisplacedDileptonAnalysis::doDecayChannel(const edm::Event& iEvent,
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

    leptons_.push_back(lepton);
  }

  // Find candidates
  // ---------------

  // loop over the saved leptons
  for (unsigned ipart1=0; ipart1<leptons_.size(); ++ipart1) {
    const LeptonContainer part1 = leptons_[ipart1];
    for (unsigned ipart2=ipart1+1; ipart2<leptons_.size(); ++ipart2) {
      const LeptonContainer part2 = leptons_[ipart2];

      // Initialize all variables
      initializeVars();
      // std::cout << "variables initialized" << std::endl;
      vars_.run = iEvent.run();
      vars_.event = iEvent.id().event();

      unsigned iHighPt=ipart1;
      unsigned iLowPt=ipart2;
      if (part1.transientTrack.track().pt()<part2.transientTrack.track().pt()) {
        iHighPt=ipart2;
        iLowPt=ipart1;
      }
      reco::Track lowPtLepton = leptons_[iLowPt].transientTrack.track();
      reco::Track highPtLepton = leptons_[iHighPt].transientTrack.track();
      vars_.leptonPtL = lowPtLepton.pt();
      vars_.leptonPtH = highPtLepton.pt();
      vars_.leptonEtaL = lowPtLepton.eta();
      vars_.leptonEtaH = highPtLepton.eta();
      vars_.leptonChargeL = lowPtLepton.charge();
      vars_.leptonChargeH = highPtLepton.charge();

      try {
        vars_.leptonD0L = leptons_[iLowPt].transientTrack.stateAtBeamLine().transverseImpactParameter().value();
        vars_.leptonAbsD0SignificanceL = fabs(leptons_[iLowPt].transientTrack.stateAtBeamLine().transverseImpactParameter().significance());
        vars_.leptonD0H = leptons_[iHighPt].transientTrack.stateAtBeamLine().transverseImpactParameter().value();
        vars_.leptonAbsD0SignificanceH = fabs(leptons_[iHighPt].transientTrack.stateAtBeamLine().transverseImpactParameter().significance());
      } catch (cms::Exception& e) {};

      vars_.leptonIsoL = trackerIsolation(iEvent,leptons_[iLowPt].transientTrack.track(),leptons_[iLowPt].transientTrack.track());
      vars_.leptonIsoH = trackerIsolation(iEvent,leptons_[iHighPt].transientTrack.track(),leptons_[iHighPt].transientTrack.track());

      // FIXME: this can be improved
      vars_.triggerMatchL = 0;
      vars_.triggerMatchH = 0;
      if( leptons_[iLowPt].triggerMatch != leptons_[iHighPt].triggerMatch ) {
        // if( leptonTriggerMatch_[iLowPt] != 0 ) vars_.triggerMatchL = 1;
        // if( leptonTriggerMatch_[iHighPt] != 0 ) vars_.triggerMatchH = 2;
        if( leptons_[iLowPt].triggerMatch != 0 ) vars_.triggerMatchL = 1;
        if( leptons_[iHighPt].triggerMatch != 0 ) vars_.triggerMatchH = 2;
      }

      reco::TrackBase::Vector mom1 = lowPtLepton.momentum();
      reco::TrackBase::Vector mom2 = highPtLepton.momentum();
      vars_.cosine = mom1.Dot(mom2)/mom1.R()/mom2.R();
      vars_.deltaR = deltaR<reco::TrackBase::Vector,reco::TrackBase::Vector>(mom1,mom2);

      // ---------- //
      // Vertex fit //
      // ---------- //

      // attempt vertex fit with our lepton pair
      std::vector<reco::TransientTrack> fitTracks;
      fitTracks.push_back(leptons_[iLowPt].transientTrack);
      fitTracks.push_back(leptons_[iHighPt].transientTrack);
      bool validTracks = (fitTracks[0].isValid() && fitTracks[1].isValid()); 
      if( !validTracks ) continue;

      vars_.mass = (*(part1.p4)+*(part2.p4)).mass();

      if (validTracks) {
        TransientVertex leptonVertex;
        try {
          leptonVertex = vertexFitter_->vertex(fitTracks);
          vars_.validVertex=leptonVertex.isValid();
        } catch(cms::Exception& e) {};

        if (leptonVertex.isValid()) {
          vars_.vertexChi2 = leptonVertex.normalisedChiSquared();
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

          vars_.corrDileptonMass = totalMom.M();
          vars_.triggerCorrDileptonMass = triggerCorrectedMom.M();
          vars_.caloCorrMass = caloCorrectedMom.M();

          // decay length and significance
          // FIXME: this could be not very appropriate for standAloneMuons with big displacements (>~ 1 m).
          VertexDistanceXY vertTool;
          Measurement1D vertexDist = vertTool.distance(leptonVertex.vertexState(),*primaryVertex_);

          vars_.decayLength = vertexDist.value();
          vars_.decayLengthSignificance = vertexDist.significance();

          // vertex flight direction vs dilepton momentum direction
          GlobalVector vertexDir = leptonVertex.position() - primaryVertex_->position();
          math::XYZTLorentzVector diLeptonMom = *(part1.p4)+*(part2.p4);
          vars_.dPhi= fabs(deltaPhi<GlobalVector, math::XYZTLorentzVector>(vertexDir, diLeptonMom));
          vars_.dPhiCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, totalMom));
          vars_.dPhiTriggerCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, triggerCorrectedMom));
          vars_.dPhiCaloCorr = fabs(deltaPhi<GlobalVector, ROOT::Math::PxPyPzMVector>(vertexDir, caloCorrectedMom));

          // general kinematic information about candidate
          vars_.eta=diLeptonMom.eta();
          vars_.phi=diLeptonMom.phi();
          vars_.etaCorr = totalMom.eta();
          vars_.phiCorr = totalMom.phi();
          vars_.etaTriggerCorr = triggerCorrectedMom.eta();
          vars_.phiTriggerCorr = triggerCorrectedMom.phi();
          vars_.etaCaloCorr = caloCorrectedMom.eta();
          vars_.phiCaloCorr = caloCorrectedMom.phi();

          // scale momenta to force dR=0 (phi component only).
          // we scale up because we assume
          // the momentum is underestimated due to energy loss
          double momScale=-(mom1.y()*vertexDir.x()-mom1.x()*vertexDir.y())/(mom2.y()*vertexDir.x()-mom2.x()*vertexDir.y());
          if (momScale>1) {
            vars_.scaleCorrMass=((*(part1.p4))+momScale*(*(part2.p4))).mass();
          } else {
            vars_.scaleCorrMass=((*(part1.p4))/momScale+(*(part2.p4))).mass();
          }

          // angle between positive lepton momentum in dilepton rest frame and dilepton momentum
          math::XYZTLorentzVector positiveMom;
          if (part1.transientTrack.track().charge()>0) positiveMom = *(part1.p4); else positiveMom = *(part2.p4);
          ROOT::Math::Boost boost;
          boost.SetComponents(-diLeptonMom/diLeptonMom.E());
          math::XYZTLorentzVector cmsMom = boost(positiveMom);
          vars_.cosThetaStar = TMath::Cos(angle<math::XYZTLorentzVector,math::XYZTLorentzVector>
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
          vars_.hitsBeforeVertexL = hitsBeforeVertex[iPtL];
          vars_.hitsBeforeVertexH = hitsBeforeVertex[iPtH];
          vars_.missedLayersAfterVertexL = missedLayersAfterVertex[iPtL];
          vars_.missedLayersAfterVertexH = missedLayersAfterVertex[iPtH];
        }
      }

      vars_.numPV = numPV_;
      if( useMCTruth_ ) {
        vars_.nvtx_m1 = nvtx_m1_;
        vars_.nvtx_0 = nvtx_0_;
        vars_.nvtx_p1 = nvtx_p1_;
      }

      vars_.isStandAloneL = leptons_[iLowPt].isStandAloneMuon;
      vars_.isStandAloneH = leptons_[iHighPt].isStandAloneMuon;
      if( leptons_[iLowPt].genPart != 0 ) vars_.pdgIdL = leptons_[iLowPt].genPart->pdgId();
      if( leptons_[iHighPt].genPart != 0 ) vars_.pdgIdH = leptons_[iHighPt].genPart->pdgId();
      if( leptons_[iLowPt].motherPart != 0 ) vars_.originPdgIdL = leptons_[iLowPt].motherPart->pdgId();
      if( leptons_[iHighPt].motherPart != 0 ) vars_.originPdgIdH = leptons_[iHighPt].motherPart->pdgId();
      if( leptons_[iLowPt].motherPart != 0 && (leptons_[iLowPt].motherPart == leptons_[iHighPt].motherPart) ) {
        GenEventProperties::DecayLengthAndType dlt(GenEventProperties::getDecayLengthAndType(*(leptons_[iLowPt].genPart)));
        vars_.genDecayLength2D = dlt.decayLength2D;
        vars_.genDecayLength3D = dlt.decayLength3D;
        vars_.genctau = dlt.ctau;
      }
      // FIXME: Improve calo matching part
      if( leptons_[iLowPt].matchedSC != 0 ) vars_.hasCaloMatchL = 1;
      if( leptons_[iHighPt].matchedSC != 0 ) vars_.hasCaloMatchH = 1;

      // Fill the tree for each candidate
      // std::cout << "Filling tree" << std::endl;
      outputTree_->Fill();
      // std::cout << "Filled outputTree" << std::endl;
    }
  }
}


const reco::Candidate* DisplacedDileptonAnalysis::signalOrigin(const reco::Candidate* part)
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


int DisplacedDileptonAnalysis::decayChannel(const reco::Candidate& part) {

  // first check direct descendants
  for (unsigned i=0; i<part.numberOfDaughters(); i++) {
    const reco::Candidate* daughter = part.daughter(i);
    int pid = abs(daughter->pdgId());
    if (pid==11 || pid==13 || pid==15) return pid;
  }
  return 0;
}

void DisplacedDileptonAnalysis::beginJob()
{
  vertexFitter_ = new KalmanVertexFitter(true);
}


void DisplacedDileptonAnalysis::endJob() 
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

void DisplacedDileptonAnalysis::initializeVars()
{
  // Single lepton variables
  vars_.leptonPtL = -999;
  vars_.leptonPtH = -999;
  vars_.leptonEtaL = -999;
  vars_.leptonEtaH = -999;
  vars_.leptonChargeL = -999;
  vars_.leptonChargeH = -999;
  vars_.leptonD0L = -999;
  vars_.leptonD0H= -999;
  vars_.leptonAbsD0SignificanceL = -999;
  vars_.leptonAbsD0SignificanceH = -999;
  vars_.leptonIsoL = -999;
  vars_.leptonIsoH = -999;
  vars_.triggerMatchL = -999;
  vars_.triggerMatchH = -999;
  // Di-lepton variables
  vars_.cosine = -999; // to reject cosmics. Leave unbiased by taking the values before refit to vertex
  vars_.deltaR = -999; // for trigger inefficiencies. Leave unbiased as above.
  vars_.cosThetaStar = -999; // angle between positive lepton momentum in dilepton rest frame and dilepton momentum

  vars_.mass = -999;
  vars_.corrDileptonMass = -999;
  vars_.triggerCorrDileptonMass = -999;
  vars_.caloCorrMass = -999;
  vars_.scaleCorrMass = -999; // corrected so that deltaR between vertex flight direction and di-lepton momentum is 0.
  vars_.eta = -999;
  vars_.phi = -999;
  vars_.etaCorr = -999;
  vars_.phiCorr = -999;
  vars_.etaTriggerCorr = -999;
  vars_.phiTriggerCorr = -999;
  vars_.etaCaloCorr = -999;
  vars_.phiCaloCorr = -999;

  vars_.decayLength = -999;
  vars_.decayLengthSignificance = -999;
  vars_.dPhi = -999; // Angle in transverse plane between vertex (secondary-primay) flight direction and di-lepton momentum
  vars_.dPhiCorr = -999; // Same as above but using leptons refitted to vertex to compute di-lepton momentum
  vars_.dPhiTriggerCorr = -999; // Same as above using trigger matches to compute di-lepton momentum
  vars_.dPhiCaloCorr = -999;

  vars_.hitsBeforeVertexL = -999;
  vars_.hitsBeforeVertexH = -999;
  vars_.missedLayersAfterVertexL = -999;
  vars_.missedLayersAfterVertexH = -999;
  vars_.isStandAloneL = -999;
  vars_.isStandAloneH = -999;
  vars_.hasCaloMatchL = -999;
  vars_.hasCaloMatchH = -999;
  vars_.originPdgIdL = -999;
  vars_.originPdgIdH = -999;
  vars_.genctau = -999;
  vars_.genDecayLength2D = -999;
  vars_.genDecayLength3D = -999;
  vars_.pdgIdL = -999;
  vars_.pdgIdH = -999;

  vars_.validVertex = -999;
  vars_.vertexChi2 = -999;

  vars_.numPV = -999;
  vars_.nvtx_m1 = -999;
  vars_.nvtx_0 = -999;
  vars_.nvtx_p1 = -999;
  vars_.run = -999;
  vars_.event = -999;

  // The trigger vector is handled by the anayzer
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DisplacedDileptonAnalysis);
