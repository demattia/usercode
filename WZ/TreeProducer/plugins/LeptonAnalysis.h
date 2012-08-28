#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TreeProducer/TreeProducer/interface/PseudoLepton.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TreeProducer/TreeProducer/interface/GenEventProperties.h"
#include "TreeProducer/CheckHitPattern/interface/CheckHitPattern.h"

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <map>

#include "TreeProducer/TreeProducer/interface/Candidates.h"
#include "TreeProducer/TreeProducer/interface/TreeCandidate.h"
#include "TreeProducer/TreeProducer/interface/TreeLepton.h"

class LeptonAnalysis : public edm::EDAnalyzer {

public:
  explicit LeptonAnalysis(const edm::ParameterSet&);
  ~LeptonAnalysis();
  
private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & setup);
  virtual void endJob() ;

  // decay channels
  typedef enum { __lepTypeElectron,
                 __lepTypeMuon,
                 __lepTypeTau,
                 __lepTypeTrack }
               leptonType;
  std::string leptonName_;

  // types of candidates that we distinguish
  typedef enum { __catSignal,             // fully reconstructed signal
                 __catSignalWrongChannel, // fully reconstructed signal, but wrong decay mode
                 __catSignalMixed,        // both leptons from signal, but wrong pairing
                 __catOneSigOneBkgWithGen,// one lepton from signal, one from unrelated GenParticle
                 __catOneSigOneBadReco,   // one lepton from signal, one without GenParticle at all
                 __catGoodButWrongType,   // fully reconstructed *other* resonance (e.g. Z)
                 __catBackground }        // both leptons are of non-signal origin
               categoryType;

  const reco::GenParticle* getGenParticleMatch(const edm::Handle<edm::View<reco::GenParticle> > & mcParticles,
                                               const reco::Particle::LorentzVector& mom);
  int decayChannel(const reco::Candidate&);
  const reco::Candidate* signalOrigin(const reco::Candidate*);

  bool doTrigger(const edm::Event&);
  template<class Lepton>
  const reco::Particle::LorentzVector * findTrigMatch(const Lepton & lepton,
                                                      const edm::Handle< pat::TriggerEvent > & triggerEvent);
  template<class Lepton>
  const reco::SuperCluster * findCaloMatch(const Lepton & lepton,
                                           const edm::Handle< std::vector<reco::SuperCluster> > & barrelSuperClusters,
                                           const edm::Handle< std::vector<reco::SuperCluster> > & endcapSuperClusters);
  template<class Lepton> void doDecayChannel(const edm::Event&, const leptonType);
  template<class Lepton> bool leptonID(const Lepton&);
  bool leptonID(const pat::Electron&);
  double trackerIsolation(const edm::Event& iEvent,
                          const reco::Track & plept,
                          const reco::Track & pveto);
  //  double trackerIsolation(const edm::Event&, const reco::Particle::LorentzVector&, const reco::Particle::LorentzVector&);
  reco::TransientTrack leptonTrack(const PseudoLepton &);
  reco::TransientTrack leptonTrack(const pat::Electron &);
  reco::TransientTrack leptonTrack(const pat::Muon &);
  reco::TransientTrack leptonTrack(const pat::Tau &);

  // detector timing
  double leptonTiming(const pat::Electron &);
  double leptonTiming(const pat::Muon &);
  double leptonTiming(const pat::Tau &);
  double leptonTiming(const PseudoLepton &);

  // Initialization
  void initializeVars();

  // beamspot
  reco::BeamSpot beamSpot_;

  // primary vertex
  VertexState* primaryVertex_;
  unsigned numPV_;
  unsigned nvtx_m1_;
  unsigned nvtx_0_;
  unsigned nvtx_p1_;

  // MET
  float MET_;
  float METPhi_;

  // for displaced vertex fit
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;
  KalmanVertexFitter* vertexFitter_;
  CheckHitPattern* checkHitPattern_;

  // container for additional signal particle information
  typedef struct {
    double rt;
    double r;
    double ndaughters;
  } GenInfo_t;
  std::map<const reco::Candidate*,GenInfo_t*> GenInfoMap_;

  // input tags  
  edm::InputTag leptonCollTag_[4];
  edm::InputTag generatorTag_;
  edm::InputTag pileupTag_;
  edm::InputTag genEventInfoTag_;
  edm::InputTag trigger_;
  edm::InputTag triggerEvent_;
  edm::InputTag barrelSuperClusters_;
  edm::InputTag endcapSuperClusters_;

  // analysis setup
  int signalPDGId_;
  int leptonPDGId_;
  leptonType thisLepton_;
  bool cutStudyMode_;
  std::vector<std::string> hltPaths_;
  int runMin_;
  int runMax_;

  // selection cuts
  double leptonPtCut_;
  bool leptonChargeCut_;

  // trigger matches storage
  std::vector<const reco::Particle::LorentzVector*> leptonTriggerMatch_;

  // calo supercluster storage
  std::vector<double> leptonCaloMatch_;

  bool useMCTruth_;

  // trigger summary storage
  unsigned numProcessedEvents_;
  unsigned numEventsPassingTrigger_;

  std::map<std::string,unsigned> trigFilterSummary_;
  std::map<std::string,unsigned> trigPathSummary_;

  // To count the total number of events on which the producer runs
  TH1F * hTotalProcessedEvents_;
  TH1F * hEventsPassingTrigger_;

  // The tree with all the variables
  edm::Service<TFileService> fs_;
  TTree* outputTree_;

  // To fill single lepton information
  struct LeptonContainer
  {
    LeptonContainer(const reco::TransientTrack & ttrack) :
      transientTrack(ttrack),
      p4(0),
      matchedSC(0),
      genPart(0),
      motherPart(0),
      triggerMatch(0),
      isStandAloneMuon(false),
      isGlobalMuon(false),
      isTrackerMuon(false),
      vx(-999.),
      vy(-999.),
      vz(-999.),
      trackIso(-999.),
      ecalIso(-999.),
      hcalIso(-999.),
      normChi2(-999.),
      numberOfValidTrackerHits(-999),
      numberOfValidPixelHits(-999),
      numberOfValidMuonHits(-999),
      numberOfMatchedStations(-999)
    {}

    reco::TransientTrack transientTrack;
    const reco::Particle::LorentzVector * p4;
    const reco::SuperCluster * matchedSC;
    const reco::GenParticle * genPart;
    const reco::Candidate * motherPart;
    const reco::Particle::LorentzVector * triggerMatch;
    bool isStandAloneMuon;
    bool isGlobalMuon;
    bool isTrackerMuon;
    float vx, vy, vz;
    float trackIso;
    float ecalIso;
    float hcalIso;
    double normChi2;
    int numberOfValidTrackerHits;
    int numberOfValidPixelHits;
    int numberOfValidMuonHits;
    int numberOfMatchedStations;
  };
  std::vector<LeptonContainer> leptons_;

  std::vector<std::string> triggers_;

  // Structure holding all the variables
  Candidates candidates_;

};
