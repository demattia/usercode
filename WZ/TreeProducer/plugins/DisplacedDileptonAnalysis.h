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


class DisplacedDileptonAnalysis : public edm::EDAnalyzer {

public:
  explicit DisplacedDileptonAnalysis(const edm::ParameterSet&);
  ~DisplacedDileptonAnalysis();
  
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
      isStandAloneMuon(false)
    {}

    reco::TransientTrack transientTrack;
    const reco::Particle::LorentzVector * p4;
    const reco::SuperCluster * matchedSC;
    const reco::GenParticle * genPart;
    const reco::Candidate * motherPart;
    const reco::Particle::LorentzVector * triggerMatch;
    bool isStandAloneMuon;
  };
  std::vector<LeptonContainer> leptons_;

  // Structure holding all the variables
  // ------------------------------------------------------------------------------------------------- //
  // IMPORTANT: when adding new variables do not forget to initalize them in the initializeVars method
  // and to add the branch in the tree.
  // ------------------------------------------------------------------------------------------------- //
  struct Variables
  {
    // Single lepton variables
    Float_t leptonPtL;
    Float_t leptonPtH;
    Float_t leptonEtaL;
    Float_t leptonEtaH;
    Int_t leptonChargeL;
    Int_t leptonChargeH;
    Float_t leptonD0L;
    Float_t leptonD0H;
    Float_t leptonAbsD0SignificanceL;
    Float_t leptonAbsD0SignificanceH;
    Float_t leptonIsoL;
    Float_t leptonIsoH;
    Int_t triggerMatchL;
    Int_t triggerMatchH;
    // Di-lepton variables
    Float_t cosine; // to reject cosmics. Leave unbiased by taking the values before refit to vertex
    Float_t deltaR; // for trigger inefficiencies. Leave unbiased as above.
    Float_t cosThetaStar; // angle between positive lepton momentum in dilepton rest frame and dilepton momentum
    // In the following:
    // - the default values are from the leptons
    // - the "corr" values are from the leptons after refit to common vertex
    // - the "triggerCorr" values are from the trigger matches
    // - the "caloCorr" values are from the calo matches
    Float_t mass;
    Float_t corrDileptonMass;
    Float_t triggerCorrDileptonMass;
    Float_t caloCorrMass;
    Float_t scaleCorrMass; // corrected so that deltaR between vertex flight direction and di-lepton momentum is 0.
    Float_t eta;
    Float_t phi;
    Float_t etaCorr;
    Float_t phiCorr;
    Float_t etaTriggerCorr;
    Float_t phiTriggerCorr;
    Float_t etaCaloCorr;
    Float_t phiCaloCorr;

    Float_t decayLength;
    Float_t decayLengthSignificance;
    Float_t dPhi; // Angle in transverse plane between vertex (secondary-primay) flight direction and di-lepton momentum
    Float_t dPhiCorr; // Same as above but using leptons refitted to vertex to compute di-lepton momentum
    Float_t dPhiTriggerCorr; // Same as above using trigger matches to compute di-lepton momentum
    Float_t dPhiCaloCorr;

    Float_t hitsBeforeVertexL;
    Float_t hitsBeforeVertexH;
    Float_t missedLayersAfterVertexL;
    Float_t missedLayersAfterVertexH;
    Int_t isStandAloneL;
    Int_t isStandAloneH;
    Int_t hasCaloMatchL;
    Int_t hasCaloMatchH;
    Int_t originPdgIdL;
    Int_t originPdgIdH;
    Int_t pdgIdL;
    Int_t pdgIdH;
    Float_t genctau;
    Float_t genDecayLength2D;
    Float_t genDecayLength3D;

    Int_t validVertex;
    Float_t vertexChi2;
    // Event related variables
    Int_t numPV;
    Int_t nvtx_m1;
    Int_t nvtx_0;
    Int_t nvtx_p1;
    Int_t run;
    Int_t event;
    // HLT paths triggered on in a particular event
    std::vector<std::string> triggers;
  } vars_;

};
