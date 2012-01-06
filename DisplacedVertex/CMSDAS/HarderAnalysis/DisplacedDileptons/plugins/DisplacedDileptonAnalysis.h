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
#include "HarderAnalysis/DisplacedDileptons/interface/PseudoLepton.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DisplacedLeptons/Samples/interface/GenEventProperties.h"
#include "HarderAnalysisTools/Histograms/interface/HistMap.h"
#include "HarderAnalysisTools/CutManager/interface/CutSummary.h"
#include "HarderAnalysisTools/CheckHitPattern/interface/CheckHitPattern.h"

#include <map>


class DisplacedDileptonAnalysis : public edm::EDAnalyzer {

public:
  explicit DisplacedDileptonAnalysis(const edm::ParameterSet&);
  ~DisplacedDileptonAnalysis();
  
private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
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

  bool findMCParticles(const edm::Event&);
  const reco::GenParticle* getGenParticleMatch(const edm::Event&, const reco::Particle::LorentzVector&);
  int decayChannel(const reco::Candidate&);
  const reco::Candidate* signalOrigin(const reco::Candidate*);
  const reco::Candidate* commonOrigin(const reco::Candidate*, const reco::Candidate*);

  bool doTrigger(const edm::Event&);
  template<class LeptonCollection> void fillTrigMatches(const edm::Event&,leptonType);
  void fillCaloMatches(const edm::Event&);
  template<class Lepton> void doDecayChannel(const edm::Event&, const leptonType,
					     GenEventProperties&);
  template<class Lepton> bool leptonID(const Lepton&);
  bool leptonID(const pat::Electron&);
  template<class Lepton> void leptonPlots(const Lepton*, const reco::GenParticle*, const leptonType, const bool, const bool);
  template<class Lepton> bool leptonCuts(const Lepton&, CutFlow &);
  double trackerIsolation(const edm::Event&, const reco::Particle::LorentzVector&, const reco::Particle::LorentzVector&);
  reco::TransientTrack leptonTrack(const PseudoLepton &);
  reco::TransientTrack leptonTrack(const pat::Electron &);
  reco::TransientTrack leptonTrack(const pat::Muon &);
  reco::TransientTrack leptonTrack(const pat::Tau &);
  void investigateTracks(const edm::Event&);

  // detector timing
  double leptonTiming(const pat::Electron &);
  double leptonTiming(const pat::Muon &);
  double leptonTiming(const pat::Tau &);
  double leptonTiming(const PseudoLepton &);

  // beamspot
  reco::BeamSpot beamSpot_;

  // primary vertex
  VertexState* primaryVertex_;
  unsigned numPV_;

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

  // histogram manager
  HistMap* mcHistos_;
  HistMap* weightHistos_;
  HistMap* trigHistos_;
  HistMap* lepHistos_;
  HistMap* dilepHistos_;

  // histogram axis ranges
  double histMaxRt_;
  double histMaxM_;
  double histMaxE_;
  double histMaxPt_;

  // histogram bin numbers
  int histNBinsRt_;

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

  // HLT paths triggered on in a particular event
  std::vector<std::string> triggers_;

  // selection cuts
  double leptonPtCut_;
  double leptonEtaCut_;
  double leptonIsolationCut_;
  bool leptonChargeCut_;
  double vetoBackToBack_;
  double minDeltaRBetweenLeptons_;
  double vertexChi2Cut_;
  double deltaPhiCut_;
  double decayLengthCut_;
  double decayLengthSignificanceCut_;
  double minD0sig_;
  int maxHitsBeforeVertex_;
  int numTrigMatches_;
  int maxNumStandAloneMuons_;
  int minNumCaloMatches_;

  // cut analysis (separately for each analysis trigger)
  std::map<std::string,CutSummary*> cutSummaryLeptonSignal_;
  std::map<std::string,CutSummary*> cutSummaryLeptonBackground_;
  std::map<std::string,CutSummary*> cutSummarySignal_;
  std::map<std::string,CutSummary*> cutSummaryPartial_;
  std::map<std::string,CutSummary*> cutSummaryBackground_;

  // trigger matches storage
  std::vector<const reco::Particle::LorentzVector*> leptonTriggerMatch_;

  // calo supercluster storage
  std::vector<double> leptonCaloMatch_;

  // trigger summary storage
  unsigned numProcessedEvents_;
  std::map<std::string,unsigned> trigFilterSummary_;
  std::map<std::string,unsigned> trigPathSummary_;

  // luminosity reweighting
  edm::LumiReWeighting* LumiWeights_;
  reweight::PoissonMeanShifter PShiftUpUp_;
  reweight::PoissonMeanShifter PShiftUp_;
  reweight::PoissonMeanShifter PShiftDown_;

  // storage for lots of correction weights
  std::map<std::string,double> weightMap_;

  // quick fix for tracking efficiency studies
  double smallestLeptonRadius_;
};
