// -*- C++ -*-
//
// Package:    TrackingEfficiencyFromCosmics
// Class:      TrackingEfficiencyFromCosmics
// 
/**\class TrackingEfficiencyFromCosmics TrackingEfficiencyFromCosmics.cc Analysis/TrackingEfficiencyFromCosmics/src/TrackingEfficiencyFromCosmics.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,40 3-B32,+41227671551,
//         Created:  Wed May 25 16:44:02 CEST 2011
// $Id: TrackingEfficiencyFromCosmics.cc,v 1.34 2011/07/26 15:30:19 demattia Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

// #include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
// #include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
// #include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/Error.h"
// #include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
// #include "DataFormats/PatCandidates/interface/TriggerEvent.h"
// #include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
// #include "SimTracker/Records/interface/TrackAssociatorRecord.h"
// #include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "Analysis/TrackingEfficiencyFromCosmics/interface/AssociatorByDeltaR.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlPlots.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlDeltaPlots.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/Efficiency.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/EfficiencyTree.h"

#include <boost/foreach.hpp>

//
// class declaration
//

class TrackingEfficiencyFromCosmics : public edm::EDAnalyzer {
public:
  explicit TrackingEfficiencyFromCosmics(const edm::ParameterSet&);
  ~TrackingEfficiencyFromCosmics();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  const reco::GenParticle * takeStableMuon(const reco::GenParticleCollection & genParticles);
  // void impactParameterForGen(const reco::GenParticle & genMuon, const math::XYZPoint & genVertex,
  //                            const int genCharge, const MagneticField * mf);
  template <class T>
  void computeGenImpactParameters(const T & genMuon, const math::XYZPoint & genVertex,
                                  const int genCharge, const MagneticField * mf);
  void computeImpactParameters( const reco::Track & track, const TransientTrackBuilder & theBuilder );
  template <class T1, class T2>
  void fillEfficiencyVsGen( const T1 & staMuons, const T2 * stableMuon,
                            Efficiency * efficiency,
                            const double & genDxy, const double & genDz, const MagneticField * mf,
                            TH1F * hMinToGenDeltaR, TH1F * hToGenDeltaDxy, TH1F * hToGenDeltaDz );
  template <class T1, class T2>
  void fillEfficiency(const T1 & staMuons, const T2 & tracks, Efficiency * efficiency,
                      TH1F * hMinDeltaR, ControlDeltaPlots * standAloneTrackDelta,
		      ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone);

  void dumpGenParticleInfo(const reco::GenParticle & genParticle);
  void dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber);

  // ----------member data ---------------------------
  // double maxDeltaR_;
  TH1F * hMinDeltaR_, * hMinCleanedDeltaR_, * hSimMinDeltaR_;
  TH1F * hMinTrackToGenDeltaR_, * hMinStaMuonToGenDeltaR_, * hMinCleanedStaMuonToGenDeltaR_;
  TH1F * hStandAloneToGenDeltaDxy_, * hStandAloneToGenDeltaDz_;
  TH1F * hCleanedStandAloneToGenDeltaDxy_, * hCleanedStandAloneToGenDeltaDz_;
  TH1F * hTrackToGenDeltaDxy_, * hTrackToGenDeltaDz_;

  TH1F * hStandAloneCounterDxy_;
  TH1F * hTrackCounterDxy_;
  TH1F * hStandAloneCounterDz_;
  TH1F * hTrackCounterDz_;

  reco::Track::TrackQuality quality_;
  bool useMCtruth_;

  std::auto_ptr<AssociatorByDeltaR> associatorByDeltaR_;
  std::auto_ptr<AssociatorByDeltaR> simAssociatorByDeltaR_;
  std::auto_ptr<ControlPlots> controlPlotsGeneralTracks_;
  std::auto_ptr<ControlPlots> controlPlotsStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsMatchedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsUnmatchedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsCleanedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsCleanedStandAloneMuonsNoDzCut_;
  std::auto_ptr<ControlPlots> controlPlotsMatchedCleanedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsUnmatchedCleanedStandAloneMuons_;
  std::auto_ptr<ControlDeltaPlots> trackDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneTrackDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneTrackDelta_;
  std::auto_ptr<Efficiency> genToStandAloneEfficiency_;
  std::auto_ptr<Efficiency> genToCleanedStandAloneEfficiency_;
  std::auto_ptr<Efficiency> genToTrackEfficiency_;
  std::auto_ptr<Efficiency> efficiency_;
  std::auto_ptr<Efficiency> efficiencyCleaned_;
  boost::shared_array<double> variables_;
  unsigned int nBins_;
  double effDxyMin_, effDxyMax_;
  double effDzMin_, effDzMax_;
  double effPtMin_, effPtMax_;
  std::string effOutputFileName_;
  std::string effCleanedOutputFileName_;
  std::string genToStandAloneEffOutputFileName_;
  std::string genToTrackEffOutputFileName_;
  AnalyticalImpactPointExtrapolator * analyticalExtrapolator_;
  TransverseImpactPointExtrapolator * transverseExtrapolator_;
  AlgebraicSymMatrix55 nullCovariance_;
  std::pair<double, double> dxy_;
  std::pair<double, double> dz_;
  std::pair<double, double> dxyz_;
  int minimumValidHits_;
  double dzMinCut_;
  double dzMaxCut_;
  double dxyCut_;
  double chi2Cut_;
  double trackPtCut_;
  double standAlonePtCut_;
  bool highPurity_;
  bool matchTwoLegs_;
  double deltaDxyCut_;
  double deltaDzCut_;
  double deltaPtCut_;
  double deltaPhiCut_;
  bool recomputeIP_;
  bool singleLegMuon_;
  edm::InputTag muonCollection_;
  edm::InputTag trackCollection_;
  edm::ESHandle<TransientTrackBuilder> theB_;
  bool useAllTracks_;
  bool useTrackParameters_;
  bool dxyErrorCut_;
  bool dzErrorCut_;
  bool cleaned_;
  unsigned int eventNum_;
  double dxyCutForNoDzCut_;
};

TrackingEfficiencyFromCosmics::TrackingEfficiencyFromCosmics(const edm::ParameterSet& iConfig) :
  useMCtruth_(iConfig.getParameter<bool>("UseMCtruth")),
  effDxyMin_(iConfig.getParameter<double>("EffDxyMin")),
  effDxyMax_(iConfig.getParameter<double>("EffDxyMax")),
  effDzMin_(iConfig.getParameter<double>("EffDzMin")),
  effDzMax_(iConfig.getParameter<double>("EffDzMax")),
  effPtMin_(iConfig.getParameter<double>("EffPtMin")),
  effPtMax_(iConfig.getParameter<double>("EffPtMax")),
  effOutputFileName_(iConfig.getParameter<std::string>("EffOutputFileName")),
  effCleanedOutputFileName_(iConfig.getParameter<std::string>("EffCleanedOutputFileName")),
  genToStandAloneEffOutputFileName_(iConfig.getParameter<std::string>("GenToStandAloneEffOutputFileName")),
  genToTrackEffOutputFileName_(iConfig.getParameter<std::string>("GenToTrackEffOutputFileName")),
  minimumValidHits_(iConfig.getParameter<int>("MinimumValidHits")),
  dzMinCut_(iConfig.getParameter<double>("DzMinCut")),
  dzMaxCut_(iConfig.getParameter<double>("DzMaxCut")),
  dxyCut_(iConfig.getParameter<double>("DxyCut")),
  chi2Cut_(iConfig.getParameter<double>("Chi2Cut")),
  trackPtCut_(iConfig.getParameter<double>("TrackPtCut")),
  standAlonePtCut_(iConfig.getParameter<double>("StandAlonePtCut")),
  highPurity_(iConfig.getParameter<bool>("HighPurity")),
  matchTwoLegs_(iConfig.getParameter<bool>("MatchTwoLegs")),
  deltaDxyCut_(iConfig.getParameter<double>("DeltaDxyCut")),
  deltaDzCut_(iConfig.getParameter<double>("DeltaDzCut")),
  deltaPtCut_(iConfig.getParameter<double>("DeltaPtCut")),
  deltaPhiCut_(iConfig.getParameter<double>("DeltaPhiCut")),
  recomputeIP_(iConfig.getParameter<bool>("RecomputeIP")),
  singleLegMuon_(iConfig.getParameter<bool>("SingleLegMuon")),
  muonCollection_(iConfig.getParameter<edm::InputTag>("MuonCollection")),
  trackCollection_(iConfig.getParameter<edm::InputTag>("TrackCollection")),
  useAllTracks_(iConfig.getParameter<bool>("UseAllTracks")),
  useTrackParameters_(iConfig.getParameter<bool>("UseTrackParameters")),
  dxyErrorCut_(iConfig.getParameter<bool>("DxyErrorCut")),
  dzErrorCut_(iConfig.getParameter<bool>("DzErrorCut")),
  cleaned_(true),
  eventNum_(0),
  dxyCutForNoDzCut_(iConfig.getParameter<double>("DxyCutForNoDzCut"))
{
  // Use the unique association for tracks to standAlone muons only when
  if( singleLegMuon_ ) associatorByDeltaR_.reset(new AssociatorByDeltaR(iConfig.getParameter<double>("MaxDeltaR"), false, true));
  else associatorByDeltaR_.reset(new AssociatorByDeltaR(iConfig.getParameter<double>("MaxDeltaR")));
  simAssociatorByDeltaR_.reset(new AssociatorByDeltaR(iConfig.getParameter<double>("SimMaxDeltaR"), false));

  nBins_ = 100;

  // Build the object to compute the efficiency
  std::vector<Efficiency::Parameters> pars;
  pars.push_back(Efficiency::Parameters("dxy", nBins_, effDxyMin_, effDxyMax_));
  pars.push_back(Efficiency::Parameters("dz", nBins_, effDzMin_, effDzMax_));
  pars.push_back(Efficiency::Parameters("pt", nBins_, effPtMin_, effPtMax_));
  genToStandAloneEfficiency_.reset(new Efficiency(pars));
  genToCleanedStandAloneEfficiency_.reset(new Efficiency(pars));
  genToTrackEfficiency_.reset(new Efficiency(pars));
  efficiency_.reset(new Efficiency(pars));
  efficiencyCleaned_.reset(new Efficiency(pars));
  variables_.reset(new double[3]);

  // maxDeltaR_ = iConfig.getParameter<double>("MaxDeltaR");

  // Initialize the nullCovariance matrix
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<5; ++j ) {
      nullCovariance_(i,j) = 0;
    }
  }
}

TrackingEfficiencyFromCosmics::~TrackingEfficiencyFromCosmics() {}

void TrackingEfficiencyFromCosmics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventNum_ = iEvent.id().event();
  edm::ESHandle<MagneticField> theMF;
  iSetup.get<IdealMagneticFieldRecord>().get(theMF);
  MagneticField * mf = const_cast<MagneticField*>(&*theMF);
  transverseExtrapolator_ = new TransverseImpactPointExtrapolator(mf);
  analyticalExtrapolator_ = new AnalyticalImpactPointExtrapolator(mf);

  edm::Handle<reco::TrackCollection> allTracks;
  iEvent.getByLabel(trackCollection_, allTracks);
  reco::TrackBase::TrackQuality trackQualityHighPurity = reco::TrackBase::qualityByName("highPurity");
  reco::TrackCollection * tracks = new reco::TrackCollection();
  reco::TrackCollection::const_iterator itTrk = allTracks->begin();
  for( ; itTrk != allTracks->end(); ++itTrk ) {
    // if( (itTrk->quality(trackQualityHighPurity)) && (fabs(itTrk->eta()) < 2.0) && (itTrk->pt() > trackPtCut_) && (itTrk->found() > 6) ) {
    if( useAllTracks_ ) {
      tracks->push_back(*itTrk);
    }
    else if( (!highPurity_ || (itTrk->quality(trackQualityHighPurity))) && (fabs(itTrk->eta()) < 2.0) && (itTrk->pt() > trackPtCut_) && (itTrk->found() > 6) ) {
      tracks->push_back(*itTrk);
    }
  }
  // Load the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB_);

  edm::Handle<reco::TrackCollection> staMuons;
  iEvent.getByLabel(muonCollection_, staMuons);

  // Select good standAloneMuons (done such that the matching with gen is good)
  reco::TrackCollection cleanedStaMuons;
  reco::TrackCollection cleanedStaMuonsNoDzCut;
  reco::TrackCollection::const_iterator it = staMuons->begin();
  for( ; it != staMuons->end(); ++it ) {
    if( (it->found() >= minimumValidHits_) && (fabs(it->dxy()) < dxyCut_) &&
	(it->pt() > standAlonePtCut_) && (fabs(it->eta()) < 2.) && (fabs(it->normalizedChi2()) < chi2Cut_) &&
	((!dxyErrorCut_) || (fabs(it->dxyError()) < utils::dxyErrMax(it->pt()))) &&
	((!dzErrorCut_) || (fabs(it->dzError()) < utils::dxyErrMax(it->pt()))) ) { // Note the use of the same function is intentional.


      if( fabs(it->dxy()) < dxyCutForNoDzCut_ ) {
        cleanedStaMuonsNoDzCut.push_back(*it);
      }

      if( (fabs(it->dz()) < dzMaxCut_) && (fabs(it->dz()) > dzMinCut_) ) {
        cleanedStaMuons.push_back(*it);
      }
    }
  }

  if( singleLegMuon_ ) {
    // Two legged muon, require only one in the event
    if( cleanedStaMuons.size() != 1 ) {
      cleanedStaMuons.clear();
    }
  }
  else {
    // Single legs, require two in the event
    if( cleanedStaMuons.size() == 2 ) {
      // Compare the dxy and dz of the two muons. Accept them only if they match within 1 sigma of the resolution taken by interpolating the MC-truth
      // resolution on SingleMuPt10 and SingleMuPt100 to a Pt of 25 GeV.
      if( matchTwoLegs_ ) {
        double dxy1 = cleanedStaMuons[0].dxy();
        double dxy2 = cleanedStaMuons[1].dxy();
        double dz1 = cleanedStaMuons[0].dz();
        double dz2 = cleanedStaMuons[1].dz();
        double deltaPhi = reco::deltaPhi(cleanedStaMuons[0].phi(), cleanedStaMuons[1].phi());
        double deltaPt = cleanedStaMuons[0].pt() - cleanedStaMuons[1].pt();

        // The dxy have opposite sign while the dz have same sign.
        if( !( (fabs(dxy1 + dxy2) < deltaDxyCut_) &&
               (fabs(dz1 - dz2)   < deltaDzCut_ ) &&
               (fabs(deltaPhi)    > deltaPhiCut_) &&
               (fabs(deltaPt)     < deltaPtCut_ ) ) ) {
          cleanedStaMuons.clear();
        }
      }
    }
    else {
      cleanedStaMuons.clear();
    }
  }

  // Fill all control plots
  controlPlotsGeneralTracks_->fillControlPlots(*tracks);
  controlPlotsStandAloneMuons_->fillControlPlots(*staMuons);
  controlPlotsCleanedStandAloneMuonsNoDzCut_->fillControlPlots(cleanedStaMuonsNoDzCut);
  controlPlotsCleanedStandAloneMuons_->fillControlPlots(cleanedStaMuons);

  // // Find the simTracks (for IP)
  // edm::Handle<edm::SimTrackContainer> simTracks;
  // iEvent.getByLabel("g4SimHits", simTracks);
  // std::map<const math::XYZTLorentzVectorD *, const reco::Track *> simMatchesMap;
  // simAssociatorByDeltaR_->fillAssociationMap(*simTracks, *tracks, simMatchesMap, hSimMinDeltaR_);
  // Compute efficiency from MC truth
  // std::map<const math::XYZTLorentzVectorD *, const reco::Track *>::const_iterator it = simMatchesMap.begin();
  // for( ; it != simMatchesMap.end(); ++it ) {
  //      variables_[0] = it->first->pt();
  //      bool found = false;
  //      if( it->second != 0 ) found = true;
  //      genEfficiency_->fill(variables_, found);
  //    }

  // Gen Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    const reco::GenParticle * stableMuon = takeStableMuon(*genParticles);

    // Compute impact parameters for generator particle
    computeGenImpactParameters(*stableMuon, stableMuon->vertex(), stableMuon->charge(), mf);

    // Gen muon values
    double genDxy = dxy_.first;
    double genDz = dz_.first;

    variables_[0] = fabs(genDxy);
    variables_[1] = fabs(genDz);
    variables_[2] = stableMuon->pt();

    // Compute efficiency for track vs MC-truth
    fillEfficiencyVsGen( *tracks, stableMuon, genToTrackEfficiency_.get(), genDxy, genDz, mf,
                         hMinTrackToGenDeltaR_, hTrackToGenDeltaDxy_, hTrackToGenDeltaDz_ );
    // Compute efficiency for standalone vs MC-truth
    fillEfficiencyVsGen( *(staMuons.product()), stableMuon, genToStandAloneEfficiency_.get(), genDxy, genDz, mf,
                         hMinStaMuonToGenDeltaR_, hStandAloneToGenDeltaDxy_, hStandAloneToGenDeltaDz_ );
    // Compute efficiency for cleaned standalone vs MC-truth
    fillEfficiencyVsGen( cleanedStaMuons, stableMuon, genToCleanedStandAloneEfficiency_.get(), genDxy, genDz, mf,
                         hMinCleanedStaMuonToGenDeltaR_, hCleanedStandAloneToGenDeltaDxy_, hCleanedStandAloneToGenDeltaDz_ );
  }

  // Deltas between the two standAloneMuons
  if( !singleLegMuon_ ) {
    if( staMuons->size() == 2 ) {
      standAloneDelta_->fillControlPlots((*staMuons)[0], (*staMuons)[1]);
    }
    if( cleanedStaMuons.size() == 2 ) {
      cleanedStandAloneDelta_->fillControlPlots(cleanedStaMuons[0], cleanedStaMuons[1]);
    }
  }
  if( tracks->size() == 2 ) {
    trackDelta_->fillControlPlots((*tracks)[0], (*tracks)[1]);
  }

  // Efficiency for tracks vs standAlone
  cleaned_ = false;
  fillEfficiency(*staMuons, tracks, efficiency_.get(), hMinDeltaR_, standAloneTrackDelta_.get(),
		 controlPlotsMatchedStandAloneMuons_.get(), controlPlotsUnmatchedStandAloneMuons_.get());
  // Efficiency for tracks vs cleanedStandAlone
  cleaned_ = true;
  fillEfficiency(cleanedStaMuons, tracks, efficiencyCleaned_.get(), hMinCleanedDeltaR_, cleanedStandAloneTrackDelta_.get(),
		 controlPlotsMatchedCleanedStandAloneMuons_.get(), controlPlotsUnmatchedCleanedStandAloneMuons_.get());

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

const reco::GenParticle * TrackingEfficiencyFromCosmics::takeStableMuon(const reco::GenParticleCollection & genParticles)
{
  for( reco::GenParticleCollection::const_iterator it = genParticles.begin();
       it != genParticles.end(); ++it ) {
    // dumpGenParticleInfo(*it);
    if( it->status() == 1 ) return( &*it );
  }
  return 0;
}

// ------------ method called once each job just before starting event loop  ------------
void TrackingEfficiencyFromCosmics::beginJob()
{
  edm::Service<TFileService> fileService;
  hMinDeltaR_ =                    utils::bookHistogram(fileService, "minDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinCleanedDeltaR_ =             utils::bookHistogram(fileService, "minCleanedDeltaR", "", "#Delta R", "", 500, 0, 5);
  hSimMinDeltaR_ =                 utils::bookHistogram(fileService, "simMinDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinTrackToGenDeltaR_ =          utils::bookHistogram(fileService, "minTrackToGenDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinStaMuonToGenDeltaR_ =        utils::bookHistogram(fileService, "minStaMuonToGenDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinCleanedStaMuonToGenDeltaR_ = utils::bookHistogram(fileService, "minCleanedStaMuonToGenDeltaR", "", "#Delta R", "", 500, 0, 5);

  hStandAloneToGenDeltaDxy_ =        utils::bookHistogram(fileService, "standAloneToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
  hStandAloneToGenDeltaDz_ =         utils::bookHistogram(fileService, "standAloneToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);
  hCleanedStandAloneToGenDeltaDxy_ = utils::bookHistogram(fileService, "cleanedStandAloneToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
  hCleanedStandAloneToGenDeltaDz_ =  utils::bookHistogram(fileService, "cleanedStandAloneToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);
  hTrackToGenDeltaDxy_ =             utils::bookHistogram(fileService, "trackToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
  hTrackToGenDeltaDz_ =              utils::bookHistogram(fileService, "trackToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);

  // hTrackCounterDxy_ =      utils::bookHistogram(fileService, "trackCounterDxy", "", "|d_{0}|", "cm", 100, 0, 100);
  // hStandAloneCounterDxy_ = utils::bookHistogram(fileService, "standAloneCounterDxy", "", "|d_{0}|", "cm", 100, 0, 100);
  // hTrackCounterDz_ =       utils::bookHistogram(fileService, "trackCounterDz", "", "|d_{z}|", "cm", 100, 0, 100);
  // hStandAloneCounterDz_ =  utils::bookHistogram(fileService, "standAloneCounterDz", "", "|d_{z}|", "cm", 100, 0, 100);

  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "generalTracks"));

  controlPlotsStandAloneMuons_.reset(new ControlPlots(fileService, "standAloneMuons"));
  controlPlotsMatchedStandAloneMuons_.reset(new ControlPlots(fileService, "matchedStandAloneMuons"));
  controlPlotsUnmatchedStandAloneMuons_.reset(new ControlPlots(fileService, "unmatchedStandAloneMuons"));

  controlPlotsCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "cleanedStandAloneMuons"));
  controlPlotsCleanedStandAloneMuonsNoDzCut_.reset(new ControlPlots(fileService, "cleanedStandAloneMuonsNoDzCut"));
  controlPlotsMatchedCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "matchedCleanedStandAloneMuons"));
  controlPlotsUnmatchedCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "unmatchedCleanedStandAloneMuons"));

  trackDelta_.reset(new ControlDeltaPlots(fileService, "tracksDelta", -1));
  standAloneDelta_.reset(new ControlDeltaPlots(fileService, "standAloneDelta", -1));
  cleanedStandAloneDelta_.reset(new ControlDeltaPlots(fileService, "cleanedStandAloneDelta", -1));
  standAloneTrackDelta_.reset(new ControlDeltaPlots(fileService, "standAloneTrackDelta"));
  cleanedStandAloneTrackDelta_.reset(new ControlDeltaPlots(fileService, "cleanedStandAloneTrackDelta"));
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackingEfficiencyFromCosmics::endJob() 
{
//  TCanvas canvasDxy;
//  canvasDxy.Draw();
//  double * eff = new double[nBins_];
//  double * x = new double[nBins_];
//  hStandAloneCounterDxy_->Rebin(4);
//  hTrackCounterDxy_->Rebin(4);
//  for( unsigned int i=0; i<nBins_/4; ++i ) {
//    double numSA = hStandAloneCounterDxy_->GetBinContent(i+1);
//    double numTk = hTrackCounterDxy_->GetBinContent(i+1);
//    if(numSA == 0) eff[i] = 0;
//    else eff[i] = numTk/numSA;
//    x[i] = hStandAloneCounterDxy_->GetBinLowEdge(i);
//  }
//  TGraph * effGraphDxy = new TGraph(nBins_/4, x, eff);
//  effGraphDxy->Draw("AP");
//  canvasDxy.SaveAs("checkEffDxy.root");


//  TCanvas canvasDz;
//  canvasDz.Draw();

//  double * effDz = new double[nBins_];
//  double * xDz = new double[nBins_];
//  hStandAloneCounterDz_->Rebin(4);
//  hTrackCounterDz_->Rebin(4);
//  for( unsigned int i=0; i<nBins_/4; ++i ) {
//    double numSA = hStandAloneCounterDz_->GetBinContent(i+1);
//    double numTk = hTrackCounterDz_->GetBinContent(i+1);
//    std::cout << "numSA = " << hStandAloneCounterDz_->GetBinContent(i+1) << std::endl;
//    std::cout << "numTk = " << hTrackCounterDz_->GetBinContent(i+1) << std::endl;
//    if(numSA == 0) effDz[i] = 0;
//    else effDz[i] = numTk/numSA;
//    xDz[i] = hStandAloneCounterDz_->GetBinLowEdge(i);
//  }
//  TGraph * effGraphDz = new TGraph(nBins_/4, xDz, effDz);
//  effGraphDz->Draw("AP");
//  canvasDz.SaveAs("checkEffDz.root");




  EfficiencyTree tree;
  tree.writeTree(effOutputFileName_, &*efficiency_);

  EfficiencyTree treeCleaned;
  treeCleaned.writeTree(effCleanedOutputFileName_, &*efficiencyCleaned_);

  EfficiencyTree genToStandAloneTree;
  genToStandAloneTree.writeTree(genToStandAloneEffOutputFileName_, &*genToStandAloneEfficiency_);

  EfficiencyTree genToTrackTree;
  genToTrackTree.writeTree(genToTrackEffOutputFileName_, &*genToTrackEfficiency_);
}

// ------------ method called when starting to processes a run  ------------
void TrackingEfficiencyFromCosmics::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TrackingEfficiencyFromCosmics::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TrackingEfficiencyFromCosmics::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TrackingEfficiencyFromCosmics::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackingEfficiencyFromCosmics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template <class T>
void TrackingEfficiencyFromCosmics::computeGenImpactParameters( const T & track, const math::XYZPoint & genVertex,
                                                                const int genCharge, const MagneticField * mf )
{
  VertexDistanceXY distXY;
  VertexDistance3D dist3D;

  TVector3 genMomentum(0,0,0);
  genMomentum.SetPtEtaPhi(track.pt(),track.eta(),track.phi());
  FreeTrajectoryState ftsAtProduction(GlobalPoint(genVertex.x(),genVertex.y(),genVertex.z()),
                                      GlobalVector(genMomentum.x(),genMomentum.y(),genMomentum.z()),
                                      TrackCharge(genCharge), mf);

  TrajectoryStateOnSurface transverseTSOS_ = transverseExtrapolator_->extrapolate(ftsAtProduction, GlobalPoint(0,0,0));
  TrajectoryStateOnSurface analyticalTSOS_ = analyticalExtrapolator_->extrapolate(ftsAtProduction, GlobalPoint(0,0,0));

  if(transverseTSOS_.isValid()) {
    TrajectoryStateOnSurface transverseTSOS(transverseTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_),
                                            transverseTSOS_.surface(), transverseTSOS_.magneticField(), transverseTSOS_.weight());
    // Reco vertex default constructed to (0,0,0)
    std::pair<bool,Measurement1D> dxy = IPTools::absoluteImpactParameter(transverseTSOS, reco::Vertex(), distXY);
    dxy_.first = dxy.second.value();
    dxy_.second = dxy.second.error();
    dz_.first = transverseTSOS.globalPosition().z();
    dz_.second = transverseTSOS.cartesianError().position().czz();
  }
  else {
    std::cout << "Invalid trajectoryStateClosestToPoint for GEN" << std::endl;
    dxy_.first = 65535;
    dxy_.second = 65535;
    dz_.first = 65535;
    dz_.second = 65535;
  }
  if(analyticalTSOS_.isValid()) {
    // analyticalTSOS_ has no errors defined. Explicitly set the errors to 0 for the genparticle state
    TrajectoryStateOnSurface analyticalTSOS(analyticalTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_), analyticalTSOS_.surface(), analyticalTSOS_.magneticField()
                                            , analyticalTSOS_.weight());
    std::pair<bool,Measurement1D> dxyz = IPTools::absoluteImpactParameter(analyticalTSOS, reco::Vertex(), dist3D);
    dxyz_.first = dxyz.second.value();
    dxyz_.second = dxyz.second.error();
  }
  else {
    dxyz_.first = 65535;
    dxyz_.second = 65535;
  }
}

void TrackingEfficiencyFromCosmics::computeImpactParameters( const reco::Track & track, const TransientTrackBuilder & theBuilder )
{
  const reco::TransientTrack transientTrack = theBuilder.build(&track);
  GlobalPoint vert(0., 0., 0.);
  TrajectoryStateClosestToPoint traj = transientTrack.trajectoryStateClosestToPoint(vert);
  if( traj.isValid() ) {
    dxy_.first = traj.perigeeParameters().transverseImpactParameter();
    dxy_.second = traj.perigeeError().transverseImpactParameterError();
    dz_.first = traj.perigeeParameters().longitudinalImpactParameter();
    dz_.second = traj.perigeeError().longitudinalImpactParameterError();
    // std::cout << "From origin dxy = " << dxy_.first << " +/ " << dxy_.second << std::endl;
  }
  else {
    std::cout << "Invalid trajectoryStateClosestToPoint" << std::endl;
    dxyz_.first = 65535;
    dxyz_.second = 65535;
    dz_.first = 65535;
    dz_.second = 65535;
  }
  //  // Taking the dxy from the beamline does not make any difference
  //  TrajectoryStateClosestToBeamLine traj2 = transientTrack.stateAtBeamLine();
  //  Measurement1D measDxy = traj2.transverseImpactParameter();
  // if( traj2.isValid() ) {
  //  dxy_.first = measDxy.value();
  //  dxy_.second = measDxy.error();
  //  // std::cout << "From beamline dxy = " << dxy_.first << " +/ " << dxy_.second << std::endl;
  // }
}

void TrackingEfficiencyFromCosmics::dumpGenParticleInfo(const reco::GenParticle & genParticle)
{
  std::cout << "pdgid = " << genParticle.pdgId() << std::endl;
  std::cout << "status = " << genParticle.status() << std::endl;
  std::cout << "pt = " << genParticle.pt() << ", eta = " << genParticle.eta() << ", phi = " << genParticle.phi() << std::endl;
}

void TrackingEfficiencyFromCosmics::dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber)
{
  track.quality(quality_);
  std::stringstream ss;
  ss << "track["<<trackNumber<<"]" << std::endl <<
        "algoName = " << track.algoName() << std::endl <<
        "pt = " << track.pt() << std::endl <<
        "eta = " << track.eta() << std::endl <<
        "phi = " << track.phi() << std::endl <<
        "number of hits = " << track.recHitsSize() << std::endl <<
        "number of valid hits = " << track.found() << std::endl <<
        "number of invalid hits = " << track.lost() << std::endl <<
        "quality = " << track.qualityName(quality_) << std::endl <<
        "d0 = " << track.d0() << std::endl <<
        "d0/d0Error = " << track.d0Error() << std::endl <<
        "dz = " << track.dz() << std::endl <<
        "dz/dzError = " << track.dzError() << std::endl <<
        "chi^2/ndof = " << track.normalizedChi2() << std::endl;
  edm::LogInfo("Demo") << ss.str();
}

template <class T1, class T2>
void TrackingEfficiencyFromCosmics::fillEfficiencyVsGen( const T1 & staMuons, const T2 * stableMuon,
                                                         Efficiency * efficiency,
                                                         const double & genDxy, const double & genDz, const MagneticField * mf,
                                                         TH1F * hMinToGenDeltaR, TH1F * hToGenDeltaDxy, TH1F * hToGenDeltaDz )
{
  if( staMuons.size() == 0 ) {
    // Do it twice as we expect two standalone muons and we reconstruct none
    efficiency->fill(variables_, false);
    if( !singleLegMuon_ ) efficiency->fill(variables_, false);
  }
  else if( staMuons.size() == 1 ) {
    if( !singleLegMuon_ ) efficiency->fill(variables_, false);
  }
  else if( staMuons.size() > 2 ) {
    std::cout << "How did we get three standAloneMuons in simulation from a single cosmic track?" << std::endl;
  }
  BOOST_FOREACH( const reco::Track & staMuon, staMuons ) {
    double standAloneDxy = staMuon.dxy();
    double standAloneDz = staMuon.dz();
    if( recomputeIP_ ) {
      computeImpactParameters(staMuon, *theB_);
      standAloneDxy = dxy_.first;
      standAloneDxy = dz_.first;
    }
    hMinToGenDeltaR->Fill(reco::deltaR(*stableMuon, staMuon));
    hToGenDeltaDxy->Fill(fabs(standAloneDxy) - fabs(genDxy));
    hToGenDeltaDz->Fill(fabs(standAloneDz) - fabs(genDz));
    efficiency->fill(variables_, true);
  }
}

template <class T1, class T2>
void TrackingEfficiencyFromCosmics::fillEfficiency(const T1 & staMuons, const T2 & tracks, Efficiency * efficiency,
                                                   TH1F * hMinDeltaR, ControlDeltaPlots * standAloneTrackDelta,
						   ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone)
{
  // Association map of StandAloneMuons and TrackerTracks
  std::map<const reco::Track *, const reco::Track *> matchesMap;
  std::map<const reco::Track *, const reco::Track *> oppositeMatchesMap;
  if( staMuons.size() > 0 ) {
    associatorByDeltaR_->fillAssociationMap(staMuons, *tracks, matchesMap, hMinDeltaR, &oppositeMatchesMap);

    bool found = false;
    std::map<const reco::Track *, const reco::Track *>::const_iterator it = matchesMap.begin();
    for( ; it != matchesMap.end(); ++it ) {

      double standAloneDxy = it->first->dxy();
      double standAloneDz = it->first->dz();
      if( recomputeIP_ ) {
	computeImpactParameters(*(it->first), *theB_);
	// std::cout << "standAloneDxy            = " << standAloneDxy << std::endl;
	standAloneDxy = dxy_.first;
	// std::cout << "recomputed standAloneDxy = " << standAloneDxy << std::endl;
	standAloneDz = dz_.first;
      }
      variables_[0] = fabs(standAloneDxy);
      variables_[1] = fabs(standAloneDz);
      variables_[2] = it->first->pt();

      if( it->second == 0 ) {
        found = false;
	if( cleaned_ ) {
	  std::cout << "No match found in this event: " << eventNum_ << std::endl;
	}
      }
      else {
        found = true;
	if( useTrackParameters_ ) {
	  variables_[0] = it->second->dxy();
	  variables_[1] = it->second->dz();
	  variables_[2] = it->second->pt();
	}
      }
      efficiency->fill(variables_, found);

      // hStandAloneCounterDxy_->Fill(variables_[0]);
      // hStandAloneCounterDz_->Fill(variables_[1]);
      if( found ) {
        // hTrackCounterDxy_->Fill(variables_[0]);
        // hTrackCounterDz_->Fill(variables_[1]);
        standAloneTrackDelta->fillControlPlots(*(it->first), *(it->second));
	matchedStandAlone->fillControlPlots(staMuons);
      }
      else {
	unmatchedStandAlone->fillControlPlots(staMuons);
      }
    }

    found = false;
    std::map<const reco::Track *, const reco::Track *>::const_iterator opIt = oppositeMatchesMap.begin();
    for( ; opIt != oppositeMatchesMap.end(); ++opIt ) {

      double standAloneDxy = opIt->first->dxy();
      double standAloneDz = opIt->first->dz();
      if( recomputeIP_ ) {
        computeImpactParameters(*(opIt->first), *theB_);
        standAloneDxy = dxy_.first;
        standAloneDz = dz_.first;
      }
      variables_[0] = fabs(standAloneDxy);
      variables_[1] = fabs(standAloneDz);
      variables_[2] = opIt->first->pt();

      if( opIt->second == 0 ) {
        found = false;
	if( cleaned_ ) {
	  std::cout << "No opposite match found in this event: " << eventNum_ << std::endl;
	}
      }
      else {
        found = true;
	if( useTrackParameters_ ) {
	  variables_[0] = opIt->second->dxy();
	  variables_[1] = opIt->second->dz();
	  variables_[2] = opIt->second->pt();
	}
      }
      efficiency->fill(variables_, found);

      if( found ) {
        standAloneTrackDelta->fillControlPlots(*(opIt->first), *(opIt->second));
	matchedStandAlone->fillControlPlots(staMuons);
      }
      else {
	unmatchedStandAlone->fillControlPlots(staMuons);
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingEfficiencyFromCosmics);
