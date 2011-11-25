
// -*- C++ -*-
//
// Package:    MuonAnalyzerTreeWriter
// Class:      MuonAnalyzerTreeWriter
//
/**\class MuonAnalyzerTreeWriter MuonAnalyzerTreeWriter.cc Analysis/MuonAnalyzerTreeWriter/src/MuonAnalyzerTreeWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco De Mattia,40 3-B32,+41227671551,
//         Created:  Wed May 25 16:44:02 CEST 2011
// $Id: MuonAnalyzerTreeWriter.cc,v 1.42 2011/11/21 16:58:00 demattia Exp $
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

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "Analysis/TrackingEfficiencyFromCosmics/interface/AssociatorByDeltaR.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlPlots.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/ControlDeltaPlots.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/Efficiency.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/EfficiencyTree.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"
#include "Analysis/Records/interface/SmartPropagatorWithIPComponentsRecord.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/RootTreeHandler.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/TreeTrack.h"

#include <boost/foreach.hpp>

//
// class declaration
//

class MuonAnalyzerTreeWriter : public edm::EDAnalyzer {
public:
  explicit MuonAnalyzerTreeWriter(const edm::ParameterSet&);
  ~MuonAnalyzerTreeWriter();

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
  template <class T>
  SmartPropagatorWithIP::IP computeGenImpactParameters(const T & genMuon, const math::XYZPoint & genVertex,
                                                       const int genCharge);
  // void computeImpactParameters( const reco::Track & track, const TransientTrackBuilder & theBuilder );
  template <class T1, class T2>
  void fillEfficiencyVsGen( const T1 & tracks, const T2 stableMuon,
                            Efficiency * efficiency,
                            // const double & genDxy, const double & genDz, const MagneticField * mf,
                            const GlobalPoint & vertex,
                            TH1F * hMinToGenDeltaR);
  template <class T1, class T2>
  void fillEfficiency(const T1 & staMuons, const T2 & tracks, Efficiency * efficiency,
                      TH1F * hMinDeltaR, ControlDeltaPlots * standAloneTrackDelta,
                      ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone);

  void fillEfficiency(const std::map<const reco::Track *, const reco::Track *> & matchesMap,
                      Efficiency * efficiency, ControlDeltaPlots * standAloneTrackDelta,
                      ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone);

  void dumpGenParticleInfo(const reco::GenParticle & genParticle);
  void dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber);
  template <class T1, class T2>
  void cleanMuons(const T1 & startingCollection, T2 & cleanedCollection);
  template <class T>
  void makeExclusive( T & collection );
  template <class T>
  void fillTreeTracks( const T & collection, const GlobalPoint & vertex, std::vector<TreeTrack> & tracks );
  void setGen(std::vector<TreeTrack> & treeTracks, const reco::GenParticle * gen, const SmartPropagatorWithIP::IP & genIP);

  // ----------member data ---------------------------
  TH1F * hMinDeltaR_, * hMinCleanedDeltaR_, * hSimMinDeltaR_;
  TH1F * hMinTrackToGenDeltaR_, * hMinStaMuonToGenDeltaR_, * hMinCleanedStaMuonToGenDeltaR_;

  TH1F * hStandAloneCounterDxy_;
  TH1F * hTrackCounterDxy_;
  TH1F * hStandAloneCounterDz_;
  TH1F * hTrackCounterDz_;

  reco::Track::TrackQuality quality_;
  bool useMCtruth_;

  std::auto_ptr<AssociatorByDeltaR> associatorByDeltaR_;
  std::auto_ptr<AssociatorByDeltaR> simAssociatorByDeltaR_;

  std::auto_ptr<ControlPlots> controlPlotsGenTracks_;
  std::auto_ptr<ControlPlots> controlPlotsGeneralTracks_;
  std::auto_ptr<ControlPlots> controlPlotsStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsMatchedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsUnmatchedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsCleanedStandAloneMuons_;
  // std::auto_ptr<ControlPlots> controlPlotsCleanedStandAloneMuonsNoDzCut_;
  std::auto_ptr<ControlPlots> controlPlotsMatchedCleanedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsUnmatchedCleanedStandAloneMuons_;

  std::auto_ptr<ControlDeltaPlots> trackDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneTrackDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneTrackDelta_;

  std::auto_ptr<ControlDeltaPlots> trackVsGenDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneVsGenDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneVsGenDelta_;

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
  const SmartPropagatorWithIP * smartPropIP_;
  bool useAllTracks_;
  bool useTrackParameters_;
  bool dxyErrorCut_;
  bool dzErrorCut_;
  bool cleaned_;
  bool countSameSide_;
  bool countOppoSide_;
  unsigned int eventNum_;
  bool phiRegion_;
  double phiMin_;
  double phiMax_;
  bool genInsideTkVol_;
  boost::shared_ptr<RootTreeHandler> treeHandlerStandAloneMuons_;
  boost::shared_ptr<RootTreeHandler> treeHandlerCleanedStandAloneMuons_;
};

MuonAnalyzerTreeWriter::MuonAnalyzerTreeWriter(const edm::ParameterSet& iConfig) :
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
  countSameSide_(iConfig.getParameter<bool>("CountSameSide")),
  countOppoSide_(iConfig.getParameter<bool>("CountOppoSide")),
  eventNum_(0),
  phiRegion_(iConfig.getParameter<bool>("PhiRegion")),
  phiMin_(iConfig.getParameter<double>("PhiMinCut")),
  phiMax_(iConfig.getParameter<double>("PhiMaxCut")),
  genInsideTkVol_(iConfig.getParameter<bool>("GenInsideTkVol"))
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

  // Initialize the nullCovariance matrix
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<5; ++j ) {
      nullCovariance_(i,j) = 0;
    }
  }
  smartPropIP_ = 0;

  treeHandlerStandAloneMuons_.reset(new RootTreeHandler(muonCollection_.label()+".root"));
  treeHandlerCleanedStandAloneMuons_.reset(new RootTreeHandler("cleaned"+muonCollection_.label()+".root"));
}

MuonAnalyzerTreeWriter::~MuonAnalyzerTreeWriter()
{
  std::cout << "Saving trees" << std::endl;
  treeHandlerStandAloneMuons_->writeTree();
  treeHandlerCleanedStandAloneMuons_->writeTree();
}

void MuonAnalyzerTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventNum_ = iEvent.id().event();
  std::cout << "event = " << eventNum_ << std::endl;

  // Vertex with respect to which do all the propagations
  GlobalPoint vertex(0,0,0);

  // Load the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB_);

  // Load the propagator and IP calculator
  edm::ESHandle<Propagator> smartPropIPHandle;
  iSetup.get<SmartPropagatorWithIPComponentsRecord>().get("SmartPropagatorWithIP", smartPropIPHandle);
  smartPropIP_ = dynamic_cast<const SmartPropagatorWithIP*>(&*smartPropIPHandle);

  reco::TrackCollection * tracks = new reco::TrackCollection();

  try {
    edm::Handle<reco::TrackCollection> allTracks;
    iEvent.getByLabel(trackCollection_, allTracks);
    reco::TrackBase::TrackQuality trackQualityHighPurity = reco::TrackBase::qualityByName("highPurity");
    reco::TrackCollection::const_iterator itTrk = allTracks->begin();
    for( ; itTrk != allTracks->end(); ++itTrk ) {
      if( useAllTracks_ ) {
        tracks->push_back(*itTrk);
      }
      else if( (!highPurity_ || (itTrk->quality(trackQualityHighPurity))) &&
               (!phiRegion_ || (itTrk->phi() > phiMin_ && itTrk->phi() < phiMax_ )) &&
               (fabs(itTrk->eta()) < 2.0) && (itTrk->pt() > trackPtCut_) && (itTrk->found() > 6) ) {
        tracks->push_back(*itTrk);
      }
    }
  }
  catch (cms::Exception & ex) {
    std::cerr << ex;
  }

  edm::Handle<reco::TrackCollection> staMuons;
  iEvent.getByLabel(muonCollection_, staMuons);

  // Fill the track collection for the tree
  std::vector<TreeTrack> standAloneMuonsTreeTracks;
  fillTreeTracks(*(staMuons.product()), vertex, standAloneMuonsTreeTracks);

  // Select good standAloneMuons (done such that the matching with gen is good)
  reco::TrackCollection cleanedStaMuons;
  cleanMuons(staMuons, cleanedStaMuons);
  // keep the muons only if there is exactly 1 (single leg) or 2 (separated legs) muons in the collection.
  // In the case of two legs also ask for matching if requested in the cfg.
  makeExclusive(cleanedStaMuons);

  // Fill the cleaned track collection for the tree
  std::vector<TreeTrack> cleanedStandAloneMuonsTreeTracks;
  fillTreeTracks(cleanedStaMuons, vertex, cleanedStandAloneMuonsTreeTracks);

  // Fill all control plots
  controlPlotsGeneralTracks_->fillControlPlots(*tracks, smartPropIP_);
  controlPlotsStandAloneMuons_->fillControlPlots(*staMuons, smartPropIP_);
  controlPlotsCleanedStandAloneMuons_->fillControlPlots(cleanedStaMuons, smartPropIP_);

  SmartPropagatorWithIP::IP tkIp1, tkIp2;
  SmartPropagatorWithIP::IP saIp1, saIp2;
  SmartPropagatorWithIP::IP cleanedSaIp1, cleanedSaIp2;
  if( tracks->size() > 0 ) {
    tkIp1 = smartPropIP_->computeImpactParameters((*tracks)[0], vertex);
  }
  if( tracks->size() > 1 ) {
    tkIp2 = smartPropIP_->computeImpactParameters((*tracks)[1], vertex);
  }
  if( staMuons->size() > 0 ) {
    saIp1 = smartPropIP_->computeImpactParameters((*staMuons)[0], vertex);
  }
  if( staMuons->size() > 1 ) {
    saIp2 = smartPropIP_->computeImpactParameters((*staMuons)[1], vertex);
  }
  if( cleanedStaMuons.size() > 0 ) {
    cleanedSaIp1 = smartPropIP_->computeImpactParameters(cleanedStaMuons[0], vertex);
  }
  if( cleanedStaMuons.size() > 1 ) {
    cleanedSaIp2 = smartPropIP_->computeImpactParameters(cleanedStaMuons[1], vertex);
  }

  // Gen Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    const reco::GenParticle * stableMuon = takeStableMuon(*genParticles);

    // Compute impact parameters for generator particle
    SmartPropagatorWithIP::IP stableMuonIP(computeGenImpactParameters(*stableMuon, stableMuon->vertex(), stableMuon->charge()));

    // Assume one gen particle per event. This is true when using cosmics or particleGun.
    setGen(standAloneMuonsTreeTracks, stableMuon, stableMuonIP);
    setGen(cleanedStandAloneMuonsTreeTracks, stableMuon, stableMuonIP);

    // Gen muon values
    variables_[0] = fabs(stableMuonIP.dxyValue);
    variables_[1] = fabs(stableMuonIP.dzValue);
    variables_[2] = stableMuonIP.pt;

    controlPlotsGenTracks_->fillControlPlots(stableMuonIP, stableMuon->vertex());

    // Compute efficiency for track vs MC-truth
    fillEfficiencyVsGen( *tracks, stableMuonIP, genToTrackEfficiency_.get(), vertex,
                         hMinTrackToGenDeltaR_ );
    // Compute efficiency for standalone vs MC-truth
    fillEfficiencyVsGen( *(staMuons.product()), stableMuonIP, genToStandAloneEfficiency_.get(), vertex,
                         hMinStaMuonToGenDeltaR_ );
    // Compute efficiency for cleaned standalone vs MC-truth
    fillEfficiencyVsGen( cleanedStaMuons, stableMuonIP, genToCleanedStandAloneEfficiency_.get(), vertex,
                         hMinCleanedStaMuonToGenDeltaR_ );

    if( tracks->size() > 0 ) {
      trackVsGenDelta_->fillControlPlots((*tracks)[0], tkIp1, *stableMuon, stableMuonIP);
    }
    if( tracks->size() > 1 ) {
      trackVsGenDelta_->fillControlPlots((*tracks)[1], tkIp2, *stableMuon, stableMuonIP);
    }
    if( staMuons->size() > 0 ) {
      standAloneVsGenDelta_->fillControlPlots((*staMuons)[0], saIp1, *stableMuon, stableMuonIP);
    }
    if( staMuons->size() > 1 ) {
      standAloneVsGenDelta_->fillControlPlots((*staMuons)[1], saIp2, *stableMuon, stableMuonIP);
    }
    if( cleanedStaMuons.size() > 0 ) {
      cleanedStandAloneVsGenDelta_->fillControlPlots(cleanedStaMuons[0], cleanedSaIp1, *stableMuon, stableMuonIP);
    }
    if( cleanedStaMuons.size() > 1 ) {
      cleanedStandAloneVsGenDelta_->fillControlPlots(cleanedStaMuons[1], cleanedSaIp2, *stableMuon, stableMuonIP);
    }
  }

  // Fill the trees
  treeHandlerStandAloneMuons_->saveToTree(standAloneMuonsTreeTracks, iEvent.id().event(), iEvent.run());
  treeHandlerCleanedStandAloneMuons_->saveToTree(cleanedStandAloneMuonsTreeTracks, iEvent.id().event(), iEvent.id().run());

  // Deltas between the two standAloneMuons
  if( !singleLegMuon_ ) {
    if( staMuons->size() == 2 ) {
      standAloneDelta_->fillControlPlots((*staMuons)[0], saIp1, (*staMuons)[1], saIp2);
    }
    if( cleanedStaMuons.size() == 2 ) {
      cleanedStandAloneDelta_->fillControlPlots(cleanedStaMuons[0], saIp1, cleanedStaMuons[1], saIp2);
    }
  }
  if( tracks->size() == 2 ) {
    trackDelta_->fillControlPlots((*tracks)[0], tkIp1, (*tracks)[1], tkIp2);
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

const reco::GenParticle * MuonAnalyzerTreeWriter::takeStableMuon(const reco::GenParticleCollection & genParticles)
{
  for( reco::GenParticleCollection::const_iterator it = genParticles.begin();
       it != genParticles.end(); ++it ) {
    // dumpGenParticleInfo(*it);
    if( it->status() == 1 ) return( &*it );
  }
  return 0;
}

// ------------ method called once each job just before starting event loop  ------------
void MuonAnalyzerTreeWriter::beginJob()
{
  edm::Service<TFileService> fileService;
  hMinDeltaR_ =                    utils::bookHistogram(fileService, "minDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinCleanedDeltaR_ =             utils::bookHistogram(fileService, "minCleanedDeltaR", "", "#Delta R", "", 500, 0, 5);
  hSimMinDeltaR_ =                 utils::bookHistogram(fileService, "simMinDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinTrackToGenDeltaR_ =          utils::bookHistogram(fileService, "minTrackToGenDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinStaMuonToGenDeltaR_ =        utils::bookHistogram(fileService, "minStaMuonToGenDeltaR", "", "#Delta R", "", 500, 0, 5);
  hMinCleanedStaMuonToGenDeltaR_ = utils::bookHistogram(fileService, "minCleanedStaMuonToGenDeltaR", "", "#Delta R", "", 500, 0, 5);

  controlPlotsGenTracks_.reset(new ControlPlots(fileService, "genTracks"));
  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "generalTracks"));
  //  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "hltL2Muons"));
  controlPlotsStandAloneMuons_.reset(new ControlPlots(fileService, "standAloneMuons"));
  controlPlotsMatchedStandAloneMuons_.reset(new ControlPlots(fileService, "matchedStandAloneMuons"));
  controlPlotsUnmatchedStandAloneMuons_.reset(new ControlPlots(fileService, "unmatchedStandAloneMuons"));

  controlPlotsCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "cleanedStandAloneMuons"));
  // controlPlotsCleanedStandAloneMuonsNoDzCut_.reset(new ControlPlots(fileService, "cleanedStandAloneMuonsNoDzCut"));
  controlPlotsMatchedCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "matchedCleanedStandAloneMuons"));
  controlPlotsUnmatchedCleanedStandAloneMuons_.reset(new ControlPlots(fileService, "unmatchedCleanedStandAloneMuons"));

  trackDelta_.reset(new ControlDeltaPlots(fileService, "tracksDelta", -1));
  standAloneDelta_.reset(new ControlDeltaPlots(fileService, "standAloneDelta", -1));
  cleanedStandAloneDelta_.reset(new ControlDeltaPlots(fileService, "cleanedStandAloneDelta", -1));
  standAloneTrackDelta_.reset(new ControlDeltaPlots(fileService, "standAloneTrackDelta"));
  cleanedStandAloneTrackDelta_.reset(new ControlDeltaPlots(fileService, "cleanedStandAloneTrackDelta"));

  trackVsGenDelta_.reset(new ControlDeltaPlots(fileService, "trackVsGenDelta", -1));
  standAloneVsGenDelta_.reset(new ControlDeltaPlots(fileService, "standAloneVsGenDelta", -1));
  cleanedStandAloneVsGenDelta_.reset(new ControlDeltaPlots(fileService, "cleanedStandAloneVsGenDelta", -1));
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonAnalyzerTreeWriter::endJob()
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
void MuonAnalyzerTreeWriter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void MuonAnalyzerTreeWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void MuonAnalyzerTreeWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuonAnalyzerTreeWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonAnalyzerTreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template <class T>
SmartPropagatorWithIP::IP MuonAnalyzerTreeWriter::computeGenImpactParameters( const T & track, const math::XYZPoint & genVertex,
                                                                                     const int genCharge )
{
  SmartPropagatorWithIP::IP ip;
  if( genInsideTkVol_ ) {
    ip = smartPropIP_->computeGenImpactParametersInsideTkVol( track, genVertex, genCharge, GlobalPoint(0,0,0) );
  }
  else {
    ip = smartPropIP_->computeGenImpactParametersOutsideTkVol( track, genVertex, genCharge, GlobalPoint(0,0,0) );
  }
  // std::cout << "gen pt before propagation = " << track.pt() << std::endl;
  // std::cout << "gen pt after propagation = " << ip.pt << std::endl;
  // std::cout << "gen dxy = " << ip.dxyValue << std::endl;
  return ip;
}

void MuonAnalyzerTreeWriter::dumpGenParticleInfo(const reco::GenParticle & genParticle)
{
  std::cout << "pdgid = " << genParticle.pdgId() << std::endl;
  std::cout << "status = " << genParticle.status() << std::endl;
  std::cout << "pt = " << genParticle.pt() << ", eta = " << genParticle.eta() << ", phi = " << genParticle.phi() << std::endl;
}

void MuonAnalyzerTreeWriter::dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber)
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
void MuonAnalyzerTreeWriter::fillEfficiencyVsGen( const T1 & tracks, const T2 stableMuonIP,
                                                  Efficiency * efficiency,
                                                  const GlobalPoint & vertex,
                                                  TH1F * hMinToGenDeltaR)
{
  if( tracks.size() == 0 ) {
    // Do it twice as we expect two standalone muons and we reconstruct none
    efficiency->fill(variables_, false);
    if( !singleLegMuon_ ) efficiency->fill(variables_, false);
  }
  else if( tracks.size() == 1 ) {
    if( !singleLegMuon_ ) efficiency->fill(variables_, false);
  }
  else if( tracks.size() > 2 ) {
    std::cout << "How did we get more than two standAloneMuons in simulation from a single cosmic track?" << std::endl;
  }
  BOOST_FOREACH( const reco::Track & track, tracks ) {
    SmartPropagatorWithIP::IP ip(track.pt(), track.ptError(),
                                 track.eta(), track.etaError(),
                                 track.phi(), track.phiError(),
                                 track.dxy(), track.dxyError(),
                                 track.dz(), track.dzError());
    if( recomputeIP_ ) {
      ip = smartPropIP_->computeImpactParameters(track, vertex);
    }
    hMinToGenDeltaR->Fill(reco::deltaR(stableMuonIP.eta, stableMuonIP.phi, ip.eta, ip.phi));
    //    hMinToGenDeltaR->Fill(std::min(reco::deltaR(stableMuonIP.eta, stableMuonIP.phi, ip.eta, ip.phi),
    //                                   reco::deltaR(stableMuonIP.eta, stableMuonIP.phi, -ip.eta, ip.phi+TMath::Pi()/2.)));
    efficiency->fill(variables_, true);
  }
}

void MuonAnalyzerTreeWriter::fillEfficiency(const std::map<const reco::Track *, const reco::Track *> & matchesMap,
                                                   Efficiency * efficiency, ControlDeltaPlots * standAloneTrackDelta,
                                                   ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone)
{
  bool found = false;
  std::map<const reco::Track *, const reco::Track *>::const_iterator it = matchesMap.begin();
  for( ; it != matchesMap.end(); ++it ) {
    found = false;
    SmartPropagatorWithIP::IP ip1(it->first->pt(), it->first->ptError(),
                                  it->first->eta(), it->first->etaError(),
                                  it->first->phi(), it->first->phiError(),
                                  it->first->dxy(), it->first->dxyError(),
                                  it->first->dz(), it->first->dzError());
    if( recomputeIP_ ) {
      // ip1 = smartPropIP_->computeImpactParametersOutsideInTkVol(*(it->first), GlobalPoint(0,0,0));
      ip1 = smartPropIP_->computeImpactParameters(*(it->first), GlobalPoint(0,0,0));
    }
    SmartPropagatorWithIP::IP ip2;
    variables_[0] = fabs(ip1.dxyValue);
    variables_[1] = fabs(ip1.dzValue);
    variables_[2] = ip1.pt;

    if( it->second == 0 ) {
      found = false;
      if( cleaned_ ) {
        std::cout << "No match found in this event: " << eventNum_ << std::endl;
      }
    }
    else {
      found = true;
      if( recomputeIP_ ) {
        // ip2 = smartPropIP_->computeImpactParametersInsideTkVol(*(it->second), GlobalPoint(0,0,0));
        ip2 = smartPropIP_->computeImpactParameters(*(it->second), GlobalPoint(0,0,0));
      }
      else {
        ip2 = SmartPropagatorWithIP::IP(it->second->pt(), it->second->ptError(),
                                        it->second->eta(), it->second->etaError(),
                                        it->second->phi(), it->second->phiError(),
                                        it->second->dxy(), it->second->dxyError(),
                                        it->second->dz(), it->second->dzError());
      }
      if( useTrackParameters_ ) {
        variables_[0] = ip2.dxyValue;
        variables_[1] = ip2.dzValue;
        variables_[2] = ip2.pt;
      }
    }
    efficiency->fill(variables_, found);

    if( found ) {
      standAloneTrackDelta->fillControlPlots(*(it->first), ip1, *(it->second), ip2);
      matchedStandAlone->fillControlPlots(*(it->first), smartPropIP_);
    }
    else {
      unmatchedStandAlone->fillControlPlots(*(it->first), smartPropIP_);
    }
  }
}

template <class T1, class T2>
void MuonAnalyzerTreeWriter::fillEfficiency(const T1 & staMuons, const T2 & tracks, Efficiency * efficiency,
                                                   TH1F * hMinDeltaR, ControlDeltaPlots * standAloneTrackDelta,
                                                   ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone)
{
  // Association map of StandAloneMuons and TrackerTracks
  std::map<const reco::Track *, const reco::Track *> matchesMap;
  std::map<const reco::Track *, const reco::Track *> oppositeMatchesMap;
  if( staMuons.size() > 0 ) {
    associatorByDeltaR_->fillAssociationMap(staMuons, *tracks, matchesMap, hMinDeltaR, &oppositeMatchesMap);
    if( countSameSide_ ) {
      fillEfficiency(matchesMap, efficiency, standAloneTrackDelta, matchedStandAlone, unmatchedStandAlone);
    }
    if( countOppoSide_ ) {
      fillEfficiency(oppositeMatchesMap, efficiency, standAloneTrackDelta, matchedStandAlone, unmatchedStandAlone);
    }
  }
}

template <class T1, class T2>
void MuonAnalyzerTreeWriter::cleanMuons( const T1 & startingCollection, T2 & cleanedCollection)
{
  reco::TrackCollection::const_iterator it = startingCollection->begin();
  for( ; it != startingCollection->end(); ++it ) {
    if( (it->found() >= minimumValidHits_) &&
        ( it->hitPattern().dtStationsWithValidHits() + it->hitPattern().cscStationsWithValidHits() > 1 ) &&
        (it->pt() > standAlonePtCut_) && (fabs(it->eta()) < 2.) && (fabs(it->normalizedChi2()) < chi2Cut_) &&
        ((!dxyErrorCut_) || (fabs(it->dxyError()) < utils::dxyErrMax(it->pt()))) &&
        ((!dzErrorCut_) || (fabs(it->dzError()) < utils::dxyErrMax(it->pt()))) ) { // Note the use of the same function is intentional.
      if( (fabs(it->dz()) < dzMaxCut_) && (fabs(it->dz()) > dzMinCut_) && (fabs(it->dxy()) < dxyCut_) ) {
        cleanedCollection.push_back(*it);
      }
    }
  }
}

template <class T>
void MuonAnalyzerTreeWriter::makeExclusive( T & collection )
{
  if( singleLegMuon_ ) {
    // Two legged muon, require only one in the event
    if( collection.size() != 1 ) {
      collection.clear();
    }
  }
  else {
    // Single legs, require two in the event
    if( collection.size() == 2 ) {
      // Compare the dxy and dz of the two muons. Accept them only if they match within 1 sigma of the resolution taken by interpolating the MC-truth
      // resolution on SingleMuPt10 and SingleMuPt100 to a Pt of 25 GeV.
      if( matchTwoLegs_ ) {
        // The dxy have opposite sign while the dz have same sign.
        if( !( (fabs(collection[0].dxy() + collection[1].dxy()) < deltaDxyCut_) &&
               (fabs(collection[0].dz() - collection[1].dz()) < deltaDzCut_ ) &&
               (fabs(reco::deltaPhi(collection[0].phi(), collection[1].phi())) > deltaPhiCut_) &&
               (fabs(collection[0].pt() - collection[1].pt()) < deltaPtCut_ ) ) ) {
          collection.clear();
        }
      }
    }
    else {
      collection.clear();
    }
  }
}

template <class T>
void MuonAnalyzerTreeWriter::fillTreeTracks( const T & collection, const GlobalPoint & vertex, std::vector<TreeTrack> & tracks )
{
  typename T::const_iterator it = collection.begin();
  for( ; it != collection.end(); ++it ) {
    SmartPropagatorWithIP::IP ip(it->pt(), it->ptError(), it->eta(), it->etaError(), it->phi(), it->phiError(),
                                 it->dxy(), it->dxyError(), it->dz(), it->dzError());
    TreeTrack treeTrack;
    if( recomputeIP_ ) {
      ip = smartPropIP_->computeImpactParameters(*it, vertex);
      utils::fillTrackToTreeTrack(treeTrack, *it, ip);
      tracks.push_back(treeTrack);
      // tracks.push_back(TreeTrack(*it, ip));
    }
    else {
      utils::fillTrackToTreeTrack(treeTrack, *it);
      tracks.push_back(treeTrack);
      // tracks.push_back(TreeTrack(*it));
    }
  }
}

void MuonAnalyzerTreeWriter::setGen(std::vector<TreeTrack> & treeTracks, const reco::GenParticle * gen, const SmartPropagatorWithIP::IP & genIP)
{
  std::vector<TreeTrack>::iterator it = treeTracks.begin();
  for( ; it != treeTracks.end(); ++it ) {
    utils::fillGenToTreeTrack(*it, *gen, genIP);
    // it->setGen(*gen, gen->vertex(), genIP);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzerTreeWriter);
