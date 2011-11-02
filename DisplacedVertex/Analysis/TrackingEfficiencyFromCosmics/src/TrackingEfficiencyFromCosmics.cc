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
// $Id: TrackingEfficiencyFromCosmics.cc,v 1.42 2011/08/03 15:42:09 demattia Exp $
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
  SmartPropagatorWithIP::IP computeGenImpactParameters(const T & genMuon, const math::XYZPoint & genVertex,
                                                       const int genCharge, const MagneticField * mf);
  // void computeImpactParameters( const reco::Track & track, const TransientTrackBuilder & theBuilder );
  template <class T1, class T2>
  void fillEfficiencyVsGen( const T1 & tracks, const T2 stableMuon,
                            Efficiency * efficiency,
                            // const double & genDxy, const double & genDz, const MagneticField * mf,
                            TH1F * hMinToGenDeltaR,
                            const ControlPlots::propType propagationType = ControlPlots::OUTSIDEIN );
  template <class T1, class T2>
  void fillEfficiency(const T1 & staMuons, const T2 & tracks, Efficiency * efficiency,
                      TH1F * hMinDeltaR, ControlDeltaPlots * standAloneTrackDelta,
		      ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone);

  void fillEfficiency(const std::map<const reco::Track *, const reco::Track *> & matchesMap,
                      Efficiency * efficiency, ControlDeltaPlots * standAloneTrackDelta,
                      ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone);

  void dumpGenParticleInfo(const reco::GenParticle & genParticle);
  void dumpTrackInfo(const reco::Track & track, const unsigned int trackNumber);

  // ----------member data ---------------------------
  // double maxDeltaR_;
  TH1F * hMinDeltaR_, * hMinCleanedDeltaR_, * hSimMinDeltaR_;
  TH1F * hMinTrackToGenDeltaR_, * hMinStaMuonToGenDeltaR_, * hMinCleanedStaMuonToGenDeltaR_;
//  TH1F * hStandAloneToGenDeltaDxy_, * hStandAloneToGenDeltaDz_;
//  TH1F * hCleanedStandAloneToGenDeltaDxy_, * hCleanedStandAloneToGenDeltaDz_;
//  TH1F * hTrackToGenDeltaDxy_, * hTrackToGenDeltaDz_;

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
  std::auto_ptr<ControlPlots> controlPlotsCleanedStandAloneMuonsNoDzCut_;
  std::auto_ptr<ControlPlots> controlPlotsMatchedCleanedStandAloneMuons_;
  std::auto_ptr<ControlPlots> controlPlotsUnmatchedCleanedStandAloneMuons_;

  std::auto_ptr<ControlDeltaPlots> trackDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneTrackDelta_;
  std::auto_ptr<ControlDeltaPlots> cleanedStandAloneTrackDelta_;

  std::auto_ptr<ControlDeltaPlots> trackVsGenDelta_;
  std::auto_ptr<ControlDeltaPlots> standAloneVsGenDelta_;

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
  const SmartPropagatorWithIP * smartPropIP_;
  bool useAllTracks_;
  bool useTrackParameters_;
  bool dxyErrorCut_;
  bool dzErrorCut_;
  bool cleaned_;
  bool countSameSide_;
  bool countOppoSide_;
  unsigned int eventNum_;
  double dxyCutForNoDzCut_;
  bool phiRegion_;
  double phiMin_;
  double phiMax_;
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
  countSameSide_(iConfig.getParameter<bool>("CountSameSide")),
  countOppoSide_(iConfig.getParameter<bool>("CountOppoSide")),
  eventNum_(0),
  dxyCutForNoDzCut_(iConfig.getParameter<double>("DxyCutForNoDzCut")),
  phiRegion_(iConfig.getParameter<bool>("PhiRegion")),
  phiMin_(iConfig.getParameter<double>("PhiMinCut")),
  phiMax_(iConfig.getParameter<double>("PhiMaxCut"))
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
  smartPropIP_ = 0;
}

TrackingEfficiencyFromCosmics::~TrackingEfficiencyFromCosmics() {}

void TrackingEfficiencyFromCosmics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventNum_ = iEvent.id().event();

  // Load the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB_);

  // Load the propagator and IP calculator
  edm::ESHandle<Propagator> smartPropIPHandle;
  iSetup.get<SmartPropagatorWithIPComponentsRecord>().get("SmartPropagatorWithIP", smartPropIPHandle);
  smartPropIP_ = dynamic_cast<const SmartPropagatorWithIP*>(&*smartPropIPHandle);

  // edm::ESHandle<MagneticField> theMF;
  // iSetup.get<IdealMagneticFieldRecord>().get(theMF);
  // MagneticField * mf = const_cast<MagneticField*>(&*theMF);
  std::cout << "event = " << eventNum_ << std::endl;
  const MagneticField * mf = smartPropIP_->magneticField();
  std::cout << "mf = " << mf << std::endl;
  transverseExtrapolator_ = new TransverseImpactPointExtrapolator(mf);
  analyticalExtrapolator_ = new AnalyticalImpactPointExtrapolator(mf);

  reco::TrackCollection * tracks = new reco::TrackCollection();

  try {
    edm::Handle<reco::TrackCollection> allTracks;
    // edm::Handle<edm::View<reco::Track> > allTracks;
    iEvent.getByLabel(trackCollection_, allTracks);
    reco::TrackBase::TrackQuality trackQualityHighPurity = reco::TrackBase::qualityByName("highPurity");
    reco::TrackCollection::const_iterator itTrk = allTracks->begin();
    // edm::View<reco::Track>::const_iterator itTrk = allTracks->begin();
    for( ; itTrk != allTracks->end(); ++itTrk ) {
      // if( (itTrk->quality(trackQualityHighPurity)) && (fabs(itTrk->eta()) < 2.0) && (itTrk->pt() > trackPtCut_) && (itTrk->found() > 6) ) {
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
    // std::cerr << ex.what();
    // std::cerr << ex.explainSelf();
  }

  edm::Handle<reco::TrackCollection> staMuons;
  iEvent.getByLabel(muonCollection_, staMuons);

  // Select good standAloneMuons (done such that the matching with gen is good)
  reco::TrackCollection cleanedStaMuons;
  reco::TrackCollection cleanedStaMuonsNoDzCut;
  reco::TrackCollection::const_iterator it = staMuons->begin();
  for( ; it != staMuons->end(); ++it ) {
    if( (it->found() >= minimumValidHits_) &&
	(it->pt() > standAlonePtCut_) && (fabs(it->eta()) < 2.) && (fabs(it->normalizedChi2()) < chi2Cut_) &&
	((!dxyErrorCut_) || (fabs(it->dxyError()) < utils::dxyErrMax(it->pt()))) &&
	((!dzErrorCut_) || (fabs(it->dzError()) < utils::dxyErrMax(it->pt()))) ) { // Note the use of the same function is intentional.


      if( fabs(it->dxy()) < dxyCutForNoDzCut_ ) {
        cleanedStaMuonsNoDzCut.push_back(*it);
      }

      if( (fabs(it->dz()) < dzMaxCut_) && (fabs(it->dz()) > dzMinCut_) && (fabs(it->dxy()) < dxyCut_) ) {
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
  controlPlotsGeneralTracks_->fillControlPlots(*tracks, smartPropIP_, ControlPlots::INSIDETKVOL);
  controlPlotsStandAloneMuons_->fillControlPlots(*staMuons, smartPropIP_, ControlPlots::OUTSIDEIN);
  controlPlotsCleanedStandAloneMuonsNoDzCut_->fillControlPlots(cleanedStaMuonsNoDzCut, smartPropIP_, ControlPlots::OUTSIDEIN);
  controlPlotsCleanedStandAloneMuons_->fillControlPlots(cleanedStaMuons, smartPropIP_, ControlPlots::OUTSIDEIN);

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

  GlobalPoint vertex(0,0,0);

  SmartPropagatorWithIP::IP tkIp1, tkIp2;
  SmartPropagatorWithIP::IP saIp1, saIp2;
  if( tracks->size() > 0 ) {
    tkIp1 = smartPropIP_->computeImpactParametersInsideTkVol((*tracks)[0], vertex);
  }
  if( tracks->size() > 1 ) {
    tkIp2 = smartPropIP_->computeImpactParametersInsideTkVol((*tracks)[1], vertex);
  }
  if( staMuons->size() > 0 ) {
    saIp1 = smartPropIP_->computeImpactParametersOutsideInTkVol((*staMuons)[0], vertex);
  }
  if( staMuons->size() > 1 ) {
    saIp2 = smartPropIP_->computeImpactParametersOutsideInTkVol((*staMuons)[1], vertex);
  }

  // Gen Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if( useMCtruth_ ) {
    iEvent.getByLabel("genParticles", genParticles);
    const reco::GenParticle * stableMuon = takeStableMuon(*genParticles);

    // Compute impact parameters for generator particle
    SmartPropagatorWithIP::IP stableMuonIP(computeGenImpactParameters(*stableMuon, stableMuon->vertex(), stableMuon->charge(), mf));

    // Gen muon values
    variables_[0] = fabs(stableMuonIP.dxyValue);
    variables_[1] = fabs(stableMuonIP.dzValue);
    variables_[2] = stableMuonIP.pt;

    controlPlotsGenTracks_->fillControlPlots(stableMuonIP, stableMuon->vertex());

    // Compute efficiency for track vs MC-truth
    fillEfficiencyVsGen( *tracks, stableMuonIP, genToTrackEfficiency_.get(),
                         hMinTrackToGenDeltaR_, ControlPlots::INSIDEOUT );
    // Compute efficiency for standalone vs MC-truth
    fillEfficiencyVsGen( *(staMuons.product()), stableMuonIP, genToStandAloneEfficiency_.get(),
                         hMinStaMuonToGenDeltaR_ );
    // Compute efficiency for cleaned standalone vs MC-truth
    fillEfficiencyVsGen( cleanedStaMuons, stableMuonIP, genToCleanedStandAloneEfficiency_.get(),
                         hMinCleanedStaMuonToGenDeltaR_ );

    if( tracks->size() > 0 ) {
      // trackVsGenDelta_->fillControlPlots((*tracks)[0], (*tracks)[0].dxy(), (*tracks)[0].dz(), *stableMuon, genDxy, genDz);
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
  }

  // Deltas between the two standAloneMuons
  if( !singleLegMuon_ ) {
    if( staMuons->size() == 2 ) {
      // standAloneDelta_->fillControlPlots((*staMuons)[0], (*staMuons)[1]);
      standAloneDelta_->fillControlPlots((*staMuons)[0], saIp1, (*staMuons)[1], saIp2);
    }
    if( cleanedStaMuons.size() == 2 ) {
      // cleanedStandAloneDelta_->fillControlPlots(cleanedStaMuons[0], cleanedStaMuons[1]);
      cleanedStandAloneDelta_->fillControlPlots(cleanedStaMuons[0], saIp1, cleanedStaMuons[1], saIp2);
    }
  }
  if( tracks->size() == 2 ) {
    // trackDelta_->fillControlPlots((*tracks)[0], (*tracks)[1]);
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

//  hStandAloneToGenDeltaDxy_ =        utils::bookHistogram(fileService, "standAloneToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
//  hStandAloneToGenDeltaDz_ =         utils::bookHistogram(fileService, "standAloneToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);
//  hCleanedStandAloneToGenDeltaDxy_ = utils::bookHistogram(fileService, "cleanedStandAloneToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
//  hCleanedStandAloneToGenDeltaDz_ =  utils::bookHistogram(fileService, "cleanedStandAloneToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);
//  hTrackToGenDeltaDxy_ =             utils::bookHistogram(fileService, "trackToGenDeltaDxy", "", "|#Delta |d_{0}||", "cm", 100, -100, 100);
//  hTrackToGenDeltaDz_ =              utils::bookHistogram(fileService, "trackToGenDeltaDz", "", "|#Delta |d_{z}||", "cm", 100, -100, 100);

  // hTrackCounterDxy_ =      utils::bookHistogram(fileService, "trackCounterDxy", "", "|d_{0}|", "cm", 100, 0, 100);
  // hStandAloneCounterDxy_ = utils::bookHistogram(fileService, "standAloneCounterDxy", "", "|d_{0}|", "cm", 100, 0, 100);
  // hTrackCounterDz_ =       utils::bookHistogram(fileService, "trackCounterDz", "", "|d_{z}|", "cm", 100, 0, 100);
  // hStandAloneCounterDz_ =  utils::bookHistogram(fileService, "standAloneCounterDz", "", "|d_{z}|", "cm", 100, 0, 100);

  controlPlotsGenTracks_.reset(new ControlPlots(fileService, "genTracks"));
  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "generalTracks"));
//  controlPlotsGeneralTracks_.reset(new ControlPlots(fileService, "hltL2Muons"));
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

  trackVsGenDelta_.reset(new ControlDeltaPlots(fileService, "trackVsGenDelta", -1));
  standAloneVsGenDelta_.reset(new ControlDeltaPlots(fileService, "standAloneVsGenDelta", -1));
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
SmartPropagatorWithIP::IP TrackingEfficiencyFromCosmics::computeGenImpactParameters( const T & track, const math::XYZPoint & genVertex,
                                                                const int genCharge, const MagneticField * mf )
{
  SmartPropagatorWithIP::IP ip(smartPropIP_->computeGenImpactParametersOutsideTkVol( track, genVertex, genCharge, GlobalPoint(0,0,0) ));

  VertexDistanceXY distXY;
  VertexDistance3D dist3D;

  TVector3 genMomentum(0,0,0);
  genMomentum.SetPtEtaPhi(track.pt(),track.eta(),track.phi());
  FreeTrajectoryState ftsAtProduction(GlobalPoint(genVertex.x(),genVertex.y(),genVertex.z()),
                                      GlobalVector(genMomentum.x(),genMomentum.y(),genMomentum.z()),
                                      TrackCharge(genCharge), mf);

  std::cout << "gen pt before propagation = " << track.pt() << std::endl;
  std::cout << "gen pt after propagation = " << ip.pt << std::endl;
//  std::cout << "genVertex (x,y,z) = (" << genVertex.x() << "," << genVertex.y() << "," << genVertex.z() << ")" << std::endl;
//  std::cout << "genMomentum (x,y,z) = (" << genMomentum.x() << "," << genMomentum.y() << "," << genMomentum.z() << ")" << std::endl;
//  std::cout << "genCharge = (" << genCharge <<  std::endl;

//  TrajectoryStateOnSurface analyticalTSOS_ = analyticalExtrapolator_->extrapolate(ftsAtProduction, GlobalPoint(0,0,0));
//  TrajectoryStateOnSurface transverseTSOS_ = transverseExtrapolator_->extrapolate(ftsAtProduction, GlobalPoint(0,0,0));

//  if(transverseTSOS_.isValid()) {
//    TrajectoryStateOnSurface transverseTSOS(transverseTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_),
//                                            transverseTSOS_.surface(), transverseTSOS_.magneticField(), transverseTSOS_.weight());
//    // Reco vertex default constructed to (0,0,0)
//    std::pair<bool,Measurement1D> dxy = IPTools::absoluteImpactParameter(transverseTSOS, reco::Vertex(), distXY);
//    dxy_.first = dxy.second.value();
//    dxy_.second = dxy.second.error();
//    dz_.first = transverseTSOS.globalPosition().z();
//    dz_.second = transverseTSOS.cartesianError().position().czz();
//    std::cout << "transverseTSOS pt = " << transverseTSOS.globalMomentum().perp() << std::endl;
//  }
//  else {
//    std::cout << "Invalid trajectoryStateClosestToPoint for GEN" << std::endl;
//    dxy_.first = 65535;
//    dxy_.second = 65535;
//    dz_.first = 65535;
//    dz_.second = 65535;
//  }
//  if(analyticalTSOS_.isValid()) {
//    std::cout << "analytical extrapolation successful" << std::endl;
//    // analyticalTSOS_ has no errors defined. Explicitly set the errors to 0 for the genparticle state
//    TrajectoryStateOnSurface analyticalTSOS(analyticalTSOS_.localParameters(), LocalTrajectoryError(nullCovariance_),
//                                            analyticalTSOS_.surface(), analyticalTSOS_.magneticField(), analyticalTSOS_.weight());
//    std::pair<bool,Measurement1D> dxyz = IPTools::absoluteImpactParameter(analyticalTSOS, reco::Vertex(), dist3D);
//    dxyz_.first = dxyz.second.value();
//    dxyz_.second = dxyz.second.error();
//  }
//  else {
//    std::cout << "Invalid trajectoryStateClosestToPoint for analyticalExtrapolator for GEN" << std::endl;
//    dxyz_.first = 65535;
//    dxyz_.second = 65535;
//  }

//  std::cout << "transverseExtrapolator dxy = " << dxy_.first << std::endl;
//  std::cout << "transverseExtrapolator dz = " << dz_.first << std::endl;
//  dxy_.first = ip.dxyValue;
//  dxy_.second = ip.dxyError;
//  dz_.first = ip.dzValue;
//  dz_.second = ip.dzError;
//  std::cout << "steppingHelix dxy = " << dxy_.first << std::endl;
//  std::cout << "steppingHelix dz = " << dz_.first << std::endl;

  return ip;
}

//void TrackingEfficiencyFromCosmics::computeImpactParameters( const reco::Track & track, const TransientTrackBuilder & theBuilder )
//{
//  const reco::TransientTrack transientTrack = theBuilder.build(&track);
//  GlobalPoint vert(0., 0., 0.);
//  TrajectoryStateClosestToPoint traj = transientTrack.trajectoryStateClosestToPoint(vert);
//  if( traj.isValid() ) {
//    dxy_.first = traj.perigeeParameters().transverseImpactParameter();
//    dxy_.second = traj.perigeeError().transverseImpactParameterError();
//    dz_.first = traj.perigeeParameters().longitudinalImpactParameter();
//    dz_.second = traj.perigeeError().longitudinalImpactParameterError();
//    std::cout << "From origin dxy = " << dxy_.first << " +/- " << dxy_.second << std::endl;
//  }
//  else {
//    std::cout << "Invalid trajectoryStateClosestToPoint" << std::endl;
//    dxyz_.first = 65535;
//    dxyz_.second = 65535;
//    dz_.first = 65535;
//    dz_.second = 65535;
//  }
//  //  // Taking the dxy from the beamline does not make any difference
//  //  TrajectoryStateClosestToBeamLine traj2 = transientTrack.stateAtBeamLine();
//  //  Measurement1D measDxy = traj2.transverseImpactParameter();
//  //  if( traj2.isValid() ) {
//  //    dxy_.first = measDxy.value();
//  //    dxy_.second = measDxy.error();
//  //    // std::cout << "From beamline dxy = " << dxy_.first << " +/ " << dxy_.second << std::endl;
//  //  }
//}

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
void TrackingEfficiencyFromCosmics::fillEfficiencyVsGen( const T1 & tracks, const T2 stableMuonIP,
                                                         Efficiency * efficiency,
                                                         TH1F * hMinToGenDeltaR,
                                                         const ControlPlots::propType propagationType )
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
    SmartPropagatorWithIP::IP ip(track.pt(), track.eta(), track.phi(),
                                 track.dxy(), track.dxyError(), track.dz(), track.dzError());
    if( recomputeIP_ ) {
      if( propagationType == ControlPlots::INSIDETKVOL ) {
        ip = smartPropIP_->computeImpactParametersInsideTkVol(track, GlobalPoint(0,0,0));
      }
      else if( propagationType == ControlPlots::OUTSIDEIN ) {
        ip = smartPropIP_->computeImpactParametersOutsideInTkVol(track, GlobalPoint(0,0,0));
      }
      else {
        std::cout << "[TrackingEfficiencyFromCosmics::fillEfficiencyVsGen]: unknown propagation type: " << propagationType << std::endl;
      }
    }
    hMinToGenDeltaR->Fill(reco::deltaR(stableMuonIP.eta, stableMuonIP.phi, ip.eta, ip.phi));
    efficiency->fill(variables_, true);
  }
}

void TrackingEfficiencyFromCosmics::fillEfficiency(const std::map<const reco::Track *, const reco::Track *> & matchesMap,
                                                   Efficiency * efficiency, ControlDeltaPlots * standAloneTrackDelta,
                                                   ControlPlots * matchedStandAlone, ControlPlots * unmatchedStandAlone)
{
  bool found = false;
  std::map<const reco::Track *, const reco::Track *>::const_iterator it = matchesMap.begin();
  for( ; it != matchesMap.end(); ++it ) {
    found = false;
    SmartPropagatorWithIP::IP ip1(it->first->pt(), it->first->eta(), it->first->phi(),
                                  it->first->dxy(), it->first->dxyError(), it->first->dz(), it->first->dzError());
    if( recomputeIP_ ) {
      ip1 = smartPropIP_->computeImpactParametersOutsideInTkVol(*(it->first), GlobalPoint(0,0,0));
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
        ip2 = smartPropIP_->computeImpactParametersInsideTkVol(*(it->second), GlobalPoint(0,0,0));
      }
      else {
        ip2 = SmartPropagatorWithIP::IP(it->second->pt(), it->second->eta(), it->second->phi(),
                                        it->second->dxy(), it->second->dxyError(), it->second->dz(), it->second->dzError());
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
      matchedStandAlone->fillControlPlots(*(it->first), smartPropIP_, ControlPlots::OUTSIDEIN);
    }
    else {
      unmatchedStandAlone->fillControlPlots(*(it->first), smartPropIP_, ControlPlots::OUTSIDEIN);
    }
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
    if( countSameSide_ ) {
      fillEfficiency(matchesMap, efficiency, standAloneTrackDelta, matchedStandAlone, unmatchedStandAlone);
    }
    if( countOppoSide_ ) {
      fillEfficiency(oppositeMatchesMap, efficiency, standAloneTrackDelta, matchedStandAlone, unmatchedStandAlone);
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingEfficiencyFromCosmics);
