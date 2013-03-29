#include "HeavyFlavorAnalysis/Onia2MuMu/interface/Onia2MuMuPAT.h"

//Headers for the data items
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//Headers for services and tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

Onia2MuMuPAT::Onia2MuMuPAT(const edm::ParameterSet& iConfig):
  muons_(iConfig.getParameter<edm::InputTag>("muons")),
  thebeamspot_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
  thePVs_(iConfig.getParameter<edm::InputTag>("primaryVertexTag")),
  higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
  lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
  dimuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
  addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
  addMuonlessPrimaryVertex_(iConfig.getParameter<bool>("addMuonlessPrimaryVertex")),
  resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),
  addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
  addThirdTrack_(iConfig.getParameter<bool>("addThirdTrack")),
  minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
  trackMass_(iConfig.getParameter<double>("trackMass")),
  diMuPlusTrackMassMax_(iConfig.getParameter<double>("diMuPlusTrackMassMax")),
  diMuPlusTrackMassMin_(iConfig.getParameter<double>("diMuPlusTrackMassMin")),
  diMuMassMax_(iConfig.getParameter<double>("diMuMassMax")),
  diMuMassMin_(iConfig.getParameter<double>("diMuMassMin")),
  prefilter_(iConfig.existsAs<std::string>("preselection") ? iConfig.getParameter<std::string>("preselection") : "")
{
  produces<pat::CompositeCandidateCollection>("DiMu");
  produces<pat::CompositeCandidateCollection>("DiMuTk");
  //produces<pat::CompositeCandidateCollection>();
}


Onia2MuMuPAT::~Onia2MuMuPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

void Onia2MuMuPAT::fillCandMuons(pat::CompositeCandidate & myCand, const pat::Muon & mu1, const pat::Muon & mu2,
                                 const reco::TransientTrack & refittedTrack1, const reco::TransientTrack & refittedTrack2)
{
  myCand.addDaughter(mu1, "muon1");
  int highPurity1 = 0;
  if(refittedTrack1.track().quality(reco::TrackBase::highPurity)) highPurity1 = 1;
  // if(mu1.track().quality(reco::TrackBase::highPurity)) highPurity1 = 1;
  myCand.addUserInt("highPurity1", highPurity1);
  myCand.addUserFloat("segComp1", muon::segmentCompatibility(mu1));
  // std::cout << "segComp1 = " << muon::segmentCompatibility(mu1) << std::endl;

  myCand.addDaughter(mu2, "muon2");
  int highPurity2 = 0;
  if(refittedTrack2.track().quality(reco::TrackBase::highPurity)) highPurity2 = 1;
  // if(mu2.track().quality(reco::TrackBase::highPurity)) highPurity2 = 1;
  myCand.addUserInt("highPurity2", highPurity2);
  myCand.addUserFloat("segComp2", muon::segmentCompatibility(mu2));
  // std::cout << "segComp2 = " << muon::segmentCompatibility(mu2) << std::endl;
}

math::XYZTLorentzVector Onia2MuMuPAT::fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi, const double & mass)
{
  double px = pt*cos(phi);
  double py = pt*sin(phi);
  double tmp = 2*atan(exp(-eta));
  double pz = pt*cos(tmp)/sin(tmp);
  double E  = sqrt(px*px+py*py+pz*pz+mass*mass);

  return math::XYZTLorentzVector(px,py,pz,E);
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
Onia2MuMuPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  vector<double> muMasses;
  muMasses.push_back( 0.1134289256 );
  muMasses.push_back( 0.1134289256 );

  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutput(new pat::CompositeCandidateCollection);
  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutputWithTrack(new pat::CompositeCandidateCollection);

  Vertex thePrimaryV;
  Vertex theBeamSpotV;

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  Handle< View<pat::Muon> > muons;
  iEvent.getByLabel(muons_,muons);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);
  TrackCollection muonLess;

  // dimuon candidates only from muons
  for(View<pat::Muon>::const_iterator it = muons->begin(), itend = muons->end(); it != itend; ++it){
    // both must pass low quality
    if(!lowerPuritySelection_(*it)) continue;
    for(View<pat::Muon>::const_iterator it2 = it+1; it2 != itend;++it2){
      // both must pass low quality
      if(!lowerPuritySelection_(*it2)) continue;
      // one must pass tight quality
      if (!(higherPuritySelection_(*it) || higherPuritySelection_(*it2))) continue;

      pat::CompositeCandidate myCand;
      pat::CompositeCandidate my3Cand;
      bool TrackFound=false;
      Vertex PVDiMuTk,SVDiMuTk;
      Track track;
      LeafCandidate tkcand;

      vector<TransientVertex> pvs;

      // Compute invariant mass to reject candidates before the vertex fit
      math::XYZTLorentzVector m1p4(fromPtEtaPhiToPxPyPz(it->track()->pt(),
                                                        it->track()->eta(),
                                                        it->track()->phi(),
                                                        muMasses[0]));
      math::XYZTLorentzVector m2p4(fromPtEtaPhiToPxPyPz(it2->track()->pt(),
                                                        it2->track()->eta(),
                                                        it2->track()->phi(),
                                                        muMasses[1]));
      LorentzVector dimuonBeforeRefit = m1p4 + m2p4;
      // std::cout << "mass = " << dimuonBeforeRefit.mass() << std::endl;

      // ---- apply the dimuon cut ----
      if(!dimuonSelection_(dimuonBeforeRefit)) continue;

      // ---- fit vertex using Tracker tracks (if they have tracks) ----
      if (it->track().isNonnull() && it2->track().isNonnull()) {

        if (it->track()->pt()==it2->track()->pt()) continue; //against split-muons
        double M2=myCand.mass();
        if (addThirdTrack_ && (it->charge()!=it2->charge()) && M2 > diMuMassMin_ && M2 < diMuMassMax_)
          TrackFound=searchForTheThirdTrack(iEvent, iSetup, my3Cand, &(*it), &(*it2),PVDiMuTk,track,resolveAmbiguity_);

        if (TrackFound){
          my3Cand.addUserData("PVwithmuons",PVDiMuTk);
          if (addCommonVertex_) {
            vector<TransientTrack> t_tks3;
            t_tks3.push_back(theTTBuilder->build(*it->track()));  // pass the reco::Track, not the reco::TrackRef (which can be transient)
            t_tks3.push_back(theTTBuilder->build(*it2->track())); // otherwise the vertex will have transient refs inside.
            t_tks3.push_back(theTTBuilder->build(track));
            TransientVertex tvtx=vtxFitter.vertex(t_tks3);
            SVDiMuTk=Vertex(tvtx);
            my3Cand.addUserData("commonVertex",SVDiMuTk);
          }
          LorentzVector vtk=LorentzVector(track.px(),track.py(),track.pz(),sqrt((track.p()*track.p())+(trackMass_*trackMass_)));
          tkcand=LeafCandidate(track.charge(),vtk);
          my3Cand.addDaughter(tkcand,"Kaon");
        }else{
          if (addCommonVertex_) my3Cand.addUserData("commonVertex",Vertex());
          my3Cand.addUserData("PVwithmuons",Vertex());
        }

        //build the dimuon secondary vertex
        vector<TransientTrack> t_tks;
        t_tks.push_back(theTTBuilder->build(*it->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(*it2->track())); // otherwise the vertex will have transient refs inside.
        TransientVertex myVertex = vtxFitter.vertex(t_tks);

        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
        Measurement1D MassWErr = massCalculator.invariantMass( VtxForInvMass, muMasses );

        myCand.addUserFloat("MassErr",MassWErr.error());

        if (myVertex.isValid()) {
          float vChi2 = myVertex.totalChiSquared();
          float vNDF  = myVertex.degreesOfFreedom();
          float vProb(TMath::Prob(vChi2,(int)vNDF));

          // Get refitted lepton tracks
          reco::TransientTrack refittedTrack1(myVertex.refittedTrack(t_tks[0]));
          reco::TransientTrack refittedTrack2(myVertex.refittedTrack(t_tks[1]));
          pat::Muon mu1(*it);
          math::XYZTLorentzVector m1p4(fromPtEtaPhiToPxPyPz(refittedTrack1.track().pt(),
                                                            refittedTrack1.track().eta(),
                                                            refittedTrack1.track().phi(),
                                                            muMasses[0]));
          mu1.setP4(m1p4);
          pat::Muon mu2(*it2);
          math::XYZTLorentzVector m2p4(fromPtEtaPhiToPxPyPz(refittedTrack2.track().pt(),
                                                            refittedTrack2.track().eta(),
                                                            refittedTrack2.track().phi(),
                                                            muMasses[1]));
          mu2.setP4(m2p4);

          // Make sure mu1 is the highest pt muon
          if( mu1.pt() > mu2.pt() ) {
            fillCandMuons(myCand, mu1, mu2, refittedTrack1, refittedTrack2);
            fillCandMuons(my3Cand, mu1, mu2, refittedTrack1, refittedTrack2);
          }
          else {
            fillCandMuons(myCand, mu2, mu1, refittedTrack2, refittedTrack1);
            fillCandMuons(my3Cand, mu2, mu1, refittedTrack2, refittedTrack1);
          }

          myCand.addUserFloat("vNChi2",vChi2/vNDF);
          myCand.addUserFloat("vProb",vProb);

          LorentzVector dimuon = m1p4 + m2p4;
          myCand.setP4(dimuon);
          myCand.setCharge(mu1.charge()+mu2.charge());

          TVector3 vtx;
          TVector3 pvtx;
          VertexDistanceXY vdistXY;
          TVector3 vtx3D;
          TVector3 pvtx3D;
          VertexDistance3D vdist3D;

          // XY
          vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
          TVector3 pperp(dimuon.px(), dimuon.py(), 0);
          AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
          // 3D
          vtx3D.SetXYZ(myVertex.position().x(),myVertex.position().y(),myVertex.position().z());
          TVector3 pCand(dimuon.px(), dimuon.py(), dimuon.pz());
          // AlgebraicVector3 vpCand(pCand.x(),pCand.y(),pCand.z());

          if (resolveAmbiguity_) {

            float minDz = 999999.;
            TwoTrackMinimumDistance ttmd;
            bool status = ttmd.calculate( GlobalTrajectoryParameters(GlobalPoint(myVertex.position().x(), myVertex.position().y(), myVertex.position().z()),
                                                                     GlobalVector(myCand.px(),myCand.py(),myCand.pz()),TrackCharge(0),&(*magneticField)),
                                          GlobalTrajectoryParameters(GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
                                                                     GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)) );
            float extrapZ=-9E20;
            if (status) extrapZ=ttmd.points().first.z();

            for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv){
              float deltaZ = fabs(extrapZ - itv->position().z()) ;
              if ( deltaZ < minDz ) {
                minDz = deltaZ;
                thePrimaryV = Vertex(*itv);
              }
            }
          }
          Vertex theOriginalPV = thePrimaryV;

          muonLess.clear();
          muonLess.reserve(thePrimaryV.tracksSize());
          if( addMuonlessPrimaryVertex_  && thePrimaryV.tracksSize()>2 ) {
            // Primary vertex matched to the dimuon, now refit it removing the two muons
            VertexReProducer revertex(priVtxs, iEvent);
            Handle<TrackCollection> pvtracks;
            iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
            if( !pvtracks.isValid()) { std::cout << "pvtracks NOT valid " << std::endl; }
            else {
              Handle<BeamSpot>        pvbeamspot;
              iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
              if (pvbeamspot.id() != theBeamSpot.id()) edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
              // I need to go back to the reco::Muon object, as the TrackRef in the pat::Muon can be an embedded ref.
              const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(it->originalObject());
              const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(it2->originalObject());
              // check that muons are truly from reco::Muons (and not, e.g., from PF Muons)
              // also check that the tracks really come from the track collection used for the BS
              if (rmu1 != 0 && rmu2 != 0 && rmu1->track().id() == pvtracks.id() && rmu2->track().id() == pvtracks.id()) {
                if( thePrimaryV.hasRefittedTracks() ) {
                  // Need to go back to the original tracks before taking the key
                  std::vector<reco::Track>::const_iterator itRefittedTrack = thePrimaryV.refittedTracks().begin();
                  std::vector<reco::Track>::const_iterator refittedTracksEnd = thePrimaryV.refittedTracks().end();
                  for( ; itRefittedTrack != refittedTracksEnd; ++itRefittedTrack ) {
                    if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu1->track().key() ) continue;
                    if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu2->track().key() ) continue;
                    muonLess.push_back(*(thePrimaryV.originalTrack(*itRefittedTrack)));
                  }
                }
                else {
                  std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV.tracks_begin();
                  for( ; itPVtrack != thePrimaryV.tracks_end(); ++itPVtrack ) if (itPVtrack->isNonnull()) {
                    if( itPVtrack->key() == rmu1->track().key() ) continue;
                    if( itPVtrack->key() == rmu2->track().key() ) continue;
                    muonLess.push_back(**itPVtrack);
                  }
                }
                if (muonLess.size()>1 && muonLess.size() < thePrimaryV.tracksSize()){
                  pvs = revertex.makeVertices(muonLess, *pvbeamspot, iSetup) ;
                  if (!pvs.empty()) {
                    Vertex muonLessPV = Vertex(pvs.front());
                    thePrimaryV = muonLessPV;
                  }
                }
              }
            }
          }

          // count the number of high Purity tracks with pT > 900 MeV attached to the chosen vertex
          double vertexWeight = 0., sumPTPV = 0.;
          int countTksOfPV = 0, Ntrk=0;
          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(it->originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(it2->originalObject());
          double minDcaDR0p7 = 9999.;
          double minDcaDR0p7PtMin0p5 = 9999.;
          double minDcaDR0p7PtMin0p9 = 9999.;
          double minDca = 9999.;
          double minDcaPtMin0p5 = 9999.;
          double minDcaPtMin0p9 = 9999.;
          double sumNdof = 0.;
          double PVndof = theOriginalPV.ndof();
          myCand.addUserFloat("PVndof",(float) PVndof);
          myCand.addUserFloat("pvw8",(float) (thePrimaryV.ndof()+2)/(2*thePrimaryV.tracksSize()));

          try{
            for(reco::Vertex::trackRef_iterator itVtx = theOriginalPV.tracks_begin(); itVtx != theOriginalPV.tracks_end(); itVtx++){
              if(itVtx->isNonnull()){
                if (rmu1 != 0 && rmu1->innerTrack().key() == itVtx->key())
                  continue;
                if (rmu2 != 0 && rmu2->innerTrack().key() == itVtx->key())
                  continue;
                const reco::Track& track = **itVtx;
                // if(!track.quality(reco::TrackBase::highPurity)) continue;
                TransientTrack tt = theTTBuilder->build(track);
                pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,myVertex);
                vertexWeight += theOriginalPV.trackWeight(*itVtx);
                if(track.pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;

                // Use all tracks
                if( tkPVdist.first ) {
                  if( tkPVdist.second.value()<minDca ) minDca = tkPVdist.second.value();
                  if( track.pt() > 0.5 && tkPVdist.second.value()<minDcaPtMin0p5 ) minDcaPtMin0p5 = tkPVdist.second.value();
                  if( track.pt() > 0.9 && tkPVdist.second.value()<minDcaPtMin0p9 ) minDcaPtMin0p9 = tkPVdist.second.value();
                }

                // Use only tracks close to the Bs momentum direction
                if(deltaR(track.eta(),track.phi(),myCand.eta(), myCand.phi()) < 0.7) {
                  if( tkPVdist.first ) {
                    if( tkPVdist.second.value()<minDcaDR0p7 ) minDcaDR0p7 = tkPVdist.second.value();
                    if( track.pt() > 0.5 && tkPVdist.second.value()<minDcaDR0p7PtMin0p5 ) minDcaDR0p7PtMin0p5 = tkPVdist.second.value();
                    if( track.pt() > 0.9 && tkPVdist.second.value()<minDcaDR0p7PtMin0p9 ) minDcaDR0p7PtMin0p9 = tkPVdist.second.value();
                  }
                  if(track.pt() < 0.9) continue; //reject all low pt tracks from counting
                  sumPTPV += track.pt();
                }
                sumNdof += track.ndof();
                countTksOfPV++;
              }
            }
          } catch (std::exception & err) {std::cout << " muon Selection%G�%@failed " << std::endl; return ; }

          // -----Count the number of the tracks not associated to any vertex--------
          Handle<TrackCollection> tkColl;
          iEvent.getByLabel("generalTracks", tkColl);

          double sumPt_noVtx = 0.;
          int countTksOfnoVtx = 0;
          // Find the tracks not belonging to any PV
          for (TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk) {
            bool hasVertex = false;
            for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv){
              for( reco::Vertex::trackRef_iterator itVtx = itv->tracks_begin(); itVtx != itv->tracks_end(); ++itVtx) {
                if( itVtx->isNonnull() ) {
                  if(tk->pt() == (*itVtx)->pt() && tk->eta()==(*itVtx)->eta()) {
                    hasVertex = true;
                    // Found a vertex, break out of all loops and go the identifier "gotoVertexFound"
                    goto gotoVertexFound;
                  }
                }
              }
            }
gotoVertexFound:
            if(hasVertex) continue;
            if(rmu1 != 0 && (*rmu1->innerTrack()).pt() == tk->pt() && (*rmu1->innerTrack()).eta() == tk->eta()) continue;
            if(rmu2 != 0 && (*rmu2->innerTrack()).pt() == tk->pt() && (*rmu2->innerTrack()).eta() == tk->eta()) continue;
            TransientTrack tt = theTTBuilder->build(*tk);
            // Count number of close tracks not associated to any PV
            pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,myVertex);
            if(tk->pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;

            // Use all tracks
            if( tkPVdist.first ) {
              if( tkPVdist.second.value()<minDca ) minDca = tkPVdist.second.value();
              if( tk->pt() > 0.5 && tkPVdist.second.value()<minDcaPtMin0p5 ) minDcaPtMin0p5 = tkPVdist.second.value();
              if( tk->pt() > 0.9 && tkPVdist.second.value()<minDcaPtMin0p9 ) minDcaPtMin0p9 = tkPVdist.second.value();
            }

            // Use only tracks close to the Bs momentum direction
            if(deltaR(tk->eta(),tk->phi(),myCand.eta(), myCand.phi()) < 0.7) {
              if( tkPVdist.first ) {
                if( tkPVdist.second.value()<minDcaDR0p7 ) minDcaDR0p7 = tkPVdist.second.value();
                if( tk->pt() > 0.5 && tkPVdist.second.value()<minDcaDR0p7PtMin0p5 ) minDcaDR0p7PtMin0p5 = tkPVdist.second.value();
                if( tk->pt() > 0.9 && tkPVdist.second.value()<minDcaDR0p7PtMin0p9 ) minDcaDR0p7PtMin0p9 = tkPVdist.second.value();
              }
              if(tk->pt() < 0.9) continue; //reject all low pt tracks from counting
              if(tkPVdist.second.value() > 0.05) continue;
              sumPTPV += tk->pt();
              countTksOfnoVtx++;
            }
          }
          double Iso = myCand.pt()/(myCand.pt()+sumPt_noVtx+sumPTPV);

          myCand.addUserFloat("Isolation", (float) Iso);
          myCand.addUserFloat("minDca", (float) minDca);
          myCand.addUserFloat("minDcaDr0p7", (float) minDcaDR0p7);
          myCand.addUserFloat("minDcaDr0p7Min0p5", (float) minDcaDR0p7PtMin0p5);
          myCand.addUserFloat("minDcaDr0p7Min0p9", (float) minDcaDR0p7PtMin0p9);
          myCand.addUserInt("Ntrk", Ntrk);
          myCand.addUserInt("countTksOfPV", countTksOfPV);
          myCand.addUserFloat("vertexWeight", (float) vertexWeight);
          myCand.addUserFloat("sumPTPV", (float) sumPTPV);
          myCand.addUserFloat("sumNdof", (float) sumNdof);

          // DCA
          // Compute both the dca in the transverse plane (dcaxy) and in 3D.
          TrajectoryStateClosestToPoint mu1TS = t_tks[0].impactPointTSCP();
          TrajectoryStateClosestToPoint mu2TS = t_tks[1].impactPointTSCP();
          float dcaxy = 1E20;
          float dca = 1E20;
          if (mu1TS.isValid() && mu2TS.isValid()) {
            ClosestApproachInRPhi cApp;
            if( cApp.calculate(mu1TS.theState(), mu2TS.theState()) ) dcaxy = cApp.distance();
            TwoTrackMinimumDistance ttmd;
            if( ttmd.calculate(mu1TS.theState(), mu2TS.theState()) ) dca = ttmd.distance();
          }
          myCand.addUserFloat("DCAXY", dcaxy);
          myCand.addUserFloat("DCA", dca);
          if (TrackFound){
            my3Cand.addUserFloat("DCAXY", dcaxy);
            my3Cand.addUserFloat("DCA", dca);
          }
          // end DCA

          // TODO: check that the error matrices can be used and combined as such. If they are variance-covariance matrices they have the
          // variance (and covariance terms), i.e. sigma^2, so they can be summed. Verify this.
          FreeTrajectoryState fts(
                GlobalPoint(myVertex.position().x(), myVertex.position().y(), myVertex.position().z()),
                GlobalVector(myCand.px(),myCand.py(),myCand.pz()),TrackCharge(0),&(*magneticField));
          fts.setCartesianError(CartesianTrajectoryError(mu1TS.theState().cartesianError().matrix()+mu2TS.theState().cartesianError().matrix()));

          TransientTrack tt = theTTBuilder->build(fts);
          pair<bool,Measurement1D> delta3d = IPTools::absoluteImpactParameter3D(tt,thePrimaryV);
          if( delta3d.first ) {
            myCand.addUserFloat("delta3d",delta3d.second.value());
            myCand.addUserFloat("delta3dErr",delta3d.second.error());
          }

          if (addMuonlessPrimaryVertex_)
            myCand.addUserData("muonlessPV",Vertex(thePrimaryV));
          else
            myCand.addUserData("PVwithmuons",thePrimaryV);

          // lifetime using PV
          pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
          TVector3 vdiff = vtx - pvtx;
          double cosAlphaXY = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
          Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
          double ctauPV = distXY.value()*cosAlphaXY * myCand.mass()/pperp.Perp();
          GlobalError v1e = (Vertex(myVertex)).error();
          GlobalError v2e = thePrimaryV.error();
          AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
          double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*myCand.mass()/(pperp.Perp2());

          pvtx3D.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),thePrimaryV.position().z());
          TVector3 vdiff3D = vtx3D - pvtx3D;
          double cosAlpha3D = vdiff3D.Dot(pCand)/(vdiff3D.Mag()*pCand.Mag());
          Measurement1D dist3D = vdist3D.distance(Vertex(myVertex), thePrimaryV);
          double l3d = dist3D.value();
          double l3dsig = l3d/dist3D.error();

          myCand.addUserFloat("ppdlPV",ctauPV);
          myCand.addUserFloat("ppdlErrPV",ctauErrPV);
          myCand.addUserFloat("cosAlphaXY",cosAlphaXY);
          myCand.addUserFloat("cosAlpha3D",cosAlpha3D);
          myCand.addUserFloat("l3d",l3d);
          myCand.addUserFloat("l3dsig",l3dsig);

          // lifetime using BS
          pvtx.SetXYZ(theBeamSpotV.position().x(),theBeamSpotV.position().y(),0);
          vdiff = vtx - pvtx;
          cosAlphaXY = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
          distXY = vdistXY.distance(Vertex(myVertex), theBeamSpotV);
          double ctauBS = distXY.value()*cosAlphaXY*myCand.mass()/pperp.Perp();
          GlobalError v1eB = (Vertex(myVertex)).error();
          GlobalError v2eB = theBeamSpotV.error();
          AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
          double ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*myCand.mass()/(pperp.Perp2());

          myCand.addUserFloat("ppdlBS",ctauBS);
          myCand.addUserFloat("ppdlErrBS",ctauErrBS);

          if (addCommonVertex_) {
            myCand.addUserData("commonVertex",Vertex(myVertex));
          }
        } else {
          myCand.addUserFloat("vNChi2",-1);
          myCand.addUserFloat("vProb", -1);
          myCand.addUserFloat("PVndof", -1.);
          myCand.addUserFloat("pvw8", -1.);
          myCand.addUserFloat("ppdlPV",-100);
          myCand.addUserFloat("ppdlErrPV",-100);
          myCand.addUserFloat("cosAlphaXY",-100);
          myCand.addUserFloat("cosAlpha3D",-100);
          myCand.addUserFloat("l3d",-100);
          myCand.addUserFloat("l3dsig",-100);
          myCand.addUserFloat("ppdlBS",-100);
          myCand.addUserFloat("ppdlErrBS",-100);
          myCand.addUserFloat("Isolation", -1);
          myCand.addUserFloat("minDca", -1);
          myCand.addUserFloat("minDcaDr0p7", -1);
          myCand.addUserFloat("minDcaDr0p7Min0p5", -1);
          myCand.addUserFloat("minDcaDr0p7Min0p9", -1);
          myCand.addUserInt("Ntrk", 9999);
          myCand.addUserInt("countTksOfPV", 9999);
          myCand.addUserFloat("vertexWeight", -1);
          myCand.addUserFloat("sumPTPV", -1);
          myCand.addUserFloat("DCAXY", -1);
          myCand.addUserFloat("DCA", -1);
          myCand.addUserFloat("delta3d",-1);
          myCand.addUserFloat("delta3dErr",-1);
          if (addCommonVertex_) {
            myCand.addUserData("commonVertex",Vertex());
          }
          if (addMuonlessPrimaryVertex_) {
            myCand.addUserData("muonlessPV",Vertex());
          } else {
            myCand.addUserData("PVwithmuons",Vertex());
          }
        }
      }

      // ---- MC Truth, if enabled ----
      if (addMCTruth_) {
        reco::GenParticleRef genMu1 = it->genParticleRef();
        reco::GenParticleRef genMu2 = it2->genParticleRef();
        if (genMu1.isNonnull() && genMu2.isNonnull()) {
          if (genMu1->numberOfMothers()>0 && genMu2->numberOfMothers()>0){
            reco::GenParticleRef mom1 = genMu1->motherRef();
            reco::GenParticleRef mom2 = genMu2->motherRef();
            if (mom1.isNonnull() && (mom1 == mom2)) {
              myCand.setGenParticleRef(mom1); // set
              myCand.embedGenParticle();      // and embed
              std::pair<int, float> MCinfo = findJpsiMCInfo(mom1);
              myCand.addUserInt("momPDGId",MCinfo.first);
              myCand.addUserFloat("ppdlTrue",MCinfo.second);
            } else {
              myCand.addUserInt("momPDGId",0);
              myCand.addUserFloat("ppdlTrue",-99.);
            }
          } else {
            Handle<GenParticleCollection>theGenParticles;
            iEvent.getByLabel("genParticles", theGenParticles);
            if (theGenParticles.isValid()){
              for(size_t iGenParticle=0; iGenParticle<theGenParticles->size();++iGenParticle) {
                const Candidate & genCand = (*theGenParticles)[iGenParticle];
                if (genCand.pdgId()==443 || genCand.pdgId()==100443 ||
                    genCand.pdgId()==553 || genCand.pdgId()==100553 || genCand.pdgId()==200553) {
                  reco::GenParticleRef mom1(theGenParticles,iGenParticle);
                  myCand.setGenParticleRef(mom1);
                  myCand.embedGenParticle();
                  std::pair<int, float> MCinfo = findJpsiMCInfo(mom1);
                  myCand.addUserInt("momPDGId",MCinfo.first);
                  myCand.addUserFloat("ppdlTrue",MCinfo.second);
                }
              }
            } else {
              myCand.addUserInt("momPDGId",0);
              myCand.addUserFloat("ppdlTrue",-99.);
            }
          }
        } else {
          myCand.addUserInt("momPDGId",0);
          myCand.addUserFloat("ppdlTrue",-99.);
        }
      }


      // Save the candidate if it passes preselection cuts
      if( prefilter_(myCand) ) {
        oniaOutput->push_back(myCand);
        if (addThirdTrack_ && TrackFound) oniaOutputWithTrack->push_back(my3Cand);
      }
    }
  }

  // std::sort(oniaOutput->begin(),oniaOutput->end(),pTComparator_);
  std::sort(oniaOutput->begin(),oniaOutput->end(),vPComparator_);
  if (addThirdTrack_) std::sort(oniaOutputWithTrack->begin(),oniaOutputWithTrack->end(),vPComparator_);

  iEvent.put(oniaOutput,"DiMu");
  if (addThirdTrack_) iEvent.put(oniaOutputWithTrack,"DiMuTk");
}

bool Onia2MuMuPAT::isAbHadron(int pdgID)
{
  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;
}

bool Onia2MuMuPAT::isAMixedbHadron(int pdgID, int momPdgID)
{
  if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
    return true;
  return false;
}

std::pair<int, float>  
Onia2MuMuPAT::findJpsiMCInfo(reco::GenParticleRef genJpsi)
{
  int momJpsiID = 0;
  float trueLife = -99.;

  if (genJpsi->numberOfMothers()>0) {

    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    bool aBhadron = false;
    reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
    if (Jpsimom.isNull()) {
      std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
      return result;
    } else {
      reco::GenParticleRef Jpsigrandmom = Jpsimom->motherRef();
      if (isAbHadron(Jpsimom->pdgId())) {
        if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {
          momJpsiID = Jpsigrandmom->pdgId();
          trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
        } else {
          momJpsiID = Jpsimom->pdgId();
          trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
        }
        aBhadron = true;
      } else {
        if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId())) {
          reco::GenParticleRef JpsiGrandgrandmom = Jpsigrandmom->motherRef();
          if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
            momJpsiID = JpsiGrandgrandmom->pdgId();
            trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
          } else {
            momJpsiID = Jpsigrandmom->pdgId();
            trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
          }
          aBhadron = true;
        }
      }
      if (!aBhadron) {
        momJpsiID = Jpsimom->pdgId();
        trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
      }
    }

    TVector3 vdiff = trueVtx - trueVtxMom;
    //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
    trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
  }
  std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
  return result;

}

bool Onia2MuMuPAT::searchForTheThirdTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup, pat::CompositeCandidate& Cand, const pat::Muon* mu1, const pat::Muon* mu2,reco::Vertex& PV, reco::Track& ktk ,bool& closestPV){

  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  vector<double> muMasses;
  muMasses.push_back( 0.1134289256 );
  muMasses.push_back( 0.1134289256 );

  bool found=false;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  Handle<TrackCollection> tkColl;
  iEvent.getByLabel("generalTracks", tkColl);

  double tmpProb=-1, errorMass=-1;
  TransientVertex SVtx;
  Track track;
  int charge=10;

  for (TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk) {

    if (!tk->quality(TrackBase::highPurity)) continue;
    if (tk->pt() < minTrackPt_) continue;
    if (tk->pt()==mu1->innerTrack()->pt() || tk->pt()==mu2->innerTrack()->pt()) continue;

    TLorentzVector dimuP4=TLorentzVector(mu1->px()+mu2->px(),mu1->py()+mu2->py(),mu1->pz()+mu2->pz(),mu1->energy()+mu2->energy());
    TLorentzVector trackP4=TLorentzVector(tk->px(),tk->py(),tk->pz(), sqrt((tk->p()*tk->p())+(trackMass_*trackMass_)));
    TLorentzVector totP4=dimuP4+trackP4;

    if (totP4.M() > diMuPlusTrackMassMax_ || totP4.M() < diMuPlusTrackMassMin_) continue;

    vector<TransientTrack> t_tks;
    track=Track(*tk);

    t_tks.push_back(theTTBuilder->build(mu1->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
    t_tks.push_back(theTTBuilder->build(mu2->track())); // otherwise the vertex will have transient refs inside.
    t_tks.push_back(theTTBuilder->build(track));

    TransientVertex tmpVtx = vtxFitter.vertex(t_tks);

    if (!tmpVtx.isValid()) continue;

    CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
    Measurement1D MassWErr = massCalculator.invariantMass( VtxForInvMass, muMasses );

    double vChi2 = tmpVtx.totalChiSquared();
    double vNDF  = tmpVtx.degreesOfFreedom();

    double vProb(TMath::Prob(vChi2,(int)vNDF));

    if (vProb < tmpProb) continue;

    found=true;
    tmpProb=vProb;
    track=Track(*tk);
    SVtx=tmpVtx;
    errorMass=MassWErr.error();

  }

  if (found){
    
    ktk=Track(track);
    LorentzVector vtk=LorentzVector(track.px(),track.py(),track.pz(),sqrt((track.p()*track.p())+(trackMass_*trackMass_)));
    LorentzVector jpsi = mu1->p4() + mu2->p4();
    Cand.setP4(jpsi+vtk);

    // LeafCandidate tkcand=LeafCandidate(track.charge(),vtk);
    //Cand.addDaughter(tkcand,"Kaon");

    charge=mu1->charge()+mu2->charge()+track.charge();
    Cand.setCharge(charge);

    Cand.addUserFloat("MassErr",errorMass);

    double vChi2 = SVtx.totalChiSquared();
    double vNDF  = SVtx.degreesOfFreedom();
    double vProb(TMath::Prob(vChi2,(int)vNDF));

    Cand.addUserFloat("vNChi2",vChi2/vNDF);
    Cand.addUserFloat("vProb",vProb);

    //now the primary vertex
    Vertex thePrimaryV;
    Vertex theBeamSpotV;

    ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

    Handle<BeamSpot> theBeamSpot;
    iEvent.getByLabel(thebeamspot_,theBeamSpot);
    BeamSpot bs = *theBeamSpot;
    theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

    Handle<VertexCollection> priVtxs;
    iEvent.getByLabel(thePVs_, priVtxs);
    if ( priVtxs->begin() != priVtxs->end() ) {
      thePrimaryV = Vertex(*(priVtxs->begin()));
    }
    else {
      thePrimaryV = Vertex(bs.position(), bs.covariance3D());
    }

    if (closestPV){

      float minDz = 999999.;
      TwoTrackMinimumDistance ttmd;
      bool status = ttmd.calculate( GlobalTrajectoryParameters(
                                      GlobalPoint(SVtx.position().x(), SVtx.position().y(), SVtx.position().z()),
                                      GlobalVector(Cand.px(),Cand.py(),Cand.pz()),TrackCharge(charge),&(*magneticField)),
                                    GlobalTrajectoryParameters(
                                      GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
                                      GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(charge),&(*magneticField)));
      float extrapZ=-9E20;
      if (status) extrapZ=ttmd.points().first.z();
      
      for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv){
        float deltaZ = fabs(extrapZ - itv->position().z()) ;
        if ( deltaZ < minDz ) {
          minDz = deltaZ;
          thePrimaryV = Vertex(*itv);
        }
      }
    }

    PV=Vertex(thePrimaryV);

    // count the number of high Purity tracks with pT > 900 MeV attached to the chosen vertex
    double vertexWeight = 0., sumPTPV = 0.;
    int countTksOfPV = 0, Ntrk=0;
    double minDca = 9999.;

    for(reco::Vertex::trackRef_iterator itVtx = thePrimaryV.tracks_begin(); itVtx != thePrimaryV.tracks_end(); itVtx++){
      if(itVtx->isNonnull()){
        const reco::Track& trk = **itVtx;
        if (trk.pt()==mu1->innerTrack()->pt() || trk.pt()==mu2->innerTrack()->pt() || trk.pt()==track.pt()) continue;
        //if(!track.quality(reco::TrackBase::highPurity)) continue;
        TransientTrack tt = theTTBuilder->build(trk);
        pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,SVtx);
        vertexWeight += thePrimaryV.trackWeight(*itVtx);
        if(track.pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;
        if(deltaR(track.eta(),track.phi(),Cand.eta(),Cand.phi()) < 0.7) {
          // cout<<"myCand.eta"<<myCand.eta()<<"myCand.phi"<<myCand.phi()<<endl;
          if (tkPVdist.first && tkPVdist.second.value()<minDca){
            minDca = tkPVdist.second.value();
          }
          if(track.pt() < 0.9) continue; //reject all rejects from counting if less than 900 MeV
          sumPTPV += track.pt();
          countTksOfPV++;
        }
      }
    }

    // now tracks outside the PV fit
    double sumPt_noVtx = 0.;
    int countTksOfnoVtx = 0;

    for (TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk) {
      bool hasVertex = false;
      for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv) {
        for( reco::Vertex::trackRef_iterator itVtx = itv->tracks_begin(); itVtx != itv->tracks_end(); ++itVtx) {
          if( itVtx->isNonnull() ) {
            if(tk->pt() == (*itVtx)->pt() && tk->eta()==(*itVtx)->eta()) {
              hasVertex=true;
              if (!hasVertex && (tk->pt()!=mu1->innerTrack()->pt() && tk->pt()!=mu2->innerTrack()->pt() && tk->pt()==track.pt())){
                TransientTrack tt = theTTBuilder->build(*tk);
                pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,SVtx);
                if(tk->pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;
                if (deltaR(tk->eta(),tk->phi(),Cand.eta(), Cand.phi())< 0.7) {
                  if (tkPVdist.first && tkPVdist.second.value()<minDca){
                    minDca = tkPVdist.second.value();
                  }
                  if(tk->pt()<0.9) continue;
                  if(tkPVdist.second.value() > 0.05) continue;
                  sumPt_noVtx += tk->pt();
                  countTksOfnoVtx++;
                }
              }
            }
          }
        }
      }
    }
    
    double Iso = Cand.pt()/(Cand.pt()+sumPt_noVtx+sumPTPV);

    Cand.addUserFloat("Isolation", (float) Iso);
    Cand.addUserFloat("minDca", (float) minDca);
    Cand.addUserInt("Ntrk", Ntrk);
    Cand.addUserInt("countTksOfPV", countTksOfPV);
    Cand.addUserFloat("vertexWeight", (float) vertexWeight);
    Cand.addUserFloat("sumPTPV", (float) sumPTPV);

    TrajectoryStateClosestToPoint mu1TS = (theTTBuilder->build(mu1->track())).impactPointTSCP();
    TrajectoryStateClosestToPoint mu2TS = (theTTBuilder->build(mu2->track())).impactPointTSCP();
    TrajectoryStateClosestToPoint tkTS = (theTTBuilder->build(track)).impactPointTSCP();

    FreeTrajectoryState fts(
          GlobalPoint(SVtx.position().x(), SVtx.position().y(), SVtx.position().z()),
          GlobalVector(Cand.px(),Cand.py(),Cand.pz()),TrackCharge(charge),&(*magneticField));
    fts.setCartesianError(CartesianTrajectoryError(mu1TS.theState().cartesianError().matrix()+mu2TS.theState().cartesianError().matrix()+tkTS.theState().cartesianError().matrix()));
    
    TransientTrack tt = theTTBuilder->build(fts);
    pair<bool,Measurement1D> delta3d = IPTools::absoluteImpactParameter3D(tt,thePrimaryV);

    if( delta3d.first ) {
      Cand.addUserFloat("delta3d",delta3d.second.value());
      Cand.addUserFloat("delta3dErr",delta3d.second.error());
    }

    //lifetime etc...

    TVector3 vtx;
    TVector3 pvtx;
    VertexDistanceXY vdistXY;
    vtx.SetXYZ(SVtx.position().x(),SVtx.position().y(),0);
    TVector3 pperp(Cand.px(), Cand.py(), 0);
    AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

    //using PV
    pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
    TVector3 vdiff = vtx - pvtx;
    double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
    Measurement1D distXY = vdistXY.distance(Vertex(SVtx), thePrimaryV);
    double ctauPV = distXY.value()*cosAlpha * Cand.mass()/pperp.Perp();
    GlobalError v1e = (Vertex(SVtx)).error();
    GlobalError v2e = thePrimaryV.error();
    AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
    double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*Cand.mass()/(pperp.Perp2());
    
    Cand.addUserFloat("ppdlPV",ctauPV);
    Cand.addUserFloat("ppdlErrPV",ctauErrPV);
    Cand.addUserFloat("cosAlpha",cosAlpha);
    
    //using BS
    pvtx.SetXYZ(theBeamSpotV.position().x(),theBeamSpotV.position().y(),0);
    vdiff = vtx - pvtx;
    cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
    distXY = vdistXY.distance(Vertex(SVtx), theBeamSpotV);
    double ctauBS = distXY.value()*cosAlpha*Cand.mass()/pperp.Perp();
    GlobalError v1eB = (Vertex(SVtx)).error();
    GlobalError v2eB = theBeamSpotV.error();
    AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
    double ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*Cand.mass()/(pperp.Perp2());

    Cand.addUserFloat("ppdlBS",ctauBS);
    Cand.addUserFloat("ppdlErrBS",ctauErrBS);

  }else{
    Cand.addUserFloat("vNChi2",-1);
    Cand.addUserFloat("vProb", -1);
    Cand.addUserFloat("ppdlPV",-100);
    Cand.addUserFloat("ppdlErrPV",-100);
    Cand.addUserFloat("cosAlpha",-100);
    Cand.addUserFloat("ppdlBS",-100);
    Cand.addUserFloat("ppdlErrBS",-100);
    Cand.addUserFloat("Isolation", -1);
    Cand.addUserFloat("minDca", -1);
    Cand.addUserInt("Ntrk", 9999);
    Cand.addUserInt("countTksOfPV", 9999);
    Cand.addUserFloat("vertexWeight", -1);
    Cand.addUserFloat("sumPTPV", -1);
    Cand.addUserFloat("delta3d",-1);
    Cand.addUserFloat("delta3dErr",-1);
    Cand.addUserFloat("MassErr",-1);
    Cand.addUserFloat("DCAXY", -1);
    Cand.addUserFloat("DCA", -1);
  }
  
  return found;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Onia2MuMuPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Onia2MuMuPAT::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(Onia2MuMuPAT);
