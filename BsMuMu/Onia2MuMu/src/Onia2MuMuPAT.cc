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
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

Onia2MuMuPAT::Onia2MuMuPAT(const edm::ParameterSet& iConfig)
  :
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
  muon_mass = 0.1056583;
  muon_sigma = 0.0000001;
  muMasses_.push_back( 0.1134289256 );
  muMasses_.push_back( 0.1134289256 );
}


Onia2MuMuPAT::~Onia2MuMuPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

void Onia2MuMuPAT::fillCandMuons(pat::CompositeCandidate & myCand, const pat::Muon & mu1, const pat::Muon & mu2)
{
  myCand.addDaughter(mu1, "muon1");
  int highPurity1 = 0;
  if(mu1.track()->quality(reco::TrackBase::highPurity)) highPurity1 = 1;
  myCand.addUserInt("highPurity1", highPurity1);
  myCand.addUserFloat("segComp1", muon::segmentCompatibility(mu1));

  myCand.addDaughter(mu2, "muon2");
  int highPurity2 = 0;
  if(mu2.track()->quality(reco::TrackBase::highPurity)) highPurity2 = 1;
  myCand.addUserInt("highPurity2", highPurity2);
  myCand.addUserFloat("segComp2", muon::segmentCompatibility(mu2));
}

math::XYZTLorentzVector Onia2MuMuPAT::fromPtEtaPhiToPxPyPz( const double & pt, const double & eta, const double & phi, const double & mass) const
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
  std::cout << iEvent.id() << std::endl;

  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutput(new pat::CompositeCandidateCollection);
  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutputWithTrack(new pat::CompositeCandidateCollection);

  std::pair<reco::Vertex, reco::Vertex> thePrimaryV;
  Vertex theBeamSpotV;

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel(thePVs_, priVtxs);
//  if ( priVtxs->begin() != priVtxs->end() ) {
//    thePrimaryV = Vertex(*(priVtxs->begin()));
//  }
//  else {
//    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
//  }

  Handle< View<pat::Muon> > muons;
  iEvent.getByLabel(muons_,muons);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  //Creating a KinematicParticleFactory
  KinematicParticleFactoryFromTransientTrack pFactory;
  KinematicParticleVertexFitter fitter;

  // dimuon candidates only from muons
  for(View<pat::Muon>::const_iterator it = muons->begin(), itend = muons->end(); it != itend; ++it){
    // both must pass low quality
    if(!lowerPuritySelection_(*it)) continue;
    if( !(it->innerTrack()->quality(TrackBase::highPurity)) ) continue;
    for(View<pat::Muon>::const_iterator it2 = it+1; it2 != itend;++it2){
      // both must pass low quality
      if(!lowerPuritySelection_(*it2)) continue;
      if( !(it2->innerTrack()->quality(TrackBase::highPurity)) ) continue;
      // one must pass tight quality
      if (!(higherPuritySelection_(*it) || higherPuritySelection_(*it2))) continue;

      // All the tracks
      edm::Handle<reco::TrackCollection> tkColl;
      iEvent.getByLabel("generalTracks", tkColl);

      // Rereco PVs removing the two muons and adding the beamspot. We only refit existing PVs.
      const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(it->originalObject());
      const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(it2->originalObject());

      // Do it only for the 3 closest
      pat::CompositeCandidate myCand;
      pat::CompositeCandidate my3Cand;
      bool TrackFound=false;
      Vertex PVDiMuTk,SVDiMuTk;
      Track track;
      LeafCandidate tkcand;

      double originalMu1P = it->track()->p();
      double originalMu2P = it2->track()->p();
      double originalMu1Eta = it->track()->eta();
      double originalMu2Eta = it2->track()->eta();
      double originalMu1Phi = it->track()->phi();
      double originalMu2Phi = it2->track()->phi();

      // Compute invariant mass to reject candidates before the vertex fit
      math::XYZTLorentzVector m1p4(fromPtEtaPhiToPxPyPz(it->track()->pt(),
                                                        it->track()->eta(),
                                                        it->track()->phi(),
                                                        muMasses_[0]));
      math::XYZTLorentzVector m2p4(fromPtEtaPhiToPxPyPz(it2->track()->pt(),
                                                        it2->track()->eta(),
                                                        it2->track()->phi(),
                                                        muMasses_[1]));
      LorentzVector dimuonBeforeRefit = m1p4 + m2p4;

      // ---- apply the dimuon cut ----
      if(!dimuonSelection_(dimuonBeforeRefit)) continue;

      // ---- fit vertex using Tracker tracks (if they have tracks) ----
      if (it->track().isNonnull() && it2->track().isNonnull()) {

        if(deltaR(m1p4.eta(),m1p4.phi(),m2p4.eta(), m2p4.phi()) < 1.e-5) continue;

        // Dimu+track
        // -------------------------------------------------------------------------------------------------------------------------------
        double M2=myCand.mass();
        if (addThirdTrack_ && (it->charge()!=it2->charge()) && M2 > diMuMassMin_ && M2 < diMuMassMax_)
          TrackFound=searchForTheThirdTrack(iEvent, iSetup, my3Cand, &(*it), &(*it2),PVDiMuTk,track,resolveAmbiguity_);

        if (TrackFound) {
          if (addMuonlessPrimaryVertex_) my3Cand.addUserData("muonlessPV",PVDiMuTk);
          else my3Cand.addUserData("PVwithmuons",PVDiMuTk);

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
        }
        else {
          if (addCommonVertex_) my3Cand.addUserData("commonVertex",Vertex());
          if (addMuonlessPrimaryVertex_) my3Cand.addUserData("muonlessPV",Vertex());
          else my3Cand.addUserData("PVwithmuons",Vertex());
        }
        // -------------------------------------------------------------------------------------------------------------------------------


        // -------------------------------------------------------------------------------------------------------------------------------
        // Secondary vertex fit
        // -------------------------------------------------------------------------------------------------------------------------------
        vector<TransientTrack> t_tks;
        t_tks.push_back(theTTBuilder->build(*it->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(*it2->track())); // otherwise the vertex will have transient refs inside.

        // -------------------------------------------------------------------------------------------------------------------------------
        // Kinematic vertex fit
        // -------------------------------------------------------------------------------------------------------------------------------
        std::cout << "starting" << std::endl;
        vector<RefCountedKinematicParticle> muonParticles;
        muonParticles.push_back(pFactory.particle (t_tks[0],muon_mass,0.,0.,muon_sigma));
        muonParticles.push_back(pFactory.particle (t_tks[1],muon_mass,0.,0.,muon_sigma));
        RefCountedKinematicTree vertexFitTree = fitter.fit(muonParticles);
        // --------------------
        if( vertexFitTree->isValid() ) {
          RefCountedKinematicVertex myVertex = vertexFitTree->currentDecayVertex();
          if( myVertex->vertexIsValid() ) {
            // Prepare a recoVertex for later use
            reco::Vertex myRecoVertex(*myVertex);
            //accessing the tree components, move pointer to top
            vertexFitTree->movePointerToTheTop();

            //We are now at the top of the decay tree getting the B_s reconstructed KinematicPartlcle
            RefCountedKinematicParticle b_s = vertexFitTree->currentParticle();

            reco::TransientTrack tt = theTTBuilder->build(b_s->currentState().freeTrajectoryState());

            // Mass from the kinematic vertex fit
            myCand.addUserFloat("mass",b_s->currentState().mass());
            myCand.addUserFloat("kinMass",b_s->currentState().mass());
            // Mass is parameter [6], the error matrix has 7 entries and parameter [6] should be the error on the mass.
            // According to https://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_3_8/doc/html/d4/de9/classKinematicParametersError.html
            myCand.addUserFloat("kinMassError",sqrt(b_s->currentState().kinematicParametersError().matrix()(6,6)));

            float vChi2 = myVertex->chiSquared();
            float vNDF  = myVertex->degreesOfFreedom();
            float vProb(TMath::Prob(vChi2,(int)vNDF));

            // Get first muon
            pat::Muon mu1(*it);
            if( vertexFitTree->movePointerToTheFirstChild() ) setP4AndCharge(vertexFitTree, mu1, theTTBuilder);
            // Get second muon
            pat::Muon mu2(*it2);
            if( vertexFitTree->movePointerToTheNextChild() ) setP4AndCharge(vertexFitTree, mu2, theTTBuilder);
            // Make sure mu1 is the highest pt muon
            if( mu1.pt() > mu2.pt() ) {
              fillCandMuons(myCand, mu1, mu2);
              fillCandMuons(my3Cand, mu1, mu2);
            }
            else {
              fillCandMuons(myCand, mu2, mu1);
              fillCandMuons(my3Cand, mu2, mu1);
            }
            std::cout << "mu1 pt = " << myCand.daughter("muon1")->pt() << std::endl;
            std::cout << "mu2 pt = " << myCand.daughter("muon2")->pt() << std::endl;

            myCand.addUserFloat("vNChi2",vChi2/vNDF);
            myCand.addUserFloat("vProb",vProb);

            LorentzVector dimuon = mu1.p4() + mu2.p4();
            myCand.setP4(dimuon);
            myCand.setCharge(mu1.charge()+mu2.charge());

            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;
            TVector3 vtx3D;
            TVector3 pvtx3D;
            VertexDistance3D vdist3D;

            // XY
            vtx.SetXYZ(myVertex->position().x(),myVertex->position().y(),0);
            TVector3 pperp(dimuon.px(), dimuon.py(), 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            // 3D
            vtx3D.SetXYZ(myVertex->position().x(),myVertex->position().y(),myVertex->position().z());
            TVector3 pCand(dimuon.px(), dimuon.py(), dimuon.pz());
            // -------------------------------------------------------------------------------------------------------------------------------

            // Refit the PVs
            std::vector<std::pair<reco::Vertex, reco::Vertex> > rePVs;
            for( auto PVit = priVtxs->begin(); PVit != priVtxs->end(); ++PVit ) {
              rePVs.push_back(buildMuonlessPV(rmu1, rmu2, muonLess_, *PVit, theTTBuilder, bs));
            }

            // Select the PV
            if (resolveAmbiguity_) {

              float minDz = 999999.;
              TwoTrackMinimumDistance ttmd;
              bool status = ttmd.calculate( GlobalTrajectoryParameters(GlobalPoint(myVertex->position().x(), myVertex->position().y(), myVertex->position().z()),
                                                                       GlobalVector(myCand.px(),myCand.py(),myCand.pz()),TrackCharge(0),&(*magneticField)),
                                            GlobalTrajectoryParameters(GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
                                                                       GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)) );
              float extrapZ=-9E20;
              if (status) extrapZ=ttmd.points().first.z();

              // Select the closest PV in z
              for( auto tPV = rePVs.begin(); tPV != rePVs.end(); ++tPV ) {
                float deltaZ = fabs(extrapZ - tPV->first.position().z()) ;
                if ( deltaZ < minDz ) {
                  minDz = deltaZ;
                  // thePrimaryV = Vertex(*tPV);
                  thePrimaryV = *tPV;
                }
              }
            }
            // Using the refitted PV
            pair<bool,Measurement1D> dist = IPTools::absoluteImpactParameter3D(tt, thePrimaryV.first);
            if( dist.first ) {
              myCand.addUserFloat("delta3d",dist.second.value());
              myCand.addUserFloat("delta3dErr",dist.second.error());
            }
            pair<bool,Measurement1D> pvlip = IPTools::signedDecayLength3D(tt,GlobalVector(0,0,1),thePrimaryV.first);
            if( pvlip.first ) {
              myCand.addUserFloat("pvlip",pvlip.second.value());
              myCand.addUserFloat("pvlipErr",pvlip.second.error());
            }

            // count the number of tracks with pT > 900 MeV attached to the chosen vertex
            double sumPTPV = 0., sumPTNoVtx = 0.;
            double sumPTPV20 = 0., sumPTNoVtx20 = 0.;
            double sumPTkMu1 = 0., sumPTkMu2 = 0.;
            double sumPTkMu1_20 = 0., sumPTkMu2_20 = 0.;
            int countTksOfPV = 0, countTksOfNoVtx = 0;
            int Ntrk = 0, Ntrkhp = 0, Ntrk20 = 0, Ntrkhp20 = 0;
            int Ntrk1sigma = 0, Ntrk1sigmahp = 0, Ntrk1sigma20 = 0, Ntrk1sigmahp20 = 0;
            int Ntrk2sigma = 0, Ntrk2sigmahp = 0, Ntrk2sigma20 = 0, Ntrk2sigmahp20 = 0;
            int Ntrk3sigma = 0, Ntrk3sigmahp = 0, Ntrk3sigma20 = 0, Ntrk3sigmahp20 = 0;
            double minDca = 9999.;
            double sumNdofPV = 0., sumNdofNoVtx = 0.;
            double PVndof = thePrimaryV.first.ndof();
            myCand.addUserFloat("PVndof",(float) PVndof);
            myCand.addUserFloat("pvw8Refit",(float) (thePrimaryV.first.ndof()+2)/(2*thePrimaryV.first.tracksSize()));
            myCand.addUserFloat("pvw8",(float) (thePrimaryV.second.ndof()+2)/(2*thePrimaryV.second.tracksSize()));

            TrajectoryStateClosestToPoint mu1TS = t_tks[0].impactPointTSCP();
            TrajectoryStateClosestToPoint mu2TS = t_tks[1].impactPointTSCP();

            // Fill a collection of all tracks and:
            // - if they belong to the PV associated to the candidate set vertexId to 1
            // - if they belong to another PV set the vertexId to 2
            // - if they belong to no PV set the vertexId to 0
            double maxDeltaR = 1.e-5;
            std::vector<SimpleTrack> simpleTracks;
            for( auto track = tkColl->begin(); track != tkColl->end(); ++track ) {
              // Exclude the two muons
              // if( rmu1 != 0 && deltaR(track->eta(),track->phi(),rmu1->eta(),rmu1->phi()) < maxDeltaR ) continue;
              // if( rmu2 != 0 && deltaR(track->eta(),track->phi(),rmu2->eta(),rmu2->phi()) < maxDeltaR ) continue;
              if( rmu1 != 0 && track->extra().key() == rmu1->track().key() ) continue;
              if( rmu2 != 0 && track->extra().key() == rmu2->track().key() ) continue;

              // Set defaults
              double doca = 999.;
              double docaSig = 999.;
              int vertexId = 999;
              bool highPurity = track->quality(reco::TrackBase::highPurity);
              // Compute doca and doca significance
              TransientTrack tt = theTTBuilder->build(*track);
              // Looking in the code it seems absoluteImpactParameter3D uses only the position of the vertex
              pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt, myRecoVertex);
              if( tkPVdist.first ) {
                doca = tkPVdist.second.value();
                if(tkPVdist.second.error() != 0) docaSig = doca/tkPVdist.second.error();
              }
              vertexId = findVertexId( thePrimaryV, rePVs, *track, maxDeltaR );

              // Muon isolation
              // --------------
              double isoMu1P = 0.;
              double deltaRMu1 = 0.;
              bool usedPVMu1 = false;
              double isoMu2P = 0.;
              double deltaRMu2 = 0.;
              bool usedPVMu2 = false;
              // std::cout << "doca = " << doca << ", p = " << track->p() << std::endl;
              // if( doca < 0.1 && track->p() > 0.5 ) {
              //   if( vertexId == 0 || vertexId == 1 ) {
              // if( (vertexId == 0 && doca < 0.1) || vertexId == 1 ) {
              if( vertexId != 2 ) {
                // if( track->p() > 0.5 ) {
                // if( deltaR(track->eta(),track->phi(),rmu1->eta(),rmu1->phi()) < 0.5 ) {
                    //                    std::cout << "vertexId = " << vertexId << std::endl;
                    //                    std::cout << "doca = " << doca << ", p = " << track->p() << std::endl;
                    //                    std::cout << "mu1_20 iso deltaR = " << deltaR(track->eta(),track->phi(),rmu1->eta(),rmu1->phi()) << std::endl;

                // deltaRMu1 = deltaR(track->eta(),track->phi(),rmu1->eta(),rmu1->phi());
                deltaRMu1 = deltaR(track->eta(),track->phi(),originalMu1Eta,originalMu1Phi);
                // deltaRMu1 = deltaR(track->eta(),track->phi(),myCand.daughter("muon1")->eta(),myCand.daughter("muon1")->phi());
                isoMu1P = track->p();
                if( deltaRMu1 < 0.5 && isoMu1P > 0.5 && vertexId == 1 ) usedPVMu1 = true;
                  //}
                  // if( deltaR(track->eta(),track->phi(),rmu2->eta(),rmu2->phi()) < 0.5 ) {
                    //                    std::cout << "vertexId = " << vertexId << std::endl;
                    //                    std::cout << "doca = " << doca << ", p = " << track->p() << std::endl;
                    //                    std::cout << "mu2_20 iso deltaR = " << deltaR(track->eta(),track->phi(),rmu2->eta(),rmu2->phi()) << std::endl;
                // deltaRMu2 = deltaR(track->eta(),track->phi(),rmu2->eta(),rmu2->phi());
                deltaRMu2 = deltaR(track->eta(),track->phi(),originalMu2Eta,originalMu2Phi);
                // deltaRMu2 = deltaR(track->eta(),track->phi(),myCand.daughter("muon2")->eta(),myCand.daughter("muon2")->phi());
                isoMu2P = track->p();
                if( deltaRMu2 < 0.5 && isoMu2P > 0.5 && vertexId == 1 ) usedPVMu2 = true;
                  // }
                // }
              }
              simpleTracks.push_back(SimpleTrack(track->pt(), track->eta(), track->phi(), track->ndof(),
                                                 doca, docaSig, vertexId, highPurity,
                                                 deltaRMu1, isoMu1P, usedPVMu1,
                                                 deltaRMu2, isoMu2P, usedPVMu2));
            }

            // Sort the simpleTracks. Smallest doca first.
            std::sort(simpleTracks.begin(), simpleTracks.end());

            // Count the number of tracks from the candidate PV or no PV and with pt > 0.5 GeV, doca < 300 microns and doca significance > 1.
            // Using only the first 20 tracks
            double maxDoca = 0.03;
            int tkCount = 0;
            for( auto tk = simpleTracks.begin(); tk != simpleTracks.end(); ++tk, ++tkCount ) {

              // From here on use only tracks from the candidate PV or no PV
              if( tk->vertexId == 2 ) continue;

              if( tk->vertexId == 1 ) ++countTksOfPV;
              else if( tk->vertexId == 0 ) ++countTksOfNoVtx;

//              if( tk->isoMu1P < 0.63 && tk->isoMu1P > 0.628 ) {
//                std::cout << "NOCUT: tk->p = " << tk->isoMu1P << ", deltaRMu1 = " << tk->deltaRMu1 << ", deltaRMu2 = " << tk->deltaRMu2 << ", doca = " << tk->doca << std::endl;
//              }

              double deltaRMuIsoCut = 0.5;

              // if( iEvent.id().event() == 238755924 || iEvent.id().event() == 238795635 ) {
                if( tk->deltaRMu1 < deltaRMuIsoCut || tk->deltaRMu2 < deltaRMuIsoCut ) {
                  if( (tk->vertexId == 0 && tk->doca < 0.1 && tkCount < 20) || (tk->vertexId == 1) ) {
                  // if( tk->vertexId == 1 ) {
                    if( tk->deltaRMu1 < deltaRMuIsoCut ) std::cout << "For mu1, ";
                    if( tk->deltaRMu2 < deltaRMuIsoCut ) std::cout << "For mu2, ";
                    if( tk->vertexId == 1 ) std::cout << "PV:   ";
                    if( tk->vertexId == 0 ) std::cout << "NOPV: ";
                    std::cout << "tk->p = " << tk->isoMu1P
                              // << ", tk->pt = " << tk->pt
                              << ", deltaRMu1 = " << tk->deltaRMu1
                              << ", deltaRMu2 = " << tk->deltaRMu2
                              << ", doca = " << tk->doca << std::endl;
                    // }
                  }
                }
              // }


              if( tk->vertexId == 1 ) {
                if( tk->deltaRMu1 < deltaRMuIsoCut && tk->isoMu1P > 0.5 ) sumPTkMu1 += tk->isoMu1P;
                // if( tk->isoMu1P > 0. ) std::cout << "mu1 iso track = " << tk->isoMu1P << std::endl;
                if( tk->deltaRMu2 < deltaRMuIsoCut && tk->isoMu2P > 0.5 ) sumPTkMu2 += tk->isoMu2P;
                // if( tk->isoMu2P > 0. ) std::cout << "mu2 iso track = " << tk->isoMu2P << std::endl;
              }
              // if( tkCount < 20 ) {
              //   sumPTkMu1_20 += tk->isoMu1P;
                // if( tk->isoMu1P > 0. ) std::cout << "mu1_20 iso track = " << tk->isoMu1P << std::endl;
              //   sumPTkMu2_20 += tk->isoMu2P;
                // if( tk->isoMu2P > 0. ) std::cout << "mu2_20 iso track = " << tk->isoMu2P << std::endl;
              // }



//              if( tkCount < 20 && tk->vertexId == 0 && tk->doca < 0.1 ) {
//                if( tk->deltaRMu1 < 0.5 && tk->isoMu1P > 0.5 ) {
//                  // if( ((tk->vertexId == 1 && tk->isoMu1P < 0.5) ||
//                  sumPTkMu1_20 += tk->isoMu1P;
//                }
//                if( tk->deltaRMu2 < 0.5 && tk->isoMu2P > 0.5 ) {
//                  sumPTkMu2_20 += tk->isoMu2P;
//                }
//              }



              if( tk->vertexId == 0 && tk->doca < 0.1 && tkCount < 20 ) {
                if( ( tk->deltaRMu1 < 0.5 && tk->isoMu1P > 0.5) )
                  sumPTkMu1_20 += tk->isoMu1P;
                if( ( tk->deltaRMu2 < 0.5 && tk->isoMu2P > 0.5) )
                  sumPTkMu2_20 += tk->isoMu2P;
              }



              // Take the minimum doca
              if( minDca > tk->doca ) minDca = tk->doca;

              // Number of close tracks
              if( tk->pt > 0.5 && tk->doca < maxDoca ) {
                ++Ntrk;
                if( tk->highPurity ) ++Ntrkhp;
                if( tkCount < 20 ) {
                  ++Ntrk20;
                  if( tk->highPurity ) ++Ntrkhp20;
                }

                // Number of close tracks with doca > 1 sigma
                if( tk->docaSignificance > 1. ) {
                  ++Ntrk1sigma;
                  if( tk->highPurity ) ++Ntrk1sigmahp;
                  if( tkCount < 20 ) {
                    ++Ntrk1sigma20;
                    if( tk->highPurity ) ++Ntrk1sigmahp20;
                  }

                  // Number of close tracks with doca > 2 sigma
                  if( tk->docaSignificance > 2. ) {
                    ++Ntrk2sigma;
                    if( tk->highPurity ) ++Ntrk2sigmahp;
                    if( tkCount < 20 ) {
                      ++Ntrk2sigma20;
                      if( tk->highPurity ) ++Ntrk2sigmahp20;
                    }

                    // Number of close tracks with doca > 3 sigma
                    if( tk->docaSignificance > 3. ) {
                      ++Ntrk3sigma;
                      if( tk->highPurity ) ++Ntrk3sigmahp;
                      if( tkCount < 20 ) {
                        ++Ntrk3sigma20;
                        if( tk->highPurity ) ++Ntrk3sigmahp20;
                      }
                    }
                  }
                }
              }

              // SumPT of tracks from candidate PV or no PV
              if( tk->pt < 0.9 ) continue;
              if( deltaR(tk->eta ,tk->phi, myCand.eta(), myCand.phi()) > 0.7 ) continue;
              if( tk->vertexId == 1 ) {
                sumNdofPV += tk->ndof;
                // std::cout << "track from PV pt = " << tk->pt << std::endl;
                sumPTPV += tk->pt;
                if( tkCount < 20 ) sumPTPV20 += tk->pt;
              }
              // In this case the doca is also required to be less than 500 microns
              else if( tk->vertexId == 0 && tk->doca < 0.05 ) {
                sumNdofNoVtx += tk->ndof;
                // std::cout << "track from no PV pt = " << tk->pt << std::endl;
                sumPTNoVtx += tk->pt;
                if( tkCount < 20 ) sumPTNoVtx20 += tk->pt;
              }
            }
            double Iso = myCand.pt()/(myCand.pt()+sumPTNoVtx+sumPTPV);
            double Iso20 = myCand.pt()/(myCand.pt()+sumPTNoVtx20+sumPTPV20);

            // Muon isolation
            double isoMu1 = mu1.p()/(mu1.p()+sumPTkMu1);
            double isoMu2 = mu2.p()/(mu2.p()+sumPTkMu2);
            double isoMu1_20 = originalMu1P/(originalMu1P+sumPTkMu1+sumPTkMu1_20);
            double isoMu2_20 = originalMu2P/(originalMu2P+sumPTkMu2+sumPTkMu2_20);
            // double isoMu1_20 = myCand.daughter("muon1")->p()/(myCand.daughter("muon1")->p()+sumPTkMu1+sumPTkMu1_20);
            // double isoMu2_20 = myCand.daughter("muon2")->p()/(myCand.daughter("muon2")->p()+sumPTkMu2+sumPTkMu2_20);
            // double isoMu1_20 = mu1.p()/(mu1.p()+sumPTkMu1+sumPTkMu1_20);
            // double isoMu2_20 = mu2.p()/(mu2.p()+sumPTkMu2+sumPTkMu2_20);

            myCand.addUserFloat("Isolation", (float) Iso);
            myCand.addUserFloat("Isolation20", (float) Iso20);
            myCand.addUserFloat("minDca", (float) minDca);
            myCand.addUserInt("Ntrk", Ntrk);
            myCand.addUserInt("Ntrkhp", Ntrkhp);
            myCand.addUserInt("Ntrk20", Ntrk20);
            myCand.addUserInt("Ntrkhp20", Ntrkhp20);
            myCand.addUserInt("Ntrk1sigma", Ntrk1sigma);
            myCand.addUserInt("Ntrk1sigmahp", Ntrk1sigmahp);
            myCand.addUserInt("Ntrk1sigma20", Ntrk1sigma20);
            myCand.addUserInt("Ntrk1sigmahp20", Ntrk1sigmahp20);
            myCand.addUserInt("Ntrk2sigma", Ntrk2sigma);
            myCand.addUserInt("Ntrk2sigmahp", Ntrk2sigmahp);
            myCand.addUserInt("Ntrk2sigma20", Ntrk2sigma20);
            myCand.addUserInt("Ntrk2sigmahp20", Ntrk2sigmahp20);
            myCand.addUserInt("Ntrk3sigma", Ntrk3sigma);
            myCand.addUserInt("Ntrk3sigmahp", Ntrk3sigmahp);
            myCand.addUserInt("Ntrk3sigma20", Ntrk3sigma20);
            myCand.addUserInt("Ntrk3sigmahp20", Ntrk3sigmahp20);
            myCand.addUserInt("countTksOfPV", countTksOfPV);
            myCand.addUserInt("countTksOfNoVtx", countTksOfNoVtx);
            myCand.addUserFloat("sumPTPV", (float) sumPTPV);
            myCand.addUserFloat("sumPTNoVtx", (float) sumPTNoVtx);
            myCand.addUserFloat("sumNdofPV", (float) sumNdofPV);
            myCand.addUserFloat("sumNdofNoVtx", (float) sumNdofNoVtx);
            myCand.addUserFloat("isoMu1", (float) isoMu1);
            myCand.addUserFloat("isoMu2", (float) isoMu2);
            myCand.addUserFloat("isoMu1_20", (float) isoMu1_20);
            myCand.addUserFloat("isoMu2_20", (float) isoMu2_20);

            // DCA
            // Compute both the dca in the transverse plane (dcaxy) and in 3D.
            float dcaxy = computeDcaXY(mu1TS, mu2TS);
            float dca = computeDca(mu1TS, mu2TS);
            myCand.addUserFloat("DCAXY", dcaxy);
            myCand.addUserFloat("DCA", dca);
            if (TrackFound){
              my3Cand.addUserFloat("DCAXY", dcaxy);
              my3Cand.addUserFloat("DCA", dca);
            }
            // end DCA

            if (addMuonlessPrimaryVertex_)
              myCand.addUserData("muonlessPV",Vertex(thePrimaryV.first));
            else
              myCand.addUserData("PVwithmuons",thePrimaryV.first);

            // lifetime using PV
            pvtx.SetXYZ(thePrimaryV.first.position().x(),thePrimaryV.first.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlphaXY = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            Measurement1D distXY = vdistXY.distance(myVertex->vertexState(), thePrimaryV.first);
            double ctauPV = distXY.value()*cosAlphaXY * myCand.mass()/pperp.Perp();
            GlobalError v1e = (myRecoVertex).error();
            GlobalError v2e = thePrimaryV.first.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*myCand.mass()/(pperp.Perp2());
            double lxy = distXY.value();
            double lxysig = lxy/distXY.error();

            pvtx3D.SetXYZ(thePrimaryV.first.position().x(),thePrimaryV.first.position().y(),thePrimaryV.first.position().z());
            TVector3 vdiff3D = vtx3D - pvtx3D;
            double cosAlpha3D = vdiff3D.Dot(pCand)/(vdiff3D.Mag()*pCand.Mag());
            Measurement1D dist3D = vdist3D.distance(myVertex->vertexState(), thePrimaryV.first);
            double l3d = dist3D.value();
            double l3dsig = l3d/dist3D.error();

            myCand.addUserFloat("ppdlPV",ctauPV);
            myCand.addUserFloat("ppdlErrPV",ctauErrPV);
            myCand.addUserFloat("cosAlphaXY",cosAlphaXY);
            myCand.addUserFloat("cosAlpha3D",cosAlpha3D);
            myCand.addUserFloat("l3d",l3d);
            myCand.addUserFloat("l3dsig",l3dsig);
            myCand.addUserFloat("lxy",lxy);
            myCand.addUserFloat("lxysig",lxysig);

            if (addCommonVertex_) {
              myCand.addUserData("commonVertex",myRecoVertex);
            }
          }
        } else {
          myCand.addUserFloat("vNChi2",-1);
          myCand.addUserFloat("vProb", -1);
          myCand.addUserFloat("PVndof", -1.);
          myCand.addUserFloat("pvw8", -1.);
          myCand.addUserFloat("pvw8Refit", -1.);
          myCand.addUserFloat("ppdlPV",-100);
          myCand.addUserFloat("ppdlErrPV",-100);
          myCand.addUserFloat("cosAlphaXY",-100);
          myCand.addUserFloat("cosAlpha3D",-100);
          myCand.addUserFloat("pvlip",-9999);
          myCand.addUserFloat("pvlipErr",-9999);
          myCand.addUserFloat("l3d",-100);
          myCand.addUserFloat("l3dsig",-100);
          myCand.addUserFloat("lxy",-100);
          myCand.addUserFloat("lxysig",-100);
          myCand.addUserFloat("ppdlBS",-100);
          myCand.addUserFloat("ppdlErrBS",-100);
          myCand.addUserFloat("Isolation", -1);
          myCand.addUserFloat("Isolation20", -1);
          myCand.addUserFloat("minDca", -1);
          myCand.addUserInt("Ntrk", 9999);
          myCand.addUserInt("Ntrkhp", 9999);
          myCand.addUserInt("Ntrk20", 9999);
          myCand.addUserInt("Ntrkhp20", 9999);
          myCand.addUserInt("Ntrk1sigma", 9999);
          myCand.addUserInt("Ntrk1sigmahp", 9999);
          myCand.addUserInt("Ntrk1sigma20", 9999);
          myCand.addUserInt("Ntrk1sigmahp20", 9999);
          myCand.addUserInt("Ntrk2sigma", 9999);
          myCand.addUserInt("Ntrk2sigmahp", 9999);
          myCand.addUserInt("Ntrk2sigma20", 9999);
          myCand.addUserInt("Ntrk2sigmahp20", 9999);
          myCand.addUserInt("Ntrk3sigma", 9999);
          myCand.addUserInt("Ntrk3sigmahp", 9999);
          myCand.addUserInt("Ntrk3sigma20", 9999);
          myCand.addUserInt("Ntrk3sigmahp20", 9999);
          myCand.addUserInt("countTksOfPV", 9999);
          myCand.addUserInt("countTksOfNoVtx", 9999);
          myCand.addUserFloat("sumPTPV", -1);
          myCand.addUserFloat("sumPTNoVtx", -1);
          myCand.addUserFloat("sumNdofPV",  -1);
          myCand.addUserFloat("sumNdofNoVtx", -1);
          myCand.addUserFloat("isoMu1", -1);
          myCand.addUserFloat("isoMu2", -1);
          myCand.addUserFloat("isoMu1_20", -1);
          myCand.addUserFloat("isoMu2_20", -1);
          myCand.addUserFloat("DCAXY", -1);
          myCand.addUserFloat("DCA", -1);
          myCand.addUserFloat("delta3d",-1);
          myCand.addUserFloat("delta3dErr",-1);
          myCand.addUserFloat("kinMass", -1.);
          myCand.addUserFloat("kinMassError", -1.);
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

//void Onia2MuMuPAT::buildMuonLessCollection(const edm::Event& iEvent,
//                                           const pat::Muon & mu1, const pat::Muon & mu2,
//                                           reco::TrackCollection & muonLess)
//{
//  edm::Handle<reco::TrackCollection> tkColl;
//  iEvent.getByLabel("generalTracks", tkColl);

//  muonLess.clear();
//  muonLess.reserve(tkColl->size());

//  TLorentzVector mu1P4=TLorentzVector(mu1.innerTrack()->px(),mu1.innerTrack()->py(),mu1.innerTrack()->pz(),mu1.energy());
//  TLorentzVector mu2P4=TLorentzVector(mu2.innerTrack()->px(),mu2.innerTrack()->py(),mu2.innerTrack()->pz(),mu2.energy());
//  TLorentzVector trackP4;

//  for( reco::TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk ) {
//    trackP4 = TLorentzVector(tk->px(), tk->py(), tk->pz(), sqrt((tk->p()*tk->p())+(trackMass_*trackMass_)));
//    if( (trackP4.DeltaR(mu1P4) < 1.e-5) || (trackP4.DeltaR(mu2P4) < 1.e-5) ) continue;
//    muonLess.push_back(*tk);
//  }
//}

bool Onia2MuMuPAT::searchForTheThirdTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                          pat::CompositeCandidate& Cand, const pat::Muon* mu1, const pat::Muon* mu2,
                                          reco::Vertex& PV, reco::Track& ktk, bool& closestPV)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  bool found=false;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  Handle<TrackCollection> tkColl;
  iEvent.getByLabel("generalTracks", tkColl);

  double tmpProb=-1, errorMass=-1;
  TransientVertex SVtx;
  Track track;
  // TrackRef tref;
  int charge=10;

  TLorentzVector mu1P4=TLorentzVector(mu1->innerTrack()->px(),mu1->innerTrack()->py(),mu1->innerTrack()->pz(),mu1->energy());
  TLorentzVector mu2P4=TLorentzVector(mu2->innerTrack()->px(),mu2->innerTrack()->py(),mu2->innerTrack()->pz(),mu2->energy());
  TLorentzVector trackP4;

  for (TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk) {

    trackP4=TLorentzVector(tk->px(),tk->py(),tk->pz(), sqrt((tk->p()*tk->p())+(trackMass_*trackMass_)));

    if (!tk->quality(TrackBase::highPurity)) continue;
    if (tk->pt() < minTrackPt_) continue;

    if ((trackP4.DeltaR(mu1P4) < 1.e-5) || (trackP4.DeltaR(mu2P4) < 1.e-5)) continue;

    TLorentzVector dimuP4=mu1P4+mu2P4;//TLorentzVector(mu1->px()+mu2->px(),mu1->py()+mu2->py(),mu1->pz()+mu2->pz(),mu1->energy()+mu2->energy());
    TLorentzVector totP4=dimuP4+trackP4;

    if (totP4.M() > diMuPlusTrackMassMax_ || totP4.M() < diMuPlusTrackMassMin_) continue;

    vector<TransientTrack> t_tks;
    track=Track(*tk);

    t_tks.push_back(theTTBuilder->build(mu1->track())); // pass the reco::Track, not  the reco::TrackRef (which can be transient)
    t_tks.push_back(theTTBuilder->build(mu2->track())); // otherwise the vertex will have transient refs inside.
    t_tks.push_back(theTTBuilder->build(track));

    TransientVertex tmpVtx = vtxFitter.vertex(t_tks);
    TransientVertex tmpaVtx = avtxFitter_.vertex(t_tks);

    if (!tmpVtx.isValid()) continue;
    if (!tmpaVtx.isValid()) continue;

    CachingVertex<5> aVtxForInvMass = avtxFitter_.vertex( t_tks );
    Measurement1D MassWErr = massCalculator.invariantMass( aVtxForInvMass, muMasses_ );

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
    trackP4=TLorentzVector(track.px(),track.py(),track.pz(), sqrt((track.p()*track.p())+(trackMass_*trackMass_)));
    LorentzVector vtk=LorentzVector(track.px(),track.py(),track.pz(),sqrt((track.p()*track.p())+(trackMass_*trackMass_)));
    LorentzVector jpsi = mu1->p4() + mu2->p4();
    Cand.setP4(jpsi+vtk);

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

    if (closestPV) {
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

    if (addMuonlessPrimaryVertex_) buildMuonlessPV(iEvent, iSetup, mu1, mu2, track, thePrimaryV);
    PV=Vertex(thePrimaryV);
    
    // count the number of high Purity tracks with pT > 900 MeV attached to the chosen vertex
    double vertexWeight = 0., sumPTPV = 0., sumPt_noVtx = 0.;
    int countTksOfPV = 0, Ntrk=0;
    double minDca = 9999.;
    double minDcaDR0p7 = 9999.;
    double minDcaDR0p7PtMin0p5 = 9999.;
    double minDcaDR0p7PtMin0p9 = 9999.;
    double minDcaPtMin0p5 = 9999.;
    double minDcaPtMin0p9 = 9999.;
    double sumNdof = 0.;
    double PVndof = thePrimaryV.ndof();
    Cand.addUserFloat("PVndof",(float) PVndof);
    Cand.addUserFloat("pvw8",(float) (thePrimaryV.ndof()+2)/(2*thePrimaryV.tracksSize()));

    // std::cout << "tracks in the PV = " << std::endl;
    for(reco::Vertex::trackRef_iterator itVtx = thePrimaryV.tracks_begin(); itVtx != thePrimaryV.tracks_end(); itVtx++){
      if(itVtx->isNonnull()){
        const reco::Track& trk = **itVtx;

        TLorentzVector tmpP4=TLorentzVector(trk.px(),trk.py(),trk.pz(),1.0);

        if ((tmpP4.DeltaR(mu1P4) < 1.e-5) || (tmpP4.DeltaR(mu2P4) < 1.e-5) || (tmpP4.DeltaR(trackP4) < 1.e-5)) continue;

        //if (trk.pt()==mu1->innerTrack()->pt() || trk.pt()==mu2->innerTrack()->pt() || trk.pt()==track.pt()) continue;

        //if(!track.quality(reco::TrackBase::highPurity)) continue;
        TransientTrack tt = theTTBuilder->build(trk);
        pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,SVtx);

        vertexWeight += thePrimaryV.trackWeight(*itVtx);
        if(track.pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;

        if( tkPVdist.first ) {
          if( tkPVdist.second.value()<minDca ) minDca = tkPVdist.second.value();
          if( track.pt() > 0.5 && tkPVdist.second.value()<minDcaPtMin0p5 ) minDcaPtMin0p5 = tkPVdist.second.value();
          if( track.pt() > 0.9 && tkPVdist.second.value()<minDcaPtMin0p9 ) minDcaPtMin0p9 = tkPVdist.second.value();
        }

        if(deltaR(track.eta(),track.phi(),Cand.eta(),Cand.phi()) < 0.7) {
          if( tkPVdist.first ) {
            if( tkPVdist.second.value()<minDcaDR0p7 ) minDcaDR0p7 = tkPVdist.second.value();
            if( track.pt() > 0.5 && tkPVdist.second.value()<minDcaDR0p7PtMin0p5 ) minDcaDR0p7PtMin0p5 = tkPVdist.second.value();
            if( track.pt() > 0.9 && tkPVdist.second.value()<minDcaDR0p7PtMin0p9 ) minDcaDR0p7PtMin0p9 = tkPVdist.second.value();
          }

          if(track.pt() < 0.9) continue; //reject all rejects from counting if less than 900 MeV
          sumPTPV += track.pt();
        }
        sumNdof += track.ndof();
        countTksOfPV++;
      }
    }
    
    // now tracks outside the PV fit
    int countTksOfNoVtx = 0;

    for (TrackCollection::const_iterator tk = tkColl->begin(); tk != tkColl->end(); ++tk) {
      bool hasVertex = false;
      for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv) {
        for( reco::Vertex::trackRef_iterator itVtx = itv->tracks_begin(); itVtx != itv->tracks_end(); ++itVtx) {
          if( itVtx->isNonnull() ) {
            if(deltaR(tk->eta(),tk->phi(),(*itVtx)->eta(),(*itVtx)->phi()) < 1.e-5) {
              hasVertex=true;
              bool isAlreadySelected=false;

              if ((deltaR(mu1->innerTrack()->eta(),mu1->innerTrack()->phi(),tk->eta(),tk->phi()) < 1.e-5) ||
                  (deltaR(mu2->innerTrack()->eta(),mu2->innerTrack()->phi(),tk->eta(),tk->phi()) < 1.e-5) ||
                  (deltaR(track.eta(),track.phi(),tk->eta(),tk->phi()) < 1.e-5)) isAlreadySelected=true;

              if (!hasVertex && !isAlreadySelected){
                TransientTrack tt = theTTBuilder->build(*tk);
                pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,SVtx);
                if(tk->pt() > 0.5 && tkPVdist.second.value()<0.03) ++Ntrk;

                if( tkPVdist.first ) {
                  if( tkPVdist.second.value()<minDca ) minDca = tkPVdist.second.value();
                  if( tk->pt() > 0.5 && tkPVdist.second.value()<minDcaPtMin0p5 ) minDcaPtMin0p5 = tkPVdist.second.value();
                  if( tk->pt() > 0.9 && tkPVdist.second.value()<minDcaPtMin0p9 ) minDcaPtMin0p9 = tkPVdist.second.value();
                }

                if(deltaR(tk->eta(),tk->phi(),Cand.eta(), Cand.phi()) < 0.7) {
                  if( tkPVdist.first ) {
                    if( tkPVdist.second.value()<minDcaDR0p7 ) minDcaDR0p7 = tkPVdist.second.value();
                    if( tk->pt() > 0.5 && tkPVdist.second.value()<minDcaDR0p7PtMin0p5 ) minDcaDR0p7PtMin0p5 = tkPVdist.second.value();
                    if( tk->pt() > 0.9 && tkPVdist.second.value()<minDcaDR0p7PtMin0p9 ) minDcaDR0p7PtMin0p9 = tkPVdist.second.value();
                  }

                  if(tk->pt()<0.9) continue;
                  if(tkPVdist.second.value() > 0.05) continue;
                  sumPt_noVtx += tk->pt();
                  countTksOfNoVtx++;
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
    Cand.addUserFloat("minDcaDr0p7", (float) minDcaDR0p7);
    Cand.addUserFloat("minDcaDr0p7Min0p5", (float) minDcaDR0p7PtMin0p5);
    Cand.addUserFloat("minDcaDr0p7Min0p9", (float) minDcaDR0p7PtMin0p9);
    Cand.addUserInt("Ntrk", Ntrk);
    Cand.addUserInt("countTksOfPV", countTksOfPV);
    Cand.addUserFloat("vertexWeight", (float) vertexWeight);
    Cand.addUserFloat("sumPTPV", (float) sumPTPV);
    Cand.addUserFloat("sumNdof", (float) sumNdof);

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

    // Longitudinal IP
    pair<bool,Measurement1D> pvlip = IPTools::signedDecayLength3D(tt,GlobalVector(0,0,1),thePrimaryV);
    if( pvlip.first ) {
      Cand.addUserFloat("pvlip",pvlip.second.value());
      Cand.addUserFloat("pvlipErr",pvlip.second.error());
    }

    //lifetime etc...

    TVector3 vtx,vtx3D;
    TVector3 pvtx,pvtx3D;
    VertexDistanceXY vdistXY;
    VertexDistance3D vdist3D;
    vtx.SetXYZ(SVtx.position().x(),SVtx.position().y(),0);
    TVector3 pperp(Cand.px(), Cand.py(), 0);
    AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

    vtx3D.SetXYZ(SVtx.position().x(),SVtx.position().y(),SVtx.position().z());
    TVector3 pCand((jpsi+vtk).px(), (jpsi+vtk).py(), (jpsi+vtk).pz());

    //using PV
    pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
    TVector3 vdiff = vtx - pvtx;
    double cosAlphaXY = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
    Measurement1D distXY = vdistXY.distance(Vertex(SVtx), thePrimaryV);
    double ctauPV = distXY.value()*cosAlphaXY * Cand.mass()/pperp.Perp();
    GlobalError v1e = (Vertex(SVtx)).error();
    GlobalError v2e = thePrimaryV.error();
    AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();
    double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*Cand.mass()/(pperp.Perp2());
    double lxy = distXY.value();
    double lxysig = lxy/distXY.error();

    pvtx3D.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),thePrimaryV.position().z());
    TVector3 vdiff3D = vtx3D - pvtx3D;
    double cosAlpha3D = vdiff3D.Dot(pCand)/(vdiff3D.Mag()*pCand.Mag());
    Measurement1D dist3D = vdist3D.distance(Vertex(SVtx), thePrimaryV);
    double l3d = dist3D.value();
    double l3dsig = l3d/dist3D.error();

    Cand.addUserFloat("ppdlPV",ctauPV);
    Cand.addUserFloat("ppdlErrPV",ctauErrPV);
    Cand.addUserFloat("cosAlphaXY",cosAlphaXY);
    Cand.addUserFloat("cosAlpha3D",cosAlpha3D);
    Cand.addUserFloat("l3d",l3d);
    Cand.addUserFloat("l3dsig",l3dsig);
    Cand.addUserFloat("lxy",lxy);
    Cand.addUserFloat("lxysig",lxysig);

  }else{
    Cand.addUserFloat("vNChi2",-1);
    Cand.addUserFloat("vProb", -1);
    Cand.addUserFloat("ppdlPV",-100);
    Cand.addUserFloat("pvw8", -1.);
    Cand.addUserFloat("ppdlErrPV",-100);
    Cand.addUserFloat("cosAlphaXY",-100);
    Cand.addUserFloat("PVndof", -1.);
    Cand.addUserFloat("cosAlpha3D",-100);
    Cand.addUserFloat("pvlip",-9999);
    Cand.addUserFloat("pvlipErr",-9999);
    Cand.addUserFloat("l3d",-100);
    Cand.addUserFloat("l3dsig",-100);
    Cand.addUserFloat("lxy",-100);
    Cand.addUserFloat("lxysig",-100);
    Cand.addUserFloat("ppdlBS",-100);
    Cand.addUserFloat("ppdlErrBS",-100);
    Cand.addUserFloat("Isolation", -1);
    Cand.addUserFloat("minDca", -1);
    Cand.addUserFloat("minDcaDr0p7", -1);
    Cand.addUserFloat("minDcaDr0p7Min0p5", -1);
    Cand.addUserFloat("minDcaDr0p7Min0p9", -1);
    Cand.addUserInt("Ntrk", 9999);
    Cand.addUserInt("countTksOfPV", 9999);
    Cand.addUserFloat("vertexWeight", -1);
    Cand.addUserFloat("sumPTPV", -1);
    Cand.addUserFloat("delta3d",-1);
    Cand.addUserFloat("delta3dErr",-1);
    Cand.addUserFloat("MassErr",-1);
    Cand.addUserFloat("DCAXY", -1);
    Cand.addUserFloat("DCA", -1);
    Cand.addUserFloat("sumNdof", -1);
  }
  return found;
}

// ---- method to build the PV without the candidate tracks included ----------
std::pair<reco::Vertex, reco::Vertex> Onia2MuMuPAT::buildMuonlessPV(const reco::Muon * rmu1, const reco::Muon * rmu2,
                                                                    reco::TrackCollection & muonLess, const reco::Vertex & thePrimaryV,
                                                                    const edm::ESHandle<TransientTrackBuilder> & theTTBuilder,
                                                                    const reco::BeamSpot &bs)
{
  muonLess.clear();
  std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV.tracks_begin();
  for( ; itPVtrack != thePrimaryV.tracks_end(); ++itPVtrack ) {
    if (itPVtrack->isNonnull()) {
      if( itPVtrack->key() == rmu1->track().key() ) continue;
      if( itPVtrack->key() == rmu2->track().key() ) continue;
      muonLess.push_back(**itPVtrack);
    }
  }
  if (muonLess.size()>1) {
    std::vector<reco::TransientTrack> t_tks;
    t_tks.reserve(muonLess.size());
    for (auto it = muonLess.begin(); it != muonLess.end(); ++it) {
      t_tks.push_back((*theTTBuilder).build(*it));
      t_tks.back().setBeamSpot(bs);
    }
    TransientVertex tpvs = avtxFitter_.vertex(t_tks, bs);
    if( tpvs.isValid() ) {
      reco::Vertex muonLessPV = reco::Vertex(tpvs);
      return std::make_pair(muonLessPV, thePrimaryV);
    }
  }
  std::cout << "PV refit failed, returning original PV" << std::endl;
  return std::make_pair(thePrimaryV, thePrimaryV);
}

// Used for the dimuon+track
void Onia2MuMuPAT::buildMuonlessPV(const edm::Event& iEvent, const edm::EventSetup& iSetup, const pat::Muon * rmu1,const pat::Muon * rmu2, const reco::Track & track, reco::Vertex & thePrimaryV)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  TrackCollection muonLess;

  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel(thePVs_, priVtxs);

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel(thebeamspot_,theBeamSpot);
  
  // Forcing to use the adaptive vertex fitter
  VertexReProducer revertex(priVtxs, iEvent, true);
  Handle<TrackCollection> pvtracks;
  iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
  if( !pvtracks.isValid()) {cout << "PV tracks not valid!!" << endl; return;}

  Handle<BeamSpot>        pvbeamspot;
  iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
  if (pvbeamspot.id() != theBeamSpot.id()) edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

  //to see if the track is in pv fit
  bool tkIn=false;

  for (TrackCollection::const_iterator tk = pvtracks->begin(); tk != pvtracks->end(); ++tk) {
    double dR=deltaR(tk->eta(),tk->phi(),track.eta(),track.phi());
    if (dR<1.e-5) tkIn=true;
  }

  if (rmu1->track().id() == pvtracks.id() || rmu2->track().id() == pvtracks.id() || tkIn) {
    if( thePrimaryV.hasRefittedTracks() ) {
      // Need to go back to the original tracks before taking the key
      std::vector<reco::Track>::const_iterator itRefittedTrack = thePrimaryV.refittedTracks().begin();
      std::vector<reco::Track>::const_iterator refittedTracksEnd = thePrimaryV.refittedTracks().end();
      for( ; itRefittedTrack != refittedTracksEnd; ++itRefittedTrack ) {
        if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu1->track().key() ) continue;
        if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu2->track().key() ) continue;
        if( deltaR(thePrimaryV.originalTrack(*itRefittedTrack)->eta(),thePrimaryV.originalTrack(*itRefittedTrack)->phi(),track.eta(),track.phi()) < 1.e-5) continue;
        //if( thePrimaryV.originalTrack(*itRefittedTrack).key() == track.key() ) continue;
        muonLess.push_back(*(thePrimaryV.originalTrack(*itRefittedTrack)));
      }
    }
    else {
      std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV.tracks_begin();
      for( ; itPVtrack != thePrimaryV.tracks_end(); ++itPVtrack ) if (itPVtrack->isNonnull()) {
        if( itPVtrack->key() == rmu1->track().key() ) continue;
        if( itPVtrack->key() == rmu2->track().key() ) continue;
        if(deltaR((*itPVtrack)->eta(),(*itPVtrack)->phi(),track.eta(),track.phi())< 1.e-5) continue;
        muonLess.push_back(**itPVtrack);
      }
    }
    if (muonLess.size()>1 && muonLess.size() < thePrimaryV.tracksSize()){
      vector<TransientVertex> pvs;
      pvs = revertex.makeVertices(muonLess, *pvbeamspot, iSetup) ;
      if (!pvs.empty()) {
        Vertex muonLessPV = Vertex(pvs.front());
        thePrimaryV = muonLessPV;
      }
    }
  }
}

int Onia2MuMuPAT::findVertexId(const std::pair<reco::Vertex, reco::Vertex> & theOriginalPV,
                               const std::vector<std::pair<reco::Vertex, reco::Vertex> >&priVtxs,
                               const reco::Track & track, const double & maxDeltaR)
{
  // Check if the track belongs to the PV associated to the candidate
  for(auto PVtk = theOriginalPV.second.tracks_begin(); PVtk != theOriginalPV.second.tracks_end(); ++PVtk) {
    if( PVtk->isNonnull() ) {
      if( PVtk->key() == track.extra().key() ) {
        // if( deltaR(track.eta(), track.phi(), (*PVtk)->eta(), (*PVtk)->phi()) < maxDeltaR ) return 1;
        return 1;
      }
    }
  }
  // Check if the track is associated to any other PV
  for(auto itv = priVtxs.begin(), itvend = priVtxs.end(); itv != itvend; ++itv){
    for( auto vtxTk = itv->second.tracks_begin(); vtxTk != itv->second.tracks_end(); ++vtxTk) {
      if( vtxTk->isNonnull() ) {
        if( vtxTk->key() == track.extra().key() ) {
          return 2;
        }
        // if( deltaR(track.eta(), track.phi(), (*vtxTk)->eta(), (*vtxTk)->phi()) < maxDeltaR ) return 2;
      }
    }
  }
  // Track is not associated to any PV
  return 0;
}

float Onia2MuMuPAT::computeDcaXY(const TrajectoryStateClosestToPoint & tt1, const TrajectoryStateClosestToPoint & tt2)
{
  float dcaxy = 1E20;
  if (tt1.isValid() && tt2.isValid()) {
    ClosestApproachInRPhi cApp;
    if( cApp.calculate(tt1.theState(), tt2.theState()) ) dcaxy = cApp.distance();
  }
  return dcaxy;
}

float Onia2MuMuPAT::computeDca(const TrajectoryStateClosestToPoint & tt1, const TrajectoryStateClosestToPoint & tt2)
{
  float dca = 1E20;
  if (tt1.isValid() && tt2.isValid()) {
    TwoTrackMinimumDistance ttmd;
    if( ttmd.calculate(tt1.theState(), tt2.theState()) ) dca = ttmd.distance();
  }
  return dca;
}

void Onia2MuMuPAT::setP4AndCharge(RefCountedKinematicTree & vertexFitTree, pat::Muon & mu,
                                  const edm::ESHandle<TransientTrackBuilder> &theTTBuilder) const
{
  reco::TransientTrack refittedTrack = theTTBuilder->build(vertexFitTree->currentParticle()->currentState().freeTrajectoryState());
  math::XYZTLorentzVector mup4(fromPtEtaPhiToPxPyPz(refittedTrack.track().pt(),
                                                    refittedTrack.track().eta(),
                                                    refittedTrack.track().phi(),
                                                    muMasses_.at(0)));
  mu.setP4(mup4);
  mu.setCharge(vertexFitTree->currentParticle()->currentState().particleCharge());
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
