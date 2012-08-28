#ifndef PseudoLepton_h
#define PseudoLepton_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class PseudoLepton;
typedef std::vector<PseudoLepton>              PseudoLeptonCollection;
typedef edm::Ref<PseudoLeptonCollection>       PseudoLeptonRef;
typedef edm::RefVector<PseudoLeptonCollection> PseudoLeptonRefVector;

class PseudoLepton : public reco::Track {

 private:

  reco::Particle::LorentzVector p4_;
  const reco::GenParticle* genPart_;
  unsigned genPartSize_;
  bool isStandAloneMuon_;
  bool isGlobalMuon_;
  bool isTrackerMuon_;
  float trackIso_;
  float ecalIso_;
  float hcalIso_;
  double normChi2_;
  int numberOfValidTrackerHits_;
  int numberOfValidPixelHits_;
  int numberOfValidMuonHits_;
  int numberOfMatchedStations_;

 public:

  PseudoLepton();
  PseudoLepton(const reco::Track & aTrack);
  virtual ~PseudoLepton();
  
  /// required reimplementation of the Candidate's clone method
  virtual PseudoLepton * clone() const { return new PseudoLepton(*this); }

  /// stuff that we need that is not in reco::Track
  const reco::Particle::LorentzVector& p4() const { return p4_; };
  const reco::GenParticle* genParticle(unsigned i=0) const { if (genPartSize_) return genPart_; else return 0; };
  unsigned genParticlesSize() const { return genPartSize_; };
  bool isTrackerMuon() const { return isTrackerMuon_; };
  bool isGlobalMuon() const { return isGlobalMuon_; };
  bool isStandAloneMuon() const { return isStandAloneMuon_; };

  void isGlobalMuon(const bool isMuon) { isTrackerMuon_=isMuon; };
  void isStandAloneMuon(const bool isMuon) { isStandAloneMuon_=isMuon; };
  void isTrackerMuon(const bool isMuon) { isGlobalMuon_=isMuon; };
  void setGenParticle(const reco::GenParticle& gen) { genPart_=&gen; genPartSize_=1; };

  void setTrackIso( const float & trackIso ) { trackIso_ = trackIso; }
  void setEcalIso( const float & ecalIso ) { ecalIso_ = ecalIso; }
  void setHcalIso( const float & hcalIso ) { hcalIso_ = hcalIso; }
  float trackIso() const { return trackIso_; }
  float ecalIso() const { return ecalIso_; }
  float hcalIso() const { return hcalIso_; }

  void setNormChi2( const double & normChi2 ) { normChi2_ = normChi2; }
  void setNumberOfValidTrackerHits( const int numberOfValidTrackerHits ) { numberOfValidTrackerHits_ = numberOfValidTrackerHits; }
  void setNumberOfValidPixelHits( const int numberOfValidPixelHits ) { numberOfValidPixelHits_ = numberOfValidPixelHits; }
  void setNumberOfValidMuonHits( const int numberOfValidMuonHits ) { numberOfValidMuonHits_ = numberOfValidMuonHits; }
  void setNumberOfMatchedStations( const int numberOfMatchedStations ) { numberOfMatchedStations_ = numberOfMatchedStations; }

  double normChi2() const { return normChi2_; }
  int numberOfValidTrackerHits() const { return numberOfValidTrackerHits_; }
  int numberOfValidPixelHits() const { return numberOfValidPixelHits_; }
  int numberOfValidMuonHits() const { return numberOfValidMuonHits_; }
  int numberOfMatchedStations() const { return numberOfMatchedStations_; }
};



#endif
