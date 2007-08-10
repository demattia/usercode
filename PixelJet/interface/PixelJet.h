#ifndef AnalysisExamples_PixelJet_h
#define AnalysisExamples_PixelJet_h

/////////////////////////////////////////////////////
//
// PixelJet class
// Author M. De Mattia - 6/8/2007
//
// Declares a PixelJet, a jet formed of pixeltracks.
//
/////////////////////////////////////////////////////

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/RefVector.h"

// Pixel tracks (collection type: recotrack)
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <vector>
#include <boost/cstdint.hpp>

class PixelJet{
 public:

  PixelJet() {

    // Initializations per PixelJet
    // ----------------------------

    pixelTracksNumber_ = 0;
    eta_ = 0.;
    phi_ = 0.;
    pt_ = 0.;
    z_ = 0.;
  }

  void PixelTrackRef( reco::TrackRef pixeltrack_ref ) {
    pt_ += pixeltrack_ref->pt();
    eta_ += pixeltrack_ref->eta()*pixeltrack_ref->pt();
    phi_ += pixeltrack_ref->phi()*pixeltrack_ref->pt();
    z_ += ( pixeltrack_ref->vertex().z() )*pixeltrack_ref->pt();
    vecRefPixelTrack.push_back(pixeltrack_ref);
    ++pixelTracksNumber_;
  }

  const edm::RefVector<std::vector<reco::Track> > GetPixelTrackRefVec() const {
    return vecRefPixelTrack;
  }

  void Close() {
    if ( pixelTracksNumber_ != 0 ) {
      eta_ = eta_/pt_;
      phi_ = phi_/pt_;
      z_ = z_/pt_;
    }
    else {
      std::cout << "Error: this pixel jet has no tracks" << std::endl;
    }
  }

  double pt() const {
    return pt_;
  }
  double eta() const {
    return eta_;
  }
  double phi() const {
    return phi_;
  }
  double z() const {
    return z_;
  }
  int NumTk() const {
    return pixelTracksNumber_;
  }

 private:
  // Vector of Refs to PixelTracks
  edm::RefVector<std::vector<reco::Track> > vecRefPixelTrack;

  int pixelTracksNumber_;
  double pt_;
  double eta_;
  double phi_;
  double z_;
};

// PixelJet collection typedef
typedef std::vector<PixelJet> PixelJetCollection;
typedef edm::Ref<PixelJetCollection> PixelJetRef;
typedef edm::RefProd<PixelJetCollection> PixelJetRefProd;
typedef edm::RefVector<PixelJetCollection> PixelJetRefVector;

#endif //AnalysiExamples_PixelJet_h
