#ifndef CONTROLDELTAPLOTS_H
#define CONTROLDELTAPLOTS_H

#include <vector>
#include <TH1F.h>
#include "DataFormats/Math/interface/deltaPhi.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/Utils.h"

/**
  * Fill a set of control plots for the differences between two tracks. 
  */
class ControlDeltaPlots
{
public:
  ControlDeltaPlots( edm::Service<TFileService> & fileService, const TString & name, const int sign = 1 ) :
    sign_(sign)
  {
    hDeltaPt_ =      utils::bookHistogram(fileService, name, "DeltaPt", "#Delta P_{T}", "[GeV/c]", 400, -50, 50);
    hDeltaEta_ =     utils::bookHistogram(fileService, name, "DeltaEta", "#Delta #eta", "", 500, -6, 6);
    hDeltaPhi_ =     utils::bookHistogram(fileService, name, "DeltaPhi", "#Delta #phi", "", 200, -3.4, 3.4);
    hDeltaDxy_ =     utils::bookHistogram(fileService, name, "DeltaDxy", "#Delta d_{0}", "cm", 1000, -50, 50);
    hDeltaDz_ =      utils::bookHistogram(fileService, name, "DeltaDz", "#Delta d_{z}", "cm", 1000, -50, 50);
    hDeltaFabsDxy_ = utils::bookHistogram(fileService, name, "DeltaFabsDxy", "#Delta |d_{0}|", "cm", 100, -50, 50);
    hDeltaFabsDz_ =  utils::bookHistogram(fileService, name, "DeltaFabsDz", "#Delta |d_{z}|", "cm", 100, -50, 50);
    hDeltaVertexX_ = utils::bookHistogram(fileService, name, "DeltaVertexX", "#Delta(Vertex X)", "cm", 120, -60, 60);
    hDeltaVertexY_ = utils::bookHistogram(fileService, name, "DeltaVertexY", "#Delta(Vertex Y)", "cm", 120, -60, 60);
    hDeltaVertexZ_ = utils::bookHistogram(fileService, name, "DeltaVertexZ", "#Delta(Vertex Z)", "cm", 200, -100, 100);
  }
  template <class T1, class T2>
    void fillControlPlots(const T1 & track1, const T2 & track2)
  {
    hDeltaPt_->Fill(track1.pt() - track2.pt());
    // With sign = -1 back-to-back tracks have DeltaEta = 0
    hDeltaEta_->Fill(track1.eta() - sign_*track2.eta());
    hDeltaPhi_->Fill(reco::deltaPhi(track1.phi(), track2.phi()));
    // With sign = -1 back-to-back tracks have deltaDxy = 0
    hDeltaDxy_->Fill(track1.dxy() - sign_*track2.dxy());
    hDeltaDz_->Fill(track1.dz() - track2.dz());
    hDeltaFabsDxy_->Fill(fabs(track1.dxy()) - fabs(track2.dxy()));
    hDeltaFabsDz_->Fill(fabs(track1.dz()) - fabs(track2.dz()));

    hDeltaVertexX_->Fill(track1.vertex().x() - track2.vertex().x());
    hDeltaVertexY_->Fill(track1.vertex().y() - track2.vertex().y());
    hDeltaVertexZ_->Fill(track1.vertex().z() - track2.vertex().z());
  }

protected:
  TH1F *hDeltaPt_, *hDeltaEta_, *hDeltaPhi_, *hDeltaDxy_, *hDeltaDz_;
  TH1F *hDeltaVertexX_, *hDeltaVertexY_, *hDeltaVertexZ_;
  TH1F *hDeltaFabsDxy_, *hDeltaFabsDz_;
  int sign_;
};

#endif // CONTROLPLOTS_H
