#ifndef CONTROLDELTAPLOTS_H
#define CONTROLDELTAPLOTS_H

#include <vector>
#include <TH1F.h>
#include "DataFormats/Math/interface/deltaPhi.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/Utils.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"

/**
  * Fill a set of control plots for the differences between two tracks. 
  */
class ControlDeltaPlots
{
public:
  ControlDeltaPlots( edm::Service<TFileService> & fileService, const TString & name, const int sign = 1 ) :
    sign_(sign)
  {
    hDeltaPt_ =       utils::bookHistogram(fileService, name, "DeltaPt", "#Delta P_{T}", "", 100, -50, 50);
    hDeltaPtOverPt_ = utils::bookHistogram(fileService, name, "DeltaPtOverPt", "#Delta P_{T}/P_{T}", "[GeV/c]", 100, -5, 5);
    hDeltaEta_ =      utils::bookHistogram(fileService, name, "DeltaEta", "#Delta #eta", "", 100, -6, 6);
    hDeltaPhi_ =      utils::bookHistogram(fileService, name, "DeltaPhi", "#Delta #phi", "", 100, -3.4, 3.4);
    hDeltaDxy_ =      utils::bookHistogram(fileService, name, "DeltaDxy", "#Delta d_{0}", "cm", 100, -100, 100);
    hDeltaDz_ =       utils::bookHistogram(fileService, name, "DeltaDz", "#Delta d_{z}", "cm", 100, -100, 100);
    hDeltaFabsDxy_ =  utils::bookHistogram(fileService, name, "DeltaFabsDxy", "#Delta |d_{0}|", "cm", 100, -100, 100);
    hDeltaFabsDz_ =   utils::bookHistogram(fileService, name, "DeltaFabsDz", "#Delta |d_{z}|", "cm", 100, -100, 100);
    hDeltaVertexX_ =  utils::bookHistogram(fileService, name, "DeltaVertexX", "#Delta(Vertex X)", "cm", 120, -60, 60);
    hDeltaVertexY_ =  utils::bookHistogram(fileService, name, "DeltaVertexY", "#Delta(Vertex Y)", "cm", 120, -60, 60);
    hDeltaVertexZ_ =  utils::bookHistogram(fileService, name, "DeltaVertexZ", "#Delta(Vertex Z)", "cm", 200, -100, 100);
  }
  template <class T1, class T2>
  void fillControlPlots(const T1 & track1, const SmartPropagatorWithIP::IP & ip1,
                        const T2 & track2, const SmartPropagatorWithIP::IP & ip2)
  {
    hDeltaPt_->Fill(ip1.pt - ip2.pt);
    if( ip2.pt != 0 ) hDeltaPtOverPt_->Fill((ip1.pt - ip2.pt)/ip2.pt);
    // With sign = -1 back-to-back tracks have DeltaEta = 0
    hDeltaEta_->Fill(ip1.eta - sign_*ip2.eta);
    hDeltaPhi_->Fill(reco::deltaPhi(ip1.phi, ip2.phi));
    // With sign = -1 back-to-back tracks have deltaDxy = 0
    hDeltaDxy_->Fill(ip1.dxyValue - sign_*ip2.dxyValue);
    hDeltaDz_->Fill(ip1.dzValue - ip2.dzValue);
    hDeltaFabsDxy_->Fill(fabs(ip1.dxyValue) - fabs(ip2.dxyValue));
    hDeltaFabsDz_->Fill(fabs(ip1.dzValue) - fabs(ip2.dzValue));

    hDeltaVertexX_->Fill(track1.vertex().x() - track2.vertex().x());
    hDeltaVertexY_->Fill(track1.vertex().y() - track2.vertex().y());
    hDeltaVertexZ_->Fill(track1.vertex().z() - track2.vertex().z());
  }

  template <class T1, class T2>
  void fillControlPlots(const T1 & track1, const double & tk1Dxy, const double & tk1Dz,
                        const T2 & track2, const double & tk2Dxy, const double & tk2Dz)
  {
    fillControlPlots(track1, SmartPropagatorWithIP::IP(track1.pt(), track1.eta(), track1.phi(), tk1Dxy, 0., tk1Dz, 0.),
                     track2, SmartPropagatorWithIP::IP(track1.pt(), track1.eta(), track1.phi(), tk1Dxy, 0., tk1Dz, 0.));
  }

  template <class T1, class T2>
  inline void fillControlPlots(const T1 & track1, const T2 & track2)
  {
    fillControlPlots(track1, track1.dxy(), track1.dz(), track2, track2.dxy(), track2.dz());
  }

protected:
  TH1F *hDeltaPt_, *hDeltaPtOverPt_, *hDeltaEta_, *hDeltaPhi_, *hDeltaDxy_, *hDeltaDz_;
  TH1F *hDeltaVertexX_, *hDeltaVertexY_, *hDeltaVertexZ_;
  TH1F *hDeltaFabsDxy_, *hDeltaFabsDz_;
  int sign_;
};

#endif // CONTROLPLOTS_H
