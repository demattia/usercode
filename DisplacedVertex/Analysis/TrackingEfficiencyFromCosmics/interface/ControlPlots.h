#ifndef CONTROLPLOTS_H
#define CONTROLPLOTS_H

#include <vector>
#include <TH1F.h>
#include <sstream>

#include "Analysis/TrackingEfficiencyFromCosmics/interface/Utils.h"
#include "Analysis/TrackingEfficiencyFromCosmics/interface/SmartPropagatorWithIP.h"

/**
  * Fill a set of control plots for all the tracks in a given collection.
  */
class ControlPlots
{
public:
  ControlPlots( edm::Service<TFileService> & fileService, const TString & name )
  {
    hPt_ =                    utils::bookHistogram(fileService, name, "pt", "P_{T}", "[GeV/c]", 500, 0, 500);
    hEta_ =                   utils::bookHistogram(fileService, name, "eta", "#eta", "", 500, -3, 3);
    hPhi_ =                   utils::bookHistogram(fileService, name, "phi", "#phi", "", 500, -3.2, 3.2);
    hNhits_ =                 utils::bookHistogram(fileService, name, "Nhits", "# hits", "", 100, 0, 100);
    hNValidHits_ =            utils::bookHistogram(fileService, name, "NValidHits", "# valid hits", "", 100, 0, 100);
    hNValidPlusInvalidHits_ = utils::bookHistogram(fileService, name, "NValidPlusInvalidHits", "# (valid + invalid )hits", "", 100, 0, 100);
    hInnermostHitRadius_ =    utils::bookHistogram(fileService, name, "innermostHitRadius", "innermost hit radius", "cm", 500, 0, 100);
    hInnermostHitZ_ =         utils::bookHistogram(fileService, name, "innermostHitZ", "innermost hit z", "cm", 500, 0, 100);
    hDxy_ =                   utils::bookHistogram(fileService, name, "d_0", "|d_{0}|", "cm", 500, 0, 500);
    hDxyErr_ =                utils::bookHistogram(fileService, name, "d_0_error", "#sigma(d_{0})", "cm", 400, -20, 20);
    hDz_ =                    utils::bookHistogram(fileService, name, "d_z", "|d_{z}|", "cm", 500, 0, 50);
    hDzErr_ =                 utils::bookHistogram(fileService, name, "d_z_error", "#sigma(d_{z})", "cm", 400, -20, 20);
    hChi2_ =                  utils::bookHistogram(fileService, name, "chi2", "#chi^{2}", "", 500, 0, 100);
    hReferencePointRadius_ =  utils::bookHistogram(fileService, name, "RefPointRadius", "radius of the reference point", "", 500, -50, 50);
    hReferencePointZ_ =       utils::bookHistogram(fileService, name, "RefPointZ", "z of the reference point", "", 500, -50, 50);
    hVertexX_ =               utils::bookHistogram(fileService, name, "VertexX", "x of the vertex of this track", "", 500, -100, 100);
    hVertexY_ =               utils::bookHistogram(fileService, name, "VertexY", "y of the vertex of this track", "", 500, -100, 100);
    hVertexZ_ =               utils::bookHistogram(fileService, name, "VertexZ", "z of the vertex of this track", "", 500, -100, 100);

    // 2D histograms
    hDxyErrVsNValidHit_ =     utils::bookHistogram(fileService, name, "d_0_errorVsNValidHit", "# valid hits", "", 80, 0, 80, "#sigma(d_{0})", "cm", 400, 0, 4);
    hDxyErrVsPt_ =            utils::bookHistogram(fileService, name, "d_0_errorVsPt", "p_{T}", "[GeV/c]", 500, 0, 500, "#sigma(d_{0})", "cm", 400, 0, 4);
    hDxyErrVsDxy_ =           utils::bookHistogram(fileService, name, "d_0_errorVsDz", "d_{0}", "cm", 500, 0, 500, "#sigma(d_{0})", "cm", 400, 0, 4);

    hDzErrVsNValidHit_ =      utils::bookHistogram(fileService, name, "d_z_errorVsNValidHit", "# valid hits", "", 80, 0, 80, "#sigma(d_{z})", "cm", 400, 0, 4);
    hDzErrVsPt_ =             utils::bookHistogram(fileService, name, "d_z_errorVsPt", "p_{T}", "[GeV/c]", 500, 0, 500, "#sigma(d_{z})", "cm", 400, 0, 4);
    hDzErrVsDz_ =             utils::bookHistogram(fileService, name, "d_z_errorVsDz", "d_{z}", "cm", 500, 0, 50, "#sigma(d_{z})", "cm", 400, 0, 4);
  }

  enum propType
  {
    INSIDETKVOL,
    INSIDEOUT,
    OUTSIDEIN
  };

  /// Used for genParticles
  void fillControlPlots(const SmartPropagatorWithIP::IP & ip, const reco::LeafCandidate::Point & vertex)
  {
    hPt_->Fill(ip.pt);
    hEta_->Fill(ip.eta);
    hPhi_->Fill(ip.phi);
    hDxy_->Fill(fabs(ip.dxyValue));
    hDxyErr_->Fill(ip.dxyError);
    hDxyErrVsPt_->Fill(ip.pt, ip.dxyError);
    hDxyErrVsDxy_->Fill(fabs(ip.dxyValue), ip.dxyError);
    hDz_->Fill(fabs(ip.dzValue));
    hDzErr_->Fill(ip.dzError);
    hDzErrVsPt_->Fill(ip.pt, ip.dzError);
    hDzErrVsDz_->Fill(fabs(ip.dzValue), ip.dzError);
    hVertexX_->Fill(vertex.x());
    hVertexY_->Fill(vertex.y());
    hVertexZ_->Fill(vertex.z());
  }

  /**
   * Propagator type:
   * 0 is inside the tracker volume
   * 1 is inside-out
   * 2 is outside-in
   */
  template <class T>
  void fillControlPlots(const T & track, const SmartPropagatorWithIP * smartPropIP = 0, const int type = 0)
  {
    hNhits_->Fill(track.recHitsSize());
    hNValidHits_->Fill(track.found());
    hNValidPlusInvalidHits_->Fill(track.found() + track.lost());
    math::XYZPoint innermostHitPosition(track.innerPosition());
    hInnermostHitRadius_->Fill(innermostHitPosition.r());
    hInnermostHitZ_->Fill(innermostHitPosition.z());
    hChi2_->Fill(track.normalizedChi2());
    hReferencePointRadius_->Fill(sqrt(track.referencePoint().x()+track.referencePoint().y()));
    hReferencePointZ_->Fill(track.referencePoint().z());
    hVertexX_->Fill(track.vertex().x());
    hVertexY_->Fill(track.vertex().y());
    hVertexZ_->Fill(track.vertex().z());

    SmartPropagatorWithIP::IP ip(track.pt(), track.eta(), track.phi(), track.dxy(), track.dxyError(), track.dz(), track.dzError());

    if( smartPropIP != 0 ) {
      if( type == INSIDETKVOL ) {
        std::cout << "propagating inside tk volume" << std::endl;
        ip = smartPropIP->computeImpactParametersInsideTkVol(track, GlobalPoint(0,0,0));
      }
      else if( type == INSIDEOUT ) {
        std::cout << "propagating inside-out" << std::endl;
        ip = smartPropIP->computeImpactParametersInsideOutTkVol(track, GlobalPoint(0,0,0));
      }
      else if( type == OUTSIDEIN ) {
        std::cout << "propagating outside-in" << std::endl;
        ip = smartPropIP->computeImpactParametersOutsideInTkVol(track, GlobalPoint(0,0,0));
      }
      else {
        std::cout << "Unknown propagation type: " << type << " reading IP from the track" << std::endl;
      }
    }
    hPt_->Fill(ip.pt);
    hEta_->Fill(ip.eta);
    hPhi_->Fill(ip.phi);
    hDxy_->Fill(fabs(ip.dxyValue));
    hDxyErr_->Fill(ip.dxyError);
    hDxyErrVsNValidHit_->Fill(track.found(), ip.dxyError);
    hDxyErrVsPt_->Fill(ip.pt, ip.dxyError);
    hDxyErrVsDxy_->Fill(fabs(ip.dxyValue), ip.dxyError);
    hDz_->Fill(fabs(ip.dzValue));
    hDzErr_->Fill(ip.dzError);
    hDzErrVsNValidHit_->Fill(track.found(), ip.dzError);
    hDzErrVsPt_->Fill(ip.pt, ip.dzError);
    hDzErrVsDz_->Fill(fabs(ip.dzValue), ip.dzError);
  }

  template <class T>
  void fillControlPlots(const std::vector<T> & collection, const SmartPropagatorWithIP * smartPropIP = 0, const int type = 0)
  {
    typename std::vector<T>::const_iterator it = collection.begin();
    for( ; it != collection.end(); ++it ) {      
      fillControlPlots(*it, smartPropIP, type);
    }
  }

protected:
  TH1F *hPt_, *hEta_, *hPhi_, *hNhits_, *hNValidHits_;
  TH1F *hNValidPlusInvalidHits_, *hInnermostHitRadius_, *hInnermostHitZ_;
  TH1F *hDxy_, *hDxyErr_, *hDz_, *hDzErr_, *hChi2_, *hReferencePointRadius_, *hReferencePointZ_;
  TH1F *hVertexX_, *hVertexY_, *hVertexZ_;


  TH2F *hDxyErrVsNValidHit_, *hDxyErrVsPt_, *hDxyErrVsDxy_;
  TH2F *hDzErrVsNValidHit_, *hDzErrVsPt_, *hDzErrVsDz_;
};

#endif // CONTROLPLOTS_H
