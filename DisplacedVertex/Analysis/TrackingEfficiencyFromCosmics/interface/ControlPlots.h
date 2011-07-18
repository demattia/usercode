#ifndef CONTROLPLOTS_H
#define CONTROLPLOTS_H

#include <vector>
#include <TH1F.h>
#include <sstream>

#include "Analysis/TrackingEfficiencyFromCosmics/interface/Utils.h"

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
    hDz_ =                    utils::bookHistogram(fileService, name, "d_z", "|d_{z}|", "cm", 500, 0, 50);
    hChi2_ =                  utils::bookHistogram(fileService, name, "chi2", "#chi^{2}", "", 500, 0, 100);
    hReferencePointRadius_ =  utils::bookHistogram(fileService, name, "RefPointRadius", "radius of the reference point", "", 500, -50, 50);
    hReferencePointZ_ =       utils::bookHistogram(fileService, name, "RefPointZ", "z of the reference point", "", 500, -50, 50);
    hVertexX_ =               utils::bookHistogram(fileService, name, "VertexX", "x of the vertex of this track", "", 500, -100, 100);
    hVertexY_ =               utils::bookHistogram(fileService, name, "VertexY", "y of the vertex of this track", "", 500, -100, 100);
    hVertexZ_ =               utils::bookHistogram(fileService, name, "VertexZ", "z of the vertex of this track", "", 500, -100, 100);
  }

  template <class T>
  void fillControlPlots(const std::vector<T> & collection)
  {
    typename std::vector<T>::const_iterator it = collection.begin();
    for( ; it != collection.end(); ++it ) {
      hPt_->Fill(it->pt());
      hEta_->Fill(it->eta());
      hPhi_->Fill(it->phi());
      hNhits_->Fill(it->recHitsSize());
      hNValidHits_->Fill(it->found());
      hNValidPlusInvalidHits_->Fill(it->found() + it->lost());
      math::XYZPoint innermostHitPosition(it->innerPosition());
      hInnermostHitRadius_->Fill(innermostHitPosition.r());
      hInnermostHitZ_->Fill(innermostHitPosition.z());
      hDxy_->Fill(fabs(it->dxy()));
      hDz_->Fill(fabs(it->dz()));
      hChi2_->Fill(it->normalizedChi2());
      hReferencePointRadius_->Fill(sqrt(it->referencePoint().x()+it->referencePoint().y()));
      hReferencePointZ_->Fill(it->referencePoint().z());
      hVertexX_->Fill(it->vertex().x());
      hVertexY_->Fill(it->vertex().y());
      hVertexZ_->Fill(it->vertex().z());
    }
  }

protected:
  TH1F *hPt_, *hEta_, *hPhi_, *hNhits_, *hNValidHits_;
  TH1F *hNValidPlusInvalidHits_, *hInnermostHitRadius_, *hInnermostHitZ_;
  TH1F *hDxy_, *hDz_, *hChi2_, *hReferencePointRadius_, *hReferencePointZ_;
  TH1F *hVertexX_, *hVertexY_, *hVertexZ_;
};

#endif // CONTROLPLOTS_H
