#ifndef CONTROLPLOTS_H
#define CONTROLPLOTS_H

#include <vector>
#include <TH1F.h>

/**
  * Fill a set of control plots for all the tracks in a given collection.
  */
class ControlPlots
{
public:
  ControlPlots( edm::Service<TFileService> & fileService, const TString & name )
  {
    hPt_ = fileService->make<TH1F>(name+"_pt",name+" pt",500,0,500);
    hEta_ = fileService->make<TH1F>(name+"_eta",name+" eta",500,-3,3);
    hPhi_ = fileService->make<TH1F>(name+"_phi",name+" phi",500,-3.2,3.2);
    hNhits_ = fileService->make<TH1F>(name+"_Nhits",name+" Nhits",30,0,30);
    hNValidHits_ = fileService->make<TH1F>(name+"_NValidHits",name+" NValidHits",30,0,30);
    hNValidPlusInvalidHits_ = fileService->make<TH1F>(name+"_NValidPlusInvalidHits",name+" NValidPlusInvalidHits",30,0,30);
    hInnermostHitRadius_ = fileService->make<TH1F>(name+"_innermostHitRadius",name+" innermost hit radius",500,0,100);
    hInnermostHitZ_ = fileService->make<TH1F>(name+"_innermostHitZ",name+" innermost hit z", 500,0,100);
    hDxy_ = fileService->make<TH1F>(name+"_dxy",name+" dxy",500,0,500);
    hDz_ = fileService->make<TH1F>(name+"_dz",name+" dz",500,0,50);
    hChi2_ = fileService->make<TH1F>(name+"_chi2",name+" chi2",500,0,100);
    hReferencePointRadius_ = fileService->make<TH1F>(name+"_RefPointRadius",name+" radius of the reference point",500,-50,50);
    hReferencePointZ_ = fileService->make<TH1F>(name+"_RefPointZ",name+" z of the reference point",500,-50,50);
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
      hDxy_->Fill(it->dxy());
      hDz_->Fill(it->dz());
      hChi2_->Fill(it->normalizedChi2());
      hReferencePointRadius_->Fill(sqrt(it->referencePoint().x()+it->referencePoint().y()));
      hReferencePointZ_->Fill(it->referencePoint().z());
    }
  }

protected:
  TH1F *hPt_, *hEta_, *hPhi_, *hNhits_, *hNValidHits_;
  TH1F *hNValidPlusInvalidHits_, *hInnermostHitRadius_, *hInnermostHitZ_;
  TH1F *hDxy_, *hDz_, *hChi2_, *hReferencePointRadius_, *hReferencePointZ_;
};

#endif // CONTROLPLOTS_H
