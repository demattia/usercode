#ifndef PLOTS
#define PLOTS

#include <TH1F.h>
#include <TCanvas.h>

class Plots
{
 public:
  Plots(const TString & name)
  {
    hZCandidateMass_ = new TH1F("ZCandidateMass"+name, "ZCandidateMass"+name, 30, 60., 120.);
    hMass_ = new TH1F("TriLeptonMass"+name, "TriLeptonMass"+name, 50, 0, 800.);
    hPtScalarSum_ = new TH1F("TriLeptonPtScalarSum"+name, "TriLeptonPtScalarSum"+name, 80, 0, 800.);
    hMET_ = new TH1F("MET"+name, "MET"+name, 50, 0, 500);
    hVertexDistance_ = new TH1F("VertexDistance"+name, "VertexDistance"+name, 100, 0, 0.02);
  }

  void Fill(const double & zCandidateMass, const double & mass, const double & ptScalarSum, const double & MET, const double & vertexDistance, const double & weight)
  {
    hZCandidateMass_->Fill(zCandidateMass, weight);
    hMass_->Fill(mass, weight);
    hPtScalarSum_->Fill(ptScalarSum, weight);
    hMET_->Fill(MET, weight);
    hVertexDistance_->Fill(vertexDistance, weight);
  }

  void Write()
  {
    TCanvas * canvas = new TCanvas();
    canvas->Divide(2,3);
    canvas->Draw();
    canvas->cd(1);
    hVertexDistance_->Draw();
    canvas->cd(2);
    hMass_->Draw();
    canvas->cd(2);
    hPtScalarSum_->Draw();
    canvas->cd(3);
    hMET_->Draw();
    canvas->cd(4);
    hZCandidateMass_->Draw();

    hZCandidateMass_->Write();
    hMass_->Write();
    hPtScalarSum_->Write();
    hMET_->Write();
    hVertexDistance_->Write();
    canvas->Write();
  }

  TH1F * hZCandidateMass_;
  TH1F * hMass_;
  TH1F * hPtScalarSum_;
  TH1F * hVertexDistance_;
  TH1F * hMET_;
};

#endif
