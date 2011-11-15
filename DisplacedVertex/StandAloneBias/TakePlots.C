#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TROOT.h"

void printHisto(const TDirectory * dir, const TString & name)
{
  TH1F * histo = (TH1F*)dir->Get(name);
  std::cout << "name = " << name << std::endl;
  TCanvas canvas("canvas", "canvas", 500, 500);
  histo->Draw();
  canvas.Print(name+".png");
}

void TakePlots()
{
  gROOT->SetBatch(true);

  TFile * inputFile = new TFile("TrackingEfficiencyFromCosmics.root", "READ");
  TDirectory * dir = (TDirectory *)inputFile->Get("demo");

  printHisto(dir, "genTracks_pt");
  printHisto(dir, "genTracks_eta");
  printHisto(dir, "genTracks_phi");
  printHisto(dir, "genTracks_d_0");
  printHisto(dir, "genTracks_d_z");

  printHisto(dir, "generalTracks_pt");
  printHisto(dir, "generalTracks_eta");
  printHisto(dir, "generalTracks_phi");
  printHisto(dir, "generalTracks_d_0");
  printHisto(dir, "generalTracks_d_z");
  printHisto(dir, "generalTracks_chi2");
  printHisto(dir, "trackVsGenDelta_PtPulls");
  printHisto(dir, "trackVsGenDelta_DxyPulls");
  printHisto(dir, "trackVsGenDelta_DzPulls");

  printHisto(dir, "standAloneMuons_pt");
  printHisto(dir, "standAloneMuons_eta");
  printHisto(dir, "standAloneMuons_phi");
  printHisto(dir, "standAloneMuons_d_0");
  printHisto(dir, "standAloneMuons_d_z");
  printHisto(dir, "standAloneMuons_Nhits");
  printHisto(dir, "standAloneMuons_NValidHits");
  printHisto(dir, "standAloneMuons_chi2");
  printHisto(dir, "standAloneVsGenDelta_PtPulls");
  printHisto(dir, "standAloneVsGenDelta_DxyPulls");
  printHisto(dir, "standAloneVsGenDelta_DzPulls");

  printHisto(dir, "cleanedStandAloneMuons_pt");
  printHisto(dir, "cleanedStandAloneMuons_eta");
  printHisto(dir, "cleanedStandAloneMuons_phi");
  printHisto(dir, "cleanedStandAloneMuons_d_0");
  printHisto(dir, "cleanedStandAloneMuons_d_z");
  printHisto(dir, "cleanedStandAloneMuons_Nhits");
  printHisto(dir, "cleanedStandAloneMuons_NValidHits");
  printHisto(dir, "cleanedStandAloneMuons_chi2");
  printHisto(dir, "cleanedStandAloneVsGenDelta_PtPulls");
  printHisto(dir, "cleanedStandAloneVsGenDelta_DxyPulls");
  printHisto(dir, "cleanedStandAloneVsGenDelta_DzPulls");

}
