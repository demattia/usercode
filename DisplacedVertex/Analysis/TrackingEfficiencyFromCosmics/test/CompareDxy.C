#include "TFile.h"
#include "TCanvas.h"

#include "tdrstyle.C"

void CompareDxy(const TString & file1, const TString & file2)
{
  setTDRStyle();

  TFile inputFileData(file1, "READ");
  TFile inputFileSim(file2, "READ");

  TDirectory * dirData = (TDirectory*)inputFileData.Get("demo");
  TH1F * histoData = (TH1F*)dirData->Get("standAloneMuons_dxy");

  TDirectory * dirSim = (TDirectory*)inputFileSim.Get("demo");
  TH1F * histoSim = (TH1F*)dirSim->Get("standAloneMuons_dxy");

  double norm = histoData->Integral(1, histoData->FindBin(40.));
  if( norm ) histoData->Scale(1/norm);

  norm = histoSim->Integral(1, histoSim->FindBin(40.));
  if( norm ) histoSim->Scale(1/norm);

  histoSim->SetLineColor(kRed);

  TCanvas * canvas = new TCanvas();

  histoData->Draw();
  histoSim->Draw("Same");

  canvas->SaveAs("compareDxy.png");
  canvas->SaveAs("compareDxy.pdf");
}
