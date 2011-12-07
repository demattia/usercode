#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <sstream>
#include "TROOT.h"

#include "Compare.C"

TH1F * getPlot(const TString & name)
{
  TFile * file = new TFile(name+".root", "READ");
  TCanvas * canvas = (TCanvas*)file->Get(name);
  return (TH1F*)canvas->GetListOfPrimitives()->FindObject(name);
}

void setLogY(TCanvas * canvas, const int index, const bool set = true)
{
  std::stringstream ss;
  ss << index;
  TPad *pad = (TPad*)canvas->GetListOfPrimitives()->FindObject(TString("canvas_")+ss.str());
  pad->cd();
  pad->SetLogy(set);
}

void drawPlot(const TString & name, TCanvas * canvas, const int index, const bool logY)
{
  TH1F * histo = getPlot(name);
  setLogY(canvas, index, logY);
  histo->Draw();
}

void PutPlots()
{
  gROOT->SetBatch(true);
  Compare();
  gROOT->SetBatch(false);

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1000, 800);
  canvas->Divide(2,3);

  drawPlot("deltaR", canvas, 1, true);
  drawPlot("deltaPt", canvas, 2, true);
  drawPlot("deltaNhits", canvas, 3, true);
  drawPlot("deltaChi2NDOF", canvas, 4, true);
  drawPlot("unmatchedMu1Chi2NDOF", canvas, 5, true);
  drawPlot("unmatchedMu2Chi2NDOF", canvas, 6, true);

  canvas->SaveAs("comparison.root");
  canvas->SaveAs("comparison.png");
}
