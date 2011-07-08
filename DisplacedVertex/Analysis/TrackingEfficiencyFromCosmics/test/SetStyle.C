#include "tdrstyle.C"

void drawStyled(const TString & name, const TString xTitle)
{
  TFile * file = new TFile(name+".root");
  TCanvas * canvas = (TCanvas*)file->Get("c1"+name);
  TGraphErrors * graph = (TGraphErrors*)canvas->GetPrimitive(name);

  TCanvas * newCanvas = new TCanvas();
  newCanvas->Draw();
  graph->Draw("AP");
  graph->GetXaxis()->SetTitle(xTitle);
  graph->GetYaxis()->SetTitle("Efficiency");
  graph->GetYaxis()->SetRangeUser(0, 1.1);
  graph->GetYaxis()->SetTitleOffset(1.5);
  newCanvas->SaveAs("styled_"+name+".root");
}

void SetStyle()
{
  setTDRStyle();

  drawStyled("EffVsDxy", "|dxy| [cm]");
  drawStyled("EffVsDz", "|dz| [cm]");
  drawStyled("EffVsPt", "Pt [GeV/c]");
}
