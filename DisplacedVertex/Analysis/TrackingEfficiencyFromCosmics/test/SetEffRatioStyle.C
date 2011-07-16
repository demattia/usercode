#include "tdrstyle.C"
#include <sstream>

void drawStyled(const TString & fileDir, const TString & name, const TString & xTitle, const TString & type)
{
  TFile * file = new TFile(fileDir+name+".root");
  TCanvas * canvas = (TCanvas*)file->Get("c1"+name);
  TGraphErrors * graph = (TGraphErrors*)canvas->GetPrimitive(name);

  TString underscoreType("");
  if( type != "" ) underscoreType = "_" + type;

  TCanvas * newCanvas = new TCanvas();
  newCanvas->Draw();
  graph->Draw("AP");
  graph->GetXaxis()->SetTitle(xTitle);
  stringstream ss;
  ss << 2*graph->GetErrorX(0);
  graph->GetYaxis()->SetTitle("Eff(data)/Eff(sim) / ( "+TString(ss.str())+" cm )");
  graph->GetYaxis()->SetRangeUser(0, 3.);
  graph->GetXaxis()->SetRangeUser(0, 40.);
  graph->GetYaxis()->SetTitleOffset(1.2);
  newCanvas->SaveAs("styled_"+name+underscoreType+".root");
  newCanvas->SaveAs("styled_"+name+underscoreType+".pdf");
  newCanvas->SaveAs("styled_"+name+underscoreType+".png");
}

void SetEffRatioStyle(const TString & fileDir = "", const TString & type = "")
{
  setTDRStyle();

  drawStyled(fileDir, "EffVsDxyRatio", "|dxy| [cm]", type);
}
