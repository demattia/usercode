#include "tdrstyle.C"
#include <sstream>

void drawStyled(const TString & fileDir, const TString & name, const TString & xTitle, const TString & type, const double xMax = 0, const TString & unit = "")
{
  TFile * file = new TFile(fileDir+name+".root");
  TCanvas * canvas = (TCanvas*)file->Get("c1"+name);
  TGraphErrors * graph = (TGraphErrors*)canvas->GetPrimitive(name);

  TString underscoreType("");
  if( type != "" ) underscoreType = "_" + type;
  TString spacedUnit = "";
  if( unit != "" ) spacedUnit = " "+unit;


  TCanvas * newCanvas = new TCanvas();
  newCanvas->Draw();
  graph->Draw("AP");
  stringstream ss;
  ss << 2*graph->GetErrorX(0);
  graph->GetXaxis()->SetTitle(xTitle);
  graph->GetYaxis()->SetTitle("Efficiency / ( "+TString(ss.str())+spacedUnit+" )" );
  if( xMax != 0 ) graph->GetXaxis()->SetRangeUser(0, xMax);
  graph->GetYaxis()->SetRangeUser(0, 1.1);
  graph->GetYaxis()->SetTitleOffset(1.2);
  newCanvas->SaveAs("styled_"+name+underscoreType+".root");
  newCanvas->SaveAs("styled_"+name+underscoreType+".pdf");
  newCanvas->SaveAs("styled_"+name+underscoreType+".png");
}

void SetStyle(const TString & fileDir = "", const TString & type = "")
{
  setTDRStyle();

  drawStyled(fileDir, "EffVsDxy", "|dxy| [cm]", type, 40., "cm");
  drawStyled(fileDir, "EffVsDz", "|dz| [cm]", type);
  drawStyled(fileDir, "EffVsPt", "Pt [GeV/c]", type);
}
