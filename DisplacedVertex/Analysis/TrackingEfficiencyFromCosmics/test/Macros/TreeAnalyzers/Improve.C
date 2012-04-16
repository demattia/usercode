#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <sstream>

void improveHisto(TFile * inputFile, const TString & histoName, const TString & histoName2 = "")
{
  TH1F * histo = (TH1F*)inputFile->FindObjectAny(histoName);
  TCanvas * canvas = new TCanvas();
  canvas->Draw();
  histo->Draw();
  histo->GetXaxis()->SetTitle("|d_{0}| [cm]");
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(0., 50.);
  double binWidth = histo->GetXaxis()->GetBinWidth(1);
  stringstream ss;
  ss << binWidth;
  histo->GetYaxis()->SetTitle(TString("Efficiency/("+ss.str()+" cm)"));
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetRangeUser(0., 1.1);
  if( histoName2 != "" ) {
    TH1F * histo2 = (TH1F*)inputFile->FindObjectAny(histoName2);
    histo2->Draw("same");
    histo2->SetMarkerColor(kRed);
    TLegend * legend = new TLegend(0.5, 0.67, 0.88, 0.88);
    legend->AddEntry(histo, "cosmicMuons1Leg");
    legend->AddEntry(histo2, "MC truth");
    // TLegend * legend = canvas->BuildLegend();
    legend->Draw("same");
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(0);
  }
  canvas->Write();
}

void improveHisto(TFile * inputFile, const TString & histoName, TFile * inputFile2, const TString & histoName2)
{
  TH1F * histo = (TH1F*)inputFile->FindObjectAny(histoName);
  TCanvas * canvas = new TCanvas();
  canvas->Draw();
  histo->Draw();
  histo->GetXaxis()->SetTitle("|d_{0}| [cm]");
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(0., 50.);
  double binWidth = histo->GetXaxis()->GetBinWidth(1);
  stringstream ss;
  ss << binWidth;
  histo->GetYaxis()->SetTitle(TString("Efficiency/("+ss.str()+" cm)"));
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetRangeUser(0., 1.1);
  TH1F * histo2 = (TH1F*)inputFile2->FindObjectAny(histoName2);
  histo2->Draw("same");
  histo2->SetMarkerColor(kRed);
  TLegend * legend = new TLegend(0.5, 0.67, 0.88, 0.88);
  legend->AddEntry(histo, "top");
  legend->AddEntry(histo2, "bottom");
  legend->Draw("same");
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  canvas->Write();
}

void Improve()
{
  gStyle->SetOptStat(0);
  TFile * inputFile = new TFile("eff.root", "READ");
  TFile * outputFile = new TFile("improvedEff.root", "RECREATE");

  improveHisto(inputFile, "effVsD0");
  improveHisto(inputFile, "effVsGenD0");

  improveHisto(inputFile, "effVsD0", "effVsGenD0");

  TFile * inputFileTop = new TFile("eff_50bins_top.root", "READ");
  TFile * inputFileBottom = new TFile("eff_50bins_bottom.root", "READ");
  improveHisto(inputFileTop, "effVsD0", inputFileBottom, "effVsD0");

  outputFile->Write();
  outputFile->Close();
}
