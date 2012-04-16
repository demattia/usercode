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
  histo->GetXaxis()->SetTitle("#eta");
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(-1., 1.);
  double binWidth = histo->GetXaxis()->GetBinWidth(1);
  stringstream ss;
  ss << binWidth;
  histo->GetYaxis()->SetTitle(TString("Efficiency/("+ss.str()+")"));
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
  histo->GetXaxis()->SetTitle("#eta");
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(-1.5, 1.5);
  double binWidth = histo->GetXaxis()->GetBinWidth(1);
  stringstream ss;
  ss << binWidth;
  histo->GetYaxis()->SetTitle(TString("Efficiency/("+ss.str()+")"));
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

void ImproveEta()
{
  gStyle->SetOptStat(0);
  TFile * inputFile = new TFile("eff.root", "READ");
  TFile * outputFile = new TFile("improvedEff.root", "RECREATE");

  improveHisto(inputFile, "effVsEta");
  improveHisto(inputFile, "effVsGenEta");

  improveHisto(inputFile, "effVsEta", "effVsGenEta");

  // TFile * inputFileTop = new TFile("eff.root", "READ");
  // TFile * inputFileBottom = new TFile("eff.root", "READ");
  // improveHisto(inputFileTop, "effVsEta", inputFileBottom, "effVsEta");

  outputFile->Write();
  outputFile->Close();
}
