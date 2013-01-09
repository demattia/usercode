#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <sstream>

void improveHisto(TFile * inputFile, TFile * inputFile2, const TString & histoName)
{
  TH1F * histo = (TH1F*)inputFile->FindObjectAny(histoName);
  TH1F * histo2 = (TH1F*)inputFile2->FindObjectAny(histoName);
  TCanvas * canvas = new TCanvas();
  canvas->Draw();
  histo->Draw();
  histo->GetXaxis()->SetTitle("|d_{0}| [cm]");
  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetRangeUser(0., 50.);
  double binWidth = histo->GetXaxis()->GetBinWidth(1);
  stringstream ss;
  ss << binWidth;
  histo->GetYaxis()->SetTitle(TString("Efficiency / ("+ss.str()+" cm)"));
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetRangeUser(0., 1.1);
  
  histo2->Draw("same");
  histo2->SetMarkerColor(kRed);
  
  TLegend * legend = new TLegend(0.5, 0.67, 0.88, 0.88);
  legend->AddEntry(histo, "MC simulation");
  legend->AddEntry(histo2, "cosmic data");
    // TLegend * legend = canvas->BuildLegend();
  legend->Draw("same");
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  canvas->Write();
}

void Compare()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile * inputFile = new TFile("eff.root", "READ");
  TFile * inputFile2 = new TFile("eff_data.root", "READ");
  TFile * outputFile = new TFile("CompareEff.root", "RECREATE");

  improveHisto(inputFile, inputFile2, "effVsD0");

  outputFile->Write();
  outputFile->Close();
}
