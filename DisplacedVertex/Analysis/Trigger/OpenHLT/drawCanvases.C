#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TROOT.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>

TString noCutsName("_NoCuts_");
TString oneValidHitName("_OneValidHit_");
TString oneValidChamberName("_OneValidChamber_");
TString parallelismCutName("_ParallelismCut_");

TH1 * getHisto(const TString & histoName, const TString & dirName, TFile * inputFile)
{
  TDirectory * dir = (TDirectory*)inputFile->Get(dirName);
  return (TH1*)dir->Get(histoName);
}

void drawCanvas(const TString & firstPart, const TString & secondPart, TFile * inputFile)
{
  TH2F * noCutsHisto = (TH2F*)getHisto(firstPart+noCutsName+secondPart, noCutsName, inputFile);
  TH2F * oneValidHitHisto = (TH2F*)getHisto(firstPart+oneValidHitName+secondPart, oneValidHitName, inputFile);
  TH2F * oneValidChamberHisto = (TH2F*)getHisto(firstPart+oneValidChamberName+secondPart, oneValidChamberName, inputFile);
  TH2F * parallelismCutHisto = (TH2F*)getHisto(firstPart+parallelismCutName+secondPart, parallelismCutName, inputFile);

  TCanvas * canvas = new TCanvas(firstPart+secondPart, firstPart+secondPart, 2000, 400);
  canvas->Divide(4,1);
  canvas->cd(1);
  noCutsHisto->Draw();
  noCutsHisto->SetMarkerStyle(0);
  canvas->cd(2);
  oneValidHitHisto->Draw();
  oneValidHitHisto->SetMarkerStyle(0);
  canvas->cd(3);
  oneValidChamberHisto->Draw();
  oneValidChamberHisto->SetMarkerStyle(0);
  canvas->cd(4);
  parallelismCutHisto->Draw();
  parallelismCutHisto->SetMarkerStyle(0);
  // parallelismCutHisto->SetMarkerStyle(20);
  // parallelismCutHisto->SetMarkerStyle(0.8);
  canvas->SaveAs(firstPart+secondPart+".pdf");
  canvas->SaveAs(firstPart+secondPart+".png");
}

void drawSuperimposed(const TString & firstPart, const TString & secondPart, TFile * inputFile)
{
  TH2F * noCutsHisto = (TH2F*)getHisto(firstPart+noCutsName+secondPart, noCutsName, inputFile);
  TH2F * oneValidHitHisto = (TH2F*)getHisto(firstPart+oneValidHitName+secondPart, oneValidHitName, inputFile);
  TH2F * oneValidChamberHisto = (TH2F*)getHisto(firstPart+oneValidChamberName+secondPart, oneValidChamberName, inputFile);
  TH2F * parallelismCutHisto = (TH2F*)getHisto(firstPart+parallelismCutName+secondPart, parallelismCutName, inputFile);
  TCanvas * canvas = new TCanvas(firstPart+secondPart);
  canvas->cd();
  noCutsHisto->Draw();
  noCutsHisto->GetYaxis()->SetRangeUser(std::min(noCutsHisto->GetMaximum(), 0.), noCutsHisto->GetMaximum()*1.3);
  oneValidHitHisto->Draw("same");
  oneValidHitHisto->SetLineColor(kRed);
  oneValidChamberHisto->Draw("same");
  oneValidChamberHisto->SetLineColor(kBlue);
  parallelismCutHisto->Draw("same");
  parallelismCutHisto->SetLineColor(8);

  TLegend *leg = new TLegend(0.4966443,0.7937063,1,1,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  std::stringstream ss;
  ss << noCutsHisto->GetEntries();
  leg->SetHeader(firstPart+secondPart+" with "+TString(ss.str())+" entries");
  leg->AddEntry(noCutsHisto,"No cuts","l");
  // leg->AddEntry(oneValidHitHisto,"At least one valid hit","l");
  leg->AddEntry(oneValidHitHisto,"Standard trigger cuts","l");
  leg->AddEntry(oneValidChamberHisto,"More than 1 valid chamber","l");
  leg->AddEntry(parallelismCutHisto,"Parallelism cut","l");
  leg->Draw();

  canvas->SaveAs(firstPart+secondPart+".pdf");
  canvas->SaveAs(firstPart+secondPart+".png");
}

void countEntries(const TString & firstPart, const TString & secondPart, TFile * inputFile)
{
  Int_t noCutsEntries = ((TH2F*)getHisto(firstPart+noCutsName+secondPart, noCutsName, inputFile))->GetEntries();
  Int_t oneValidHitEntries = ((TH2F*)getHisto(firstPart+oneValidHitName+secondPart, oneValidHitName, inputFile))->GetEntries();
  Int_t oneValidChamberEntries = ((TH2F*)getHisto(firstPart+oneValidChamberName+secondPart, oneValidChamberName, inputFile))->GetEntries();
  Int_t parallelismCutEntries = ((TH2F*)getHisto(firstPart+parallelismCutName+secondPart, parallelismCutName, inputFile))->GetEntries();

  std::cout << "Total number of entries = " << noCutsEntries << std::endl;
  std::cout << "Entries after default trigger cuts = " << oneValidHitEntries << ", ratio with previous cut = "
            << oneValidHitEntries/double(noCutsEntries) << std::endl;
  std::cout << "Entries after > 1 valid chamber = " << oneValidChamberEntries << ", ratio with previous cut = "
            << oneValidChamberEntries/double(oneValidHitEntries) << ", total eff = " << oneValidChamberEntries/double(noCutsEntries) << std::endl;
  std::cout << "Entries after parallelism cut = " << parallelismCutEntries << ", ratio with previous cut = "
            << parallelismCutEntries/double(oneValidChamberEntries) << ", total eff with respect to default trigger cuts = "
            << parallelismCutEntries/double(oneValidHitEntries) << ", total eff = " << parallelismCutEntries/double(noCutsEntries) << std::endl;

  ofstream countEntriesFile;
  countEntriesFile.open("entriesFile.txt");
  countEntriesFile << "Total number of entries = " << noCutsEntries << "<br/>"  << std::endl;
  countEntriesFile << "Entries after default trigger cuts = " << oneValidHitEntries << ", ratio with previous cut = "
                   << oneValidHitEntries/double(noCutsEntries) << "<br/>" << std::endl;
  countEntriesFile << "Entries after > 1 valid chamber = " << oneValidChamberEntries << ", ratio with previous cut = "
                   << oneValidChamberEntries/double(oneValidHitEntries) << ", total eff = "
                   << oneValidChamberEntries/double(noCutsEntries) << "<br/>"  << std::endl;
  countEntriesFile << "Entries after parallelism cut = " << parallelismCutEntries << ", ratio with previous cut = "
                   << parallelismCutEntries/double(oneValidChamberEntries) << ", total eff with respect to default trigger cuts = "
                   << parallelismCutEntries/double(oneValidHitEntries) << ", total eff = "
                   << parallelismCutEntries/double(noCutsEntries) << "<br/>"  << std::endl;
  countEntriesFile.close();
}

void drawCanvases()
{
  gStyle->SetOptStat(0);

  TFile * inputFile = new TFile("CheckOpenHLT.root", "READ");
  drawSuperimposed("pt", "First", inputFile);
  drawSuperimposed("eta", "First", inputFile);
  drawSuperimposed("pt", "Second", inputFile);
  drawSuperimposed("eta", "Second", inputFile);
  drawSuperimposed("nhits", "First", inputFile);
  drawSuperimposed("nhits", "Second", inputFile);
  drawSuperimposed("nchambers", "First", inputFile);
  drawSuperimposed("nchambers", "Second", inputFile);
  drawSuperimposed("Parallelism", "", inputFile);
  drawCanvas("eta", "Correlation_0_1", inputFile);
  drawCanvas("pt", "Correlation_0_1", inputFile);
  drawCanvas("phi", "Correlation_0_1", inputFile);
  countEntries("eta", "Correlation_0_1", inputFile);
}
