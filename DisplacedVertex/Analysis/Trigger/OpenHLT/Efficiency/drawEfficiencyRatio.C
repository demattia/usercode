#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"

void fillEff(const TString & fileName, TH2F *& effSignal, TString histoName, TFile * outputFile, int & numSignalSamples)
{

  TFile * inputFile = new TFile(fileName, "READ");
  TH2F * effMap = (TH2F*)inputFile->Get("effMap");
  int entries = ((TH1F*)inputFile->Get("entriesHisto"))->GetEntries();

  Int_t nBinsX = effMap->GetNbinsX();
  Int_t nBinsY = effMap->GetNbinsY();

  if( effSignal == 0 ) {
    outputFile->cd();
    effSignal = new TH2F(histoName, histoName,
                         nBinsX, effMap->GetXaxis()->GetXmin(), effMap->GetXaxis()->GetXmax(),
                         nBinsY, effMap->GetYaxis()->GetXmin(), effMap->GetYaxis()->GetXmax());
  }

  for( Int_t i=1; i<=nBinsX; ++i ) {
    for( Int_t j=1; j<=nBinsY; ++j ) {
      effSignal->SetBinContent(i, j, effSignal->GetBinContent(i, j) + double(effMap->GetBinContent(i, j))/entries);
    }
  }
  ++numSignalSamples;
}

void normalize(TH2F * effHisto, const int numSamples)
{
  for( Int_t i=1; i<=effHisto->GetNbinsX(); ++i ) {
    for( Int_t j=1; j<=effHisto->GetNbinsY(); ++j ) {
      effHisto->SetBinContent(i, j, double(effHisto->GetBinContent(i, j))/numSamples);
    }
  }
}

void drawEfficiencyRatio()
{
  TFile * outputFile = new TFile("effRatio.root", "RECREATE");

  TH2F * effSignal = 0;
  int numSignalSamples = 0;
  fillEff("effMap_MH120MFF50.root", effSignal, "effSignal", outputFile, numSignalSamples);
  fillEff("effMap_MH1000MFF20.root", effSignal, "effSignal", outputFile, numSignalSamples);
  fillEff("effMap_MH400MFF50.root", effSignal, "effSignal", outputFile, numSignalSamples);

  normalize(effSignal, numSignalSamples);

  TH2F * effBackground = 0;
  int numBackroundSamples = 0;
  fillEff("effMap_part2.root", effBackground, "effBackground", outputFile, numBackroundSamples);

  if( effSignal != 0 && effBackground != 0 ) {
    Int_t nBinsX = effSignal->GetNbinsX();
    Int_t nBinsY = effSignal->GetNbinsY();

    TH2F * effRatio = new TH2F("effRatio", "errRatio",
                               nBinsX, effSignal->GetXaxis()->GetXmin(), effSignal->GetXaxis()->GetXmax(),
                               nBinsY, effSignal->GetYaxis()->GetXmin(), effSignal->GetYaxis()->GetXmax());

    for( Int_t i=1; i<=nBinsX; ++i ) {
      for( Int_t j=1; j<=nBinsY; ++j ) {
        effRatio->SetBinContent(i, j, double(effSignal->GetBinContent(i, j))/sqrt(double(effBackground->GetBinContent(i, j))));
      }
    }

    effSignal->Write();
    effBackground->Write();
    effRatio->SetXTitle("pt [GeV/c]");
    effRatio->GetXaxis()->SetTitleOffset(1.4);
    effRatio->SetYTitle("parallel cut");
    effRatio->GetYaxis()->SetTitleOffset(1.4);
    effRatio->Write();
  }


  outputFile->Write();
}
