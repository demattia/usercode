#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <TH1F.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TROOT.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>

TString getTypeName(const bool electrons)
{
  if( electrons ) return "_electrons";
  return "_muons";
}

void addHistogram(const TString & fileName, const TString & histoName, const int color, THStack & hs, const TString & leg, TLegend *legend, const bool electrons)
{
  TFile * inputFile= new TFile (fileName, "READ");
  TH1F * histo = (TH1F*)inputFile->Get(histoName);
  TH1::SetDefaultSumw2();
  if ( electrons ) {
    histo->Scale(5016.8312988830003);
  }
  else {
    histo->Scale(3398.51040391+698.80137458);
    // histo->Scale(698.80137458);
    // histo->Scale(5000);
  }
  histo->SetFillColor(color);
  hs.Add(histo);

  legend->AddEntry(histo, leg, "f");
}

void makePlot(const TString & histoName, TFile * outputFile, const TString & xTitle,  const TString & yTitle, const double & yMin = 0., const double & yMax = 0., const bool electrons = false)
{
  TH1::SetDefaultSumw2();

  THStack hs("hs","CMS #sqrt{s} = 7 TeV, L = 1.3fb^{-1}");
  TString latexText = "#mu^{+}#mu^{-}";
  if (electrons) {
    hs.SetTitle("CMS #sqrt{s} = 7 TeV, L = 3.9fb^{-1}");
    latexText = "e^{+}e^{-}";
  }
  //  cout <<histoName <<endl;
  TString type(getTypeName(electrons));
  TFile * inputFile= new TFile ("Data_combined"+type+".root", "READ");
  TH1F * histoData = (TH1F*)inputFile->Get(histoName);
  TFile * inputSignalFile= new TFile ("Signal_combined"+type+".root", "READ");
  TH1F * histoSignal = (TH1F*)inputSignalFile->Get(histoName);
  TLegend *legend= new TLegend(0.65,0.6,0.85,0.85);
  addHistogram("ZZ_combined"+type+".root", histoName, 8, hs,"ZZ", legend, electrons);
  addHistogram("WW_combined"+type+".root", histoName, 7, hs,"WW", legend, electrons);
  addHistogram("WZ_combined"+type+".root", histoName, 6, hs,"WZ", legend, electrons);
  addHistogram("DYJets_combined"+type+".root", histoName, kRed, hs,"DYJets", legend, electrons);
  addHistogram("TTJets_combined"+type+".root", histoName, 3, hs,"tt+jets", legend, electrons);
  addHistogram("QCD_combined"+type+".root", histoName, kBlue, hs,"QCD", legend, electrons);
  
  outputFile->cd();
  TCanvas canvas("canvas_"+histoName);
  canvas.cd();
  gPad->SetLogy();
  canvas.Draw();
  canvas.SetFillColor(kWhite);
  canvas.SetBorderMode(0);
  hs.Draw("HISTE");
  if( yMin != 0. && yMax != 0. ) {
    hs.SetMinimum(yMin);
    hs.SetMaximum(yMax);
  }
  hs.GetXaxis()->SetTitle(xTitle);
  hs.GetYaxis()->SetTitle(yTitle);
  histoData->SetMarkerStyle(20);
  histoData->SetMarkerSize(1);
  histoData->Draw("same,P,E,X0");

  // histoSignal->SetMarkerStyle(20);
  histoSignal->SetLineColor(kRed);
  histoSignal->Draw("sameL");

  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(histoData,"data","p");
  legend->Draw();

  canvas.Write();

  outputFile->Write();
}

void makeAllPlots(const bool electrons = false)
{
  TString type(getTypeName(electrons));
  TFile * outputFile = new TFile("output"+type+".root", "RECREATE");

  makePlot("hZLeptonPtL", outputFile, "p_{T} [GeV/c]", "Entries / 1 GeV", 0, 0, electrons);
  makePlot("hZLeptonPtH", outputFile, "p_{T} [GeV/c]", "Entries / 1 GeV", 0, 0, electrons);
  makePlot("hThirdLeptonPt", outputFile, "p_{T} [GeV/c]", "Entries / 1 GeV", 0, 0, electrons);

  makePlot("ZCandidateMass_beforeAllCuts", outputFile, "mass [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("VertexDistance_beforeAllCuts", outputFile, "[cm]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonMass_beforeAllCuts", outputFile, "mass [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonPtScalarSum_beforeAllCuts", outputFile, "mass [GeV/c]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("MET_beforeAllCuts", outputFile, "ME_{T} [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);

  makePlot("ZCandidateMass_afterTriggerAndMETCut", outputFile, "mass [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("VertexDistance_afterTriggerAndMETCut", outputFile, "[cm]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonMass_afterTriggerAndMETCut", outputFile, "mass[Gev/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonPtScalarSum_afterTriggerAndMETCut", outputFile, "mass[Gev/c]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("MET_afterTriggerAndMETCut", outputFile, "ME_{T} [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);

  makePlot("ZCandidateMass_afterAllCuts", outputFile, "mass [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("VertexDistance_afterAllCuts", outputFile, "[cm]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonMass_afterAllCuts", outputFile, "mass[Gev/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("TriLeptonPtScalarSum_afterAllCuts", outputFile, "mass[Gev/c]", "mu^{+}mu^{-}", 0, 0, electrons);
  makePlot("MET_afterAllCuts", outputFile, "ME_{T} [GeV/c^{2}]", "mu^{+}mu^{-}", 0, 0, electrons);

  outputFile->Close();
}

