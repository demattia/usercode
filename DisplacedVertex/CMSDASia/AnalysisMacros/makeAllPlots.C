#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TH1F.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TROOT.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>

TString getTypeName(const bool electrons, const bool trackAnalysis, const bool useSAMuons)
{
  TString type("");
  if ( trackAnalysis ) {
    if( electrons ) type = "_electrons";
    else type = "_muons";
  }
  else {
    if ( useSAMuons ) type = "_saMuons";
    else  type = "_globalMuons";
  }

  return type;
}

float getLumi(const bool electrons)
{
  TString lumiFile = "";
  if( electrons ) lumiFile = "LumiFiles/lumi_electron.txt";
  else lumiFile = "LumiFiles/lumi_muon.txt";

  // Read file to get lumi (in /pb)
  ifstream openFile;
  openFile.open(lumiFile.Data());
  // Only one line in file...
  string lumiString = "0";
  if (openFile.is_open()) {
    getline(openFile, lumiString);
    openFile.close();
  }
  return ::atof(lumiString.c_str());
}

//TLegend *legend= new TLegend(0.7,0.6,0.9,0.85);
void addHistogram(const TString & fileName, const TString & histoName, const int color, THStack & hs, const TString & leg, TLegend *legend, const float lumi, double & integral, double & integralError2 )
{
  TFile * inputFile= new TFile (fileName, "READ");
  TH1F * histo = (TH1F*)inputFile->Get(histoName);
  TH1::SetDefaultSumw2();
  histo->Scale(lumi);

  // Sum total integral of background MC
  double integralError=0;
  integral += histo->IntegralAndError(1,histo->GetNbinsX(),integralError);
  integralError2 += integralError * integralError;

  histo->SetFillColor(color);
  //  histo->Rebin(5);
  hs.Add(histo);

  legend->AddEntry(histo, leg, "f");
}

void makePlot(const TString & histoName, TFile * outputFile, TString type, const TString & xTitle,  const double & yMin = 0., const double & yMax = 0., const bool electrons = false, const float lumi = 0)
{
  std::cout << "----> Making plots for : " << histoName.Data() << std::endl;
  TH1::SetDefaultSumw2();

  // Make title for stacks
  TString stackTitle = "CMS #sqrt{s} = 8 TeV, L = ";
  char s[32];
  std::cout << "lumi = " << lumi/1000 << std::endl;
  sprintf(s, "%.1f",lumi/1000);
  stackTitle += s;
  stackTitle +="fb^{-1}";

  THStack hs("hs",stackTitle);
  TString latexText = "#mu^{+}#mu^{-}";
  if (electrons) {
    hs.SetTitle(stackTitle);
    // hs = THStack("hs",);
    latexText = "e^{+}e^{-}";
  }
  //  cout <<histoName <<endl;
//  TString type(getTypeName(electrons, trackAnalysis, useSAMuons));
  std::cout << "Getting data file" << std::endl;
  TFile * inputFile= new TFile ("CombinedFiles/Data_combined"+type+".root", "READ");
  std::cout << "Getting data histogram" << std::endl;
  TH1F * histoData = (TH1F*)inputFile->Get(histoName);

  std::cout << "Calculating integral of data histo" << std::endl;
  double dataIntegralError = 0;
  double dataIntegral = histoData->IntegralAndError( 1, histoData->GetNbinsX(), dataIntegralError );
  //histoData->Scale(150);
  // histoData->SetMarkerStyle(20);
  // histoData->SetMarkerSize(1);
  // histoData->SetMarkerStyle(21);
  // histoData->SetMarkerSize(0.5);
  //histoData->Rebin(5);

  double mcIntegral = 0;
  double mcIntegralError2 = 0;
  std::cout << "Doing MC" << std::endl;
  TLegend *legend= new TLegend(0.65,0.6,0.85,0.85);
  addHistogram("CombinedFiles/ZZ_combined"+type+".root", histoName, 8, hs,"ZZ", legend, lumi, mcIntegral, mcIntegralError2);
  addHistogram("CombinedFiles/WZ_combined"+type+".root", histoName, 7, hs,"WZ", legend, lumi, mcIntegral, mcIntegralError2);
  addHistogram("CombinedFiles/WW_combined"+type+".root", histoName, 6, hs,"WW", legend, lumi, mcIntegral, mcIntegralError2);
  // addHistogram("CombinedFiles/WJets_combined"+type+".root", histoName, 5, hs,"WJets", legend, lumi, mcIntegral, mcIntegralError2);
  //  addHistogram("Ztautau_combined"+type+".root", histoName, 97, hs,"Z/#gamma*->#tau#tau",legend, , lumi);
  //  if( electrons ) {
  //    addHistogram("Zee_combined"+type+".root", histoName, kRed, hs,"Z/#gamma*->ee", legend, , lumi);
  //  }
  //  else {
  //    addHistogram("Zmumu_combined"+type+".root", histoName, kRed, hs,"Z/#gamma*->#mu#mu", legend, lumi);
  //  }
  addHistogram("CombinedFiles/DYJets_combined"+type+".root", histoName, kRed, hs,"DYJets", legend, lumi, mcIntegral, mcIntegralError2);
  // addHistogram("TTbar_combined"+type+".root", histoName, 3, hs,"tt", legend, lumi);
  addHistogram("CombinedFiles/TTJets_combined"+type+".root", histoName, 3, hs,"tt+jets", legend, lumi, mcIntegral, mcIntegralError2);
  addHistogram("CombinedFiles/QCD_combined"+type+".root", histoName, kBlue, hs,"QCD", legend, lumi, mcIntegral, mcIntegralError2);

  
  std::cout << "Integral of Data : " << dataIntegral << " +/- " << dataIntegralError << endl;
  std::cout << "Integral of MC : " << mcIntegral << " +/- " << sqrt(mcIntegralError2) << endl;

  outputFile->cd();
  TCanvas canvas(histoName+"_canvas");
  canvas.cd();
  gPad->SetLogy();
  canvas.Draw();
  canvas.SetFillColor(kWhite);
  canvas.SetBorderMode(0);
  hs.Draw("HISTE");
  hs.GetXaxis()->SetLimits(0, 500);
  double xlatex = 20;
  if(xTitle=="mass[Gev/c^{2}]") {
    hs.GetXaxis()->SetLimits(0, 500);
    xlatex = 20;
  }
  else if (xTitle=="Sum pt, low pt lepton" || xTitle=="Sum pt, high pt lepton") {
    hs.GetXaxis()->SetLimits(0, 25);
    xlatex = 20;
  }
  else if (xTitle=="N reco PV" || xTitle=="N vtx true") {
    gPad->SetLogy(0);
    hs.GetXaxis()->SetLimits(0, 60);
    xlatex = 20;
  }
  else if (xTitle=="cos{#alpha}") {
    gPad->SetLogy(1);
    hs.GetXaxis()->SetLimits(-1.1, 1.1);
    xlatex = 20;
  }
  else if (xTitle=="Vertex #chi^{2}") {
    hs.GetXaxis()->SetLimits(0, 20);
    xlatex = 20;
  }
  else if (xTitle=="L_{xy}") {
    hs.GetXaxis()->SetLimits(-20, 20);
    xlatex = 20;
    legend->SetX1(0.2);
    legend->SetX2(0.35);
  }
  else if (xTitle=="pt, low pt lepton" || xTitle=="pt, high pt lepton") {
    hs.GetXaxis()->SetLimits(0, 200);
    xlatex = 20;
  }
  else if (xTitle=="d0 H" ||xTitle=="d0 L" ) {
    hs.GetXaxis()->SetLimits(0, 0.5);
    xlatex = 20;
  }
  else {
    hs.GetXaxis()->SetLimits(-40, 40);
    xlatex = 1;
    legend->SetX1(0.2);
    legend->SetX2(0.35);
  }
  if( yMin != 0. && yMax != 0. ) {
    hs.SetMinimum(yMin);
    hs.SetMaximum(yMax);
  }
  hs.GetXaxis()->SetTitle(xTitle);
  hs.GetYaxis()->SetTitle("Entries");
  histoData->SetMarkerStyle(20);
  histoData->SetMarkerSize(1);
  histoData->Draw("same,P,E,X0");

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histoData,"data","p");
  legend->Draw();

  TLatex *l1 = new TLatex(xlatex, yMax*7/10, latexText);
  l1->Draw();

  canvas.Write();

  outputFile->Write();
  //  outputFile->Close();
}

void makeAllPlots(const bool electrons = false, const bool trackAnalysis = true, const bool useSAMuons = false)
{
  //const int electon = 0;
  TString type(getTypeName(electrons, trackAnalysis, useSAMuons));

  // Get lumi
  float lumi = getLumi(electrons);

  std::cout << "Total lumi in data : " << lumi << std::endl;

  TFile * outputFile = new TFile("output"+type+".root", "RECREATE");
  // makePlot("caloCorrMass_nolifetime_inverted", outputFile, type, "mass[Gev/c^{2}]", "e^{+}e^{-}", 1, 0.07, 1500000);
  makePlot("Mass_nolifetime_inverted", outputFile, type, "mass[Gev/c^{2}]", 0.07, 10000000, electrons, lumi);
  makePlot("signedDecayLengthSignificance_lifetime_nodecaylength",outputFile, type,"L_{xy}/#sigma", 0.001, 1000, electrons, lumi);
  makePlot("Mass_lifetime_notinverted",outputFile, type,"mass[Gev/c^{2}]", 0.001, 1000, electrons, lumi);
  makePlot("signedDecayLengthSignificance_nolifetime_nodecaylength",outputFile, type,"L_{xy}/#sigma", 0.01, 10000000, electrons, lumi);
  makePlot("signedDecayLength_nolifetime_nodecaylength",outputFile, type,"L_{xy}", 0.01, 10000000, electrons, lumi);
  makePlot("cosine_nolifetime_inverted",outputFile, type,"cos{#alpha}",0,10000, electrons, lumi );
  makePlot("vertexChi2_nolifetime_inverted",outputFile, type,"Vertex #chi^{2}",1000,10000000,electrons, lumi);
  makePlot("nPV_nolifetime_nodecaylength",outputFile, type,"N reco PV",0,10000, electrons, lumi );
  makePlot("nvtx_true_nolifetime_nodecayLength",outputFile, type,"N vtx true", 0, 10000, electrons, lumi);
  makePlot("isolationPtH_nolifetime_inverted",outputFile, type,"Sum pt, high pt lepton",1.0,10000000, electrons, lumi );
  makePlot("isolationPtL_nolifetime_inverted",outputFile, type,"Sum pt, low pt lepton",1.0,10000000, electrons, lumi );
  makePlot("ptH_nolifetime_inverted",outputFile, type,"pt, high pt lepton",1.0,10000000, electrons, lumi );
  makePlot("ptL_nolifetime_inverted",outputFile, type,"pt, low pt lepton",1.0,10000000, electrons, lumi );
//  makePlot("D0H_nolifetime_nodecaylength",outputFile, type,"d0 H",1.0,10000000, electrons, lumi );
//  makePlot("D0L_nolifetime_nodecaylength",outputFile, type,"d0 L",1.0,10000000, electrons, lumi );
  outputFile->Close();
}

