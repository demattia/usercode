#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TGraph.h>

#include "tdrstyle.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooVoigtian.h"
#include "RooProdPdf.h"
#include "RooBifurGauss.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"

using namespace RooFit ;

// gSystem->Load("libRooFit");
// gSystem->Load("libRooFitCore");

using namespace RooFit;

#define FIT1M4

void fitBackgroundPDFBWWinter(void)
{
  //input parameters

  // float lowBound=20; // was 50
  // float highBound=596;// was 200

  // float lowBound=78; // was 50
  // float highBound=200;// was 200

  float lowBound=15; // was 50
  float highBound=596;// was 200

  float rebinFactor=2; // was 2
  float scaleFactor=1; // was not used
  bool isLog = true;
  
  // Select here!
  
  bool isMu = true;
  bool isData = false;
  bool isLoose1 = true; // false for loose_2
  bool isPrompt = true; // set for prompt; overrides the above
  
  std::string selection = isLoose1 ? "loose_1" : "loose_2";
  std::string lepTag = isMu ? "mu" : "elec";
  if (isPrompt) {
    lepTag = isMu ? "muon" : "electron";
    selection = "_noLifetimeCuts";
  }
  std::string typeTag = isData ? "data" : "backgroundMC";

  std::string fileNameRoot = "/afs/cern.ch/user/p/plujan/public/LLlimits/WinterSelection/masses_"+typeTag+"_"+lepTag+selection+"_rebin.root";
  std::string fileNameTxt = "/afs/cern.ch/user/p/plujan/public/LLlimits/WinterSelection/masses_"+typeTag+"_"+lepTag+selection+".txt";

  const char *fitTag = (isMu && !isPrompt) ? "BW" : "BWPlus1M4TurnOn";

  // gROOT->ProcessLine(".L ./tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptFit(1);
  
  RooRealVar *rooMass_;
  if (isMu)
    rooMass_ = new RooRealVar("mass","M_{#mu#mu}",lowBound, highBound, "GeV");
  else
    rooMass_ = new RooRealVar("mass","M_{ee}",lowBound, highBound, "GeV");
  RooRealVar Mass = *rooMass_;

  TFile* fout = new TFile("bwFit.root","RECREATE");
  gROOT->cd();

  TFile* file_data =  new TFile(fileNameRoot.c_str(), "READ");
  TH1* data_hist = (TH1D*) file_data->Get("mass");
  data_hist->Sumw2();
  data_hist->Rebin(rebinFactor);
  data_hist->Scale(scaleFactor);
  RooDataHist* data_binned = new RooDataHist("rdh_data","", *rooMass_, data_hist);

  RooRealVar weight("weight", "weight", 0, 50);
  RooRealVar lxySig("lxySig", "lxySig", 0, 100);
  RooDataSet *data_unbinned_tmp = RooDataSet::read(fileNameTxt.c_str(), RooArgList(Mass, lxySig, weight));
  RooDataSet *data_unbinned = new RooDataSet("points", data_unbinned_tmp->GetTitle(), data_unbinned_tmp,
  					     RooArgList(Mass, weight), NULL, weight.GetName());
  // or this version to skip using the weights
  //				             RooArgList(Mass), NULL, 0);


  // Base shape: Breit-Wigner
  RooRealVar mean("BWMean","BWMean", 70.0, 110.0 );
  RooRealVar width("BWWidth","BWWidth", 1.0, 10.0 );
  RooRealVar sigma("VoigtSigma","VoigtSigma", 1.0, 10.0 );
  // RooBreitWigner baseShape("bw","bw",*rooMass_,mean,width);
  RooVoigtian baseShape("bw", "bw", *rooMass_, mean, width, sigma);
  // RooVoigtian signalPdf("bw", "bw", *rooMass_, mean, width, sigma);


  // // Crystal Ball
  // RooRealVar cbmean("CBMean", "CBMean" , 70.0, 110.0) ;
  // RooRealVar cbsigma("CBSigma", "CBSigma" , 1.0, 20.0) ;
  // RooRealVar n("n", "n", 0., 10.);
  RooRealVar alpha("alpha", "alpha", -3., 0.);
  // RooCBShape baseShape("cball", "crystal ball", *rooMass_, cbmean, cbsigma, alpha, n);

  RooRealVar expCoeff("expCoeff", "expCoeff", 0., 10.);


  mean.setConstant(kFALSE);
  width.setConstant(kFALSE);
  sigma.setConstant(kFALSE);

  // 1/m^4 (or other power, check the formula to see what I'm actually doing!) with turnon
  RooRealVar C("C", "C", 0, 1);
  RooRealVar C2("C2", "C2", 0, 1);
  // RooRealVar C3("C3", "C3", 0, 1);
  RooRealVar turnOn("turnOn", "turnOn", 80, 1, 200);
  RooRealVar turnOn2("turnOn2", "turnOn2", 10, 1, 100);
  RooRealVar turnOn3("turnOn3", "turnOn3", 10, 1, 100);
  RooRealVar turnOnWidth("turnOnWidth", "turnOnWidth", 10, 1, 95);
  RooRealVar turnOnWidth2("turnOnWidth2", "turnOnWidth2", 15, 1, 100);
  RooRealVar turnOnWidth3("turnOnWidth3", "turnOnWidth3", 15, 1, 100);
  RooRealVar power("power", "power", 1, 5);
  RooRealVar power2("power2", "power2", 1, 5);
  // RooRealVar power3("power3", "power3", 2, 5);
  RooRealVar powerMinCut("powerMinCut", "powerMinCut", 200, 0, 596);



  // // Convolve with bifurcated Gaussian
  // RooRealVar mean_bg("mean_bg", "mean_bg", 0);
  // RooRealVar sigma_l("sigma_l", "sigma_l", 5.0);
  // RooRealVar sigma_r("sigma_r", "sigma_r", 50.0);
  // mean_bg.setConstant(kFALSE);
  // sigma_l.setConstant(kFALSE);
  // sigma_r.setConstant(kFALSE);
  // RooBifurGauss convShape("bifurGauss", "bifurcated gaussian", Mass, mean_bg, sigma_l, sigma_r);
  // RooFFTConvPdf baseConvShape("sigModel","final signal shape", *rooMass_,  baseShape, convShape);


  // RooRealVar mean2("mean2", "mean2", 70, 150);
  // mean2.setConstant(kFALSE);
  // RooRealVar sigma2("sigma2", "sigma2", 0, 40);
  // RooGaussian G2("G2", "G2", Mass, mean2, sigma2);
  // RooRealVar gaus_rel_frac("C1", "C1", 0, 1);
  // RooAddPdf baseShapePdf("sigModel", "intermediate signal shape", baseShape, G2, gaus_rel_frac);



  // RooGenericPdf powerWithTurnOn("powerWithTurnOn", "(TMath::Erf((@0-@1)/@2)+1)/(pow(@0,@3))", RooArgSet(Mass, turnOn, turnOnWidth, power));
  // RooGenericPdf powerWithDoubleTurnOn("powerWithDoubleTurnOn", "(TMath::Erf((@0-@1)/@2)+1)*(TMath::Erf((@0-@3)/@4)+1)/(pow(@0,@5))", RooArgSet(Mass, turnOn, turnOnWidth, turnOn2, turnOnWidth2, power));
  // RooGenericPdf powerWithTurnOn2("powerWithTurnOn2", "(TMath::Erf((@0-@1)/@2)+1)/(pow(@0,@3))", RooArgSet(Mass, turnOn2, turnOnWidth2, power2));
  // RooGenericPdf powerWithTurnOn3("powerWithTurnOn3", "(TMath::Erf((@0-@1)/@2)+1)/(pow(@0,@3))", RooArgSet(Mass, turnOn3, turnOnWidth3, power3));
  // RooAddPdf signalPdf1("sigModel1", "final signal shape 1", baseShape, powerWithTurnOn, C);


  // RooGenericPdf turnOffPdf("turnOnPdf", "(TMath::Erf((@1-@0)/@2)+1)", RooArgSet(Mass, turnOn2, turnOnWidth2));
  // RooRealVar constant("constant", "constant", 0.5);
  // RooRealVar coefficient("coefficient", "coefficient", 0.);
  // constant.setConstant(kFALSE);
  // coefficient.setConstant(kTRUE);
  // RooGenericPdf linear("linear", "@0*@1+@2", RooArgSet(Mass, constant, coefficient));
  // RooAddPdf signalPdf1("sigModel1", "final signal shape 1", baseShape, linear, C);
  // RooAddPdf signalPdf("sigModel", "final signal shape", signalPdf1, turnOffPdf, C2);


  // // A linear + Voigt*TurnOn
  // RooGenericPdf turnOnPdf("turnOnPdf", "(TMath::Erf((@0-@1)/@2)+1)", RooArgSet(Mass, turnOn2, turnOnWidth2));
  // RooProdPdf voigtWithTurnOn("powerWithTurnOn", "sig", RooArgSet(turnOnPdf, baseShape));
  // // RooRealVar constant("constant", "constant", 0.5, 0, 100);
  // RooRealVar constant("constant", "constant", 60, 0, 100);
  // RooRealVar coefficient("coefficient", "coefficient", 0.);
  // constant.setConstant(kFALSE);
  // coefficient.setConstant(kFALSE);
  // RooGenericPdf linear("linear", "@0*@1+@2", RooArgSet(Mass, coefficient, constant));
  // // RooGenericPdf exponential("expo", "exp((@1-@0)*@2)", RooArgSet(rooMass_, constant, coefficient));
  // RooAddPdf signalPdf1("sigModel1", "final signal shape 1", voigtWithTurnOn, linear, C2);
  // 
  // RooGenericPdf turnOffPdf2("turnOnPdf2", "(TMath::Erf((@1-@0)/@2)+1)", RooArgSet(Mass, turnOn2, turnOnWidth2));
  // RooGenericPdf expoPdf("expo", "exp((@1-@0)*@2)", RooArgSet(Mass, constant, coefficient));
  // RooProdPdf expoTurnOff("expoTurnOff", "sig", RooArgSet(turnOffPdf2, expoPdf));
  // // RooProdPdf signalPdf("expoTurnOff", "sig", RooArgSet(turnOffPdf2, expoPdf));
  // 
  // RooAddPdf signalPdf("sigModel", "final signal shape", signalPdf1, expoTurnOff, C);
  // 
  // // RooGenericPdf signalPdf("linear", "@0*@1+@2", RooArgSet(Mass, coefficient, constant));
  // // RooAddPdf signalPdf("sigModel", "final signal shape", voigtWithTurnOn, linear, C2);


  // ---------------------------- //
  // THIS IS THE CODE FOR THE FIT //
  // ---------------------------- //


  // Linear + turnOff for the rightmost side
  RooRealVar turnOff1("turnOff1", "turnOff1", 80, 20, 80);
  RooRealVar turnOffWidth1("turnOffWidth1", "turnOffWidth1", 5, 1, 15);
  RooGenericPdf turnOffPdf1("turnOffPdf1", "(TMath::Erf((@1-@0)/@2)+1)", RooArgSet(Mass, turnOff1, turnOffWidth1));

  RooRealVar constant("constant", "constant", 60, 0, 100);
  RooRealVar coefficient("coefficient", "coefficient", 0.);
  RooGenericPdf linear("linear", "@0*@1+@2", RooArgSet(Mass, coefficient, constant));

  RooProdPdf leftPdf("leftPdf", "final signal shape 1", RooArgSet(turnOffPdf1, linear));

  // Exponential + turnOff for the radiative tail
  RooRealVar turnOff2("turnOff2", "turnOff2", 80, 1, 200);
  RooRealVar turnOffWidth2("turnOffWidth2", "turnOffWidth2", 10, 1, 95);
  RooGenericPdf turnOffPdf2("turnOffPdf2", "(TMath::Erf((@1-@0)/@2)+1)", RooArgSet(Mass, turnOff2, turnOffWidth2));

  RooRealVar constant1("constant1", "constant1", 60, 0, 100);
  RooRealVar coefficient1("coefficient1", "coefficient1", 0.);
  RooGenericPdf expoPdf("expo", "exp((@1-@0)*@2)", RooArgSet(Mass, constant1, coefficient1));

  RooProdPdf expoTurnOffPdf("expoTurnOff", "sig", RooArgSet(turnOffPdf2, expoPdf));

  // turnOn + Voigtian
  RooGenericPdf turnOnPdf("turnOnPdf", "(TMath::Erf((@0-@1)/@2)+1)", RooArgSet(Mass, turnOn, turnOnWidth));
  RooProdPdf voigtWithTurnOn("powerWithTurnOn", "sig", RooArgSet(turnOnPdf, baseShape));

  // Combine all
  RooAddPdf signalPdf1("sigModel1", "final signal shape 1", leftPdf, expoTurnOffPdf, C);
  RooAddPdf signalPdf("sigModel", "final signal shape", signalPdf1, voigtWithTurnOn, C2);





  // RooAddPdf signalPdf("sigModel", "final signal shape", baseShape, turnOffPdf, C);


  // RooRealVar constant("constant", "constant", 0.5);
  // RooRealVar coefficient("coefficient", "coefficient", 0.5);
  // constant.setConstant(kFALSE);
  // coefficient.setConstant(kFALSE);
  // RooGenericPdf linear("linear", "@0*@1+@2", RooArgSet(Mass, constant, coefficient));

  // RooAddPdf signalPdf1("sigModel1", "final signal shape 1", baseShapePdf, powerWithTurnOn, C);
  // RooAddPdf signalPdf1("sigModel1", "final signal shape 1", baseConvShape, powerWithTurnOn, C);
  // RooAddPdf signalPdf("sigModel2", "final signal shape 2", signalPdf1, powerWithTurnOn2, C2);
  // RooAddPdf signalPdf("sigModel2", "final signal shape 2", signalPdf1, turnOnPdf, C2);
  // RooAddPdf signalPdf("sigModel2", "final signal shape 2", signalPdf1, linear, C2);
  // RooAddPdf signalPdf("sigModel", "final signal shape", signalPdf2, powerWithTurnOn3, C3);

  // RooAddPdf signalPdf("sigModel", "final signal shape", baseShape, powerWithTurnOn, C);
  // RooAddPdf signalPdf("sigModel", "final signal shape", baseShape, powerWithDoubleTurnOn, C);

  // RooGenericPdf turnOnPdf("turnOnPdf", "(TMath::Erf((@0-@1)/@2)+1)", RooArgSet(Mass, turnOn2, turnOnWidth2));
  // RooProdPdf voigtWithTurnOn("powerWithTurnOn", "sig", RooArgSet(turnOnPdf, baseShape));
  // RooProdPdf signalPdf("powerWithTurnOn", "sig", RooArgSet(turnOnPdf, baseShape));
  // RooAddPdf signalPdf("sigModel", "final signal shape", voigtWithTurnOn, powerWithTurnOn, C);


  // RooGenericPdf massPower("massPower", "1/(pow(@0,@1))", RooArgSet(Mass, power));
  // RooGenericPdf massPower("massPower", "(@0>@2)/(pow(@0,@1))", RooArgSet(Mass, power, powerMinCut));
  // RooAddPdf signalPdf("sigModel", "final signal shape", baseShape, massPower, C);

  // Multiple power laws
  // RooGenericPdf massPower2("massPower2", "1/(pow(@0,@1))", RooArgSet(Mass, power2));
  // RooGenericPdf massPower3("massPower3", "1/(pow(@0,@1))", RooArgSet(Mass, power3));
  // RooAddPdf signalMass1("sigModel1", "final signal shape 1", baseShape, massPower, C);
  // RooAddPdf signalMass2("sigModel2", "final signal shape 2", signalMass1, massPower2, C2);
  // RooAddPdf signalPdf("sigModel", "final signal shape", signalMass2, massPower3, C3);


  // RooGenericPdf signalPdf("massPower", "1/(pow(@0,@1))", RooArgSet(Mass, power));
  // RooGenericPdf signalPdf("massPower", "(@0>@2)/(pow(@0,@1))", RooArgSet(Mass, power, powerMinCut));



  // if (!isMu && isData && !isPrompt) {
  //   turnOnWidth.setVal(15.0);
  // }
  // if (isPrompt) {
  //   turnOn.setVal(90.0);
  //   turnOn.setConstant(kTRUE);
  //   //turnOnWidth.setVal(5.0);
  // }



  // Do the fit to the data. You can use either data_binned or data_binned here.

  // RooFitResult *fitResult = signalPdf.fitTo(*data_unbinned, Save(true), 
  RooFitResult *fitResult = signalPdf.fitTo(*data_binned, Save(true), 
					    Extended(false), 
					    Minos(false),
					    Optimize(2),
					    PrintEvalErrors(-1),
					    Warnings(false), SumW2Error(true)
					    );
  fitResult->Print("v");

  // tdrStyle->UseCurrentStyle();
  setTDRStyle();
  // tdrStyle->cd();

  // Create binning object with range (10,500)
  RooBinning tbins(15,596);
  // tbins.addUniform(100,20,60);
  // tbins.addUniform(100,60,150);
  // tbins.addUniform(100,150,500);
  tbins.addUniform(581,15,596);


  // Try using Andrzej's binning
  // double boundaries[] =  {20.0, 29.6, 39.2, 48.8, 58.4, 68.0, 77.6, 87.2, 96.8, 106.4, 117.818, 128.904, 141.205, 155.715, 174.239, 200.061, 239.145, 302.336, 409.346, 596.0};
  // RooBinning tbins(19, boundaries);


  TString cname = Form("fit");
  TCanvas* c1 = new TCanvas(cname,cname,500,500);
  RooPlot* frame1 = Mass.frame(lowBound, highBound);
  data_unbinned->plotOn(frame1,DataError(RooAbsData::SumW2),Binning(tbins),Rescale(1.0/tbins.averageBinWidth()));
  signalPdf.plotOn(frame1,ProjWData(*data_unbinned),Normalization(1.0/tbins.averageBinWidth(), RooAbsPdf::Relative));
  signalPdf.paramOn(frame1,Layout(0.5, 0.95, 0.92));
  std::cout<<"chi2/dof = " << frame1->chiSquare(4) << std::endl;
  frame1->Draw("e0");
  frame1->GetYaxis()->SetTitle("Events / GeV");
  if (isLog) {
    if (isPrompt)
      frame1->SetMinimum(1*scaleFactor);
    else {
      if (isMu) {
	frame1->SetMinimum(0.05*scaleFactor);
      } else {
	if (isData)
	  frame1->SetMinimum(0.3*scaleFactor);
	else
	frame1->SetMinimum(0.01*scaleFactor);
      }
    }
    gPad->SetLogy();
  }
  frame1->SetMaximum(frame1->GetMaximum()/tbins.averageBinWidth());
  frame1->SetMinimum(frame1->GetMinimum()/tbins.averageBinWidth());
  if (!isMu)
    frame1->SetTitle("CMS #sqrt{s}=7 TeV L=4.1 fb^{-1}");
  else
    frame1->SetTitle("CMS #sqrt{s}=7 TeV L=5.1 fb^{-1}");
  
  fout->cd();
  c1->Write();
  fout->Write();
  char buf[512];
  sprintf(buf, "BackgroundPDFWinter_%s%s%s_%s%s_New.png",
	  isMu ? "muon" : "electron", isData ? "Data" : "MC", isPrompt ? "Prompt" : isLoose1 ? "L1" : "L2", fitTag);
  c1->Print(buf);
  sprintf(buf, "BackgroundPDFWinter_%s%s%s_%s%s_New.pdf",
	  isMu ? "muon" : "electron", isData ? "Data" : "MC", isPrompt ? "Prompt" : isLoose1 ? "L1" : "L2", fitTag);
  c1->Print(buf);

}
