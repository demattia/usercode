//example of fit where the model is histogram + function
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h" 
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TChain.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>

// function to minimize
void chi2function(int &npar, double *grad, double &fval, double *par, int flag);
void likelihoodfunction(int &npar, double *grad, double &fval, double *par, int flag);
// function to fit data
// 2 components: sgn and bkg
// 1 parameter: sgn fraction
double func(double sgn, double bkg, double *par);

// construct a smeared template from the original one
// smear within xsigma each bin
TH1D* smearedHistoConstruction( TH1D *histo, double xsigma );

void pseudoDataGeneration(TH1D* bkgHisto, TH1D* sgnHisto);

void fitHistoGeneration(TH1D* bkgHisto, TH1D* sgnHisto, double par_val, double par_err);

void drawCanvas(TCanvas* c, TH1D* fithisto, TH1D* data, TH1D* bkg, TH1D* sgn, double k2, double fittedBkgvalue, double fittedBkgerror, double fittedSgnvalue, double fittedSgnerror, int nf, int ipseudoexp);

void computePull(double expectedParvalue, double constrainedParerror, double fittedParvalue,double fittedParerror, TH1D* significance, TH1D* pull);

// global variables
// histograms
TH1D *backgroundHisto;
TH1D *signalHisto;
TH1D *pseudoData;
TH1D *totHisto; // used only for plotting fit result
                // totHisto[i] = par*paseudoData[i]+(1-par)*pseudoData[i]
                // what about its error???
static int   nbin_;
static double xmin_;
static double xmax_;

// number of entries for the different histograms
double bkgEntries = 0; // input background template
double sgnEntries = 0; // input signal template
double pseudoDataEntries        = 0; // generated pseudo data
double bkg_in_pseudoDataEntries = 0; // background used in pseudo data
double sgn_in_pseudoDataEntries = 0; // signal used in pseudo data

// number of bin content for the different histograms
double bkg_content  = 0.; 
double sgn_content  = 0.; 
double data_content = 0.;

// number of bin error
double data_error = 0.00000001;  

bool   flag_bkgSmear  = kFALSE;
double sigma_bkgSmear = 1.0;
bool   flag_sgnSmear  = kFALSE;
double sigma_sgnSmear = 1.0;
double bkgScale = 1.0;
double sgnScale = 1.0;

// variables used by fitting procedure
int    ndof      =  0;
double chi       =  0.; 
double like      = -9.;
double uncert_bkgEntries = 0.01; // uncertainty on background knowledge
                                 // used as parameter in chi2 constraint
                                 // be carefull in pull computation!!!
                                 // ..hope to have done in the right way

// parameters used by minimization procedure 
const int ilun = 5;
const int olun = 6;
const int elun = 7;
int flag_fitMethod = 0;  // flag to specify minimation method
  	                 // 0: chi2 [default]
	                 // 1: likelihood

// number of fit parameters
const int nparam = 1;

// fit function ftotal to signal + background
void fit(TString  inputFileName,                  // input file name 
	 TString  bkHisto,                        // input background template name
	 TString  sgHisto,                        // input signal template name
	 TString  outputFileName,                 // root name for output files
	 int      fitMethod             = 0,      // flag to specify minimation method
	                                          // 0: chi2 [default]
	                                          // 1: likelihood
	 double   signalfrac            = 0.01,   // parameter starting value [used by minuit]
	 double   stepPar               = 0.001,  // parameter step in minization [used by minuit]
	 double   minPar                = 0.0,    // parameter minimum value [used by minuit]
	 double   maxPar                = 1.0,    // parameter maximum value [used by minuit]
	 double   bkgEntriesUncertainty = 0.01,   // uncertainty on background knowledge
                                                  // used as parameter in chi2 constraint
	 bool     bkgSmear              = kFALSE, // flag to specify smearing on background template
	                                          // kFALSE: off [default]
	                                          // kTRUE: on
	 double   bkgSigmaSmear         = 1.0,    // sigma in smearing background template 
	 bool     sgnSmear              = kFALSE, // flag to specify smearing on signal template
	                                          // kFALSE: off [default]
	                                          // kTRUE: on
	 double   sgnSigmaSmear         = 1.0,    // sigma in smearing signal template 
	 double   bkgScaleFactor        = 1.0,    // factor to scale input background template
	 double   sgnScaleFactor        = 1.0,    // factor to scale input signal template
	 /*	 
		double inputLum              = 100., // fb^-1
		double minLum                = 0.1,
		double maxLum                = 10.,
		double stepLum               = 0.1,
	 */
	 int      nPseudoExp            = 1000,   // number of pseudo experiment
	 int      migrad_call           = 500     // number of minuit call 
	                                          // - does it work??
	                                          // - is it necessary??
	 ) { 
  gROOT->Reset(); 
  // avoid canvas pop up
  gROOT->SetBatch(kTRUE); // kTRUE: turn off pop up
                          // kFALSE: turn on pop up
  // general root setting
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  
  // set input parameters
  flag_fitMethod = fitMethod; // 0: chi2
                              // 1: likelihood
  uncert_bkgEntries = bkgEntriesUncertainty;
  flag_bkgSmear  = bkgSmear;
  sigma_bkgSmear = bkgSigmaSmear;
  flag_sgnSmear  = sgnSmear;
  sigma_sgnSmear = sgnSigmaSmear;
  bkgScale = bkgScaleFactor;
  sgnScale = sgnScaleFactor;

  // open input file and get bkg and sgn histograms
  TFile *inputFile = new TFile(inputFileName);
  backgroundHisto = (TH1D*)inputFile->Get(bkHisto);  //global pointer used in ftotal
  signalHisto     = (TH1D*)inputFile->Get(sgHisto);  //global pointer used in ftotal

  // get input histograms bin parameters
  // background and signal templlate shoud have the same parameters
  nbin_ = backgroundHisto->GetNbinsX();
  xmin_ = backgroundHisto->GetXaxis()->GetXmin();
  xmax_ = backgroundHisto->GetXaxis()->GetXmax();

  //get input background and signal entries
  // it is necessary to use this loop 
  // since sometimes the input histograms could have GetEntries()==nbin_!!
  double bkgEntries_inSgnBins = 0;
  double significanceEstimated = 0.; 
  for ( int ibin=1; ibin <= nbin_; ++ibin ) {
    bkgEntries += backgroundHisto->GetBinContent(ibin);  
    sgnEntries += signalHisto->GetBinContent(ibin);         
    if ( signalHisto->GetBinContent(ibin) != 0. ) 
      bkgEntries_inSgnBins+=backgroundHisto->GetBinContent(ibin);  
  }
  // std::cout << "bkgEntries: " << bkgEntries << std::endl;
  // std::cout << "sgnEntries: " << sgnEntries << std::endl;
  
  /*
  // luminosity calculation:
  // in integrated luminosity L = inputLum (100 fb^-1) 
  // MC background correspond to ... events
  // MC signal     correspond to ... events
  double factorLum = iLum / ( double()/ * inputLum );
  // Fixed signal fraction in each kin zone
  //---------------------------------------
  totentries = (float) (nhist_ref2_mcbg->Integral() * lumfactor); // with lum factor
  */
  
  int binEstimated = 4*int(TMath::Sqrt(nPseudoExp));
  double sgn = sgnEntries*sgnScale;
  double bkg = bkgEntries*bkgScale;
  double sgnfracvalueEstimated = (sgn+10.*TMath::Sqrt(sgn))/
    (sgn-10.*TMath::Sqrt(sgn)+bkg-10.*TMath::Sqrt(bkg));
  double sgnevtsvalueEstimated = sgn+10.*TMath::Sqrt(sgn);
  double sgnfracerrorEstimated = sgnfracvalueEstimated/2;
  double sgnevtserrorEstimated = sgnevtsvalueEstimated/2;
  sgnfracvalueEstimated *= 3;
  sgnevtsvalueEstimated *= 3;
  if(bkgEntries_inSgnBins!=0.) significanceEstimated = sgnEntries/TMath::Sqrt(bkgEntries_inSgnBins);
  significanceEstimated *= 3;

  // create output file
  TFile *outputFile = new TFile(outputFileName+".root","recreate");
  
  // define pull histograms
  TH1D * minChi2                  = new TH1D("minChi2",                 "minimum #chi^{2}", binEstimated,0.,2.);
  TH2D * fittedSgnFracVSminChi2   = new TH2D("fittedSgnFracVSminChi2",  "fitted signal fraction VS minimum #chi^{2}",binEstimated,0.,2.,binEstimated,0.,sgnfracvalueEstimated);
  TH2D * fittedSgnEventsVSminChi2 = new TH2D("fittedSgnEventsVSminChi2","fitted signal number of events VS minimum #chi^{2}",binEstimated,0.,2.,binEstimated,0.,sgnevtsvalueEstimated);
  TH1D * fittedSgnFrac_value   = new TH1D("fittedSgnFrac_value",  "fitted signal fraction",                  binEstimated,0.,sgnfracvalueEstimated);
  TH1D * fittedSgnFrac_error   = new TH1D("fittedSgnFrac_error",  "error on fitted signal fraction",         2*binEstimated,0.,sgnfracerrorEstimated);
  TH1D * fittedSgnEvents_value = new TH1D("fittedSgnEvents_value","fitted signal number of events",          binEstimated,0.,sgnevtsvalueEstimated);
  TH1D * fittedSgnEvents_error = new TH1D("fittedSgnEvents_error","error on fitted signal number of events", 2*binEstimated,0.,sgnevtserrorEstimated);
  TH1D * fittedSgnFrac_significance   = new TH1D("fittedSgnFrac_significance","significance on fitted signal fraction",2*binEstimated,  0., significanceEstimated);
  TH1D * fittedSgnFrac_pull           = new TH1D("fittedSgnFrac_pull",        "pull on fitted signal fraction",        binEstimated,-5., 5.);
  
  // create pseudo experiment canvas
  TCanvas *canvas[nPseudoExp];    
  stringstream canvasName;
  for ( int i = 0; i < nPseudoExp; ++i ) {
    canvasName << "fitCanvasPseudoExp_" << i;
    canvas[i] = new TCanvas(TString(canvasName.str()),"Fitting Demo",900,800);
    canvasName.str("");
    canvas[i]->UseCurrentStyle();
    canvas[i]->SetGrid();
  }
  // create signal fraction canvas
  TCanvas *canvas_fittedSgnFrac = new TCanvas("sgnFrac_CanvasPseudoExp","pull on signal fraction Fitting Demo",900,800);    
  canvas_fittedSgnFrac->Divide(2,2);
  canvas_fittedSgnFrac->SetGrid();
  TCanvas *canvas_fittedSgn = new TCanvas("sgn_CanvasPseudoExp","pull on signal number of events Fitting Demo",900,800);
  canvas_fittedSgn->Divide(2,2);
  canvas_fittedSgn->SetGrid();
  // create pull canvas
  TCanvas *canvas_pull = new TCanvas("pullCanvas","pull on signal fraction Fitting Demo",900,800);
  canvas_pull->Divide(2,2);
  canvas_pull->SetGrid();
  
  // figure output
  // type = 111   portrait  ps
  // type = 112   landscape ps
  // type = 113   eps
  TPostScript *pseudoexpPostscript = new TPostScript(outputFileName+".ps",112);

  // minimization algorithm
  // fit params
  double fitparams[nPseudoExp];
  double fiterr[nPseudoExp];
  double chi2[nPseudoExp];
  int ndf[nPseudoExp];
  for ( int i = 0; i < nPseudoExp; ++i ) {
    fitparams[i] = 0.;
    fiterr[i]    = 0.;
    chi2[i]      = 0.;
    ndf[i]       = 0;
  }
  
  pseudoData = new TH1D("pseudoData","pseudoData",nbin_,xmin_,xmax_);
  totHisto   = new TH1D("totHisto","total fit templates",nbin_,xmin_,xmax_);

  // loop on pseudo experiments
  for ( int i = 0; i < nPseudoExp; ++i ) {
    std::cout << "*** " << i << " pseudo-exp ***" << std::endl;
    
    // clean pseudo-data and fit result histograms
    pseudoData->Reset();
    totHisto->Reset();  
    
    // generation of pseudo data 
    // from input background and signal template
    pseudoDataGeneration(backgroundHisto,signalHisto);
    
    ndof = 0;
    TMinuit *gMinuit = new TMinuit(1);
    // set fit function
    if ( flag_fitMethod == 0 ) gMinuit->SetFCN(chi2function);
    if ( flag_fitMethod == 1 ) gMinuit->SetFCN(likelihoodfunction);
    
    int ipar = 0;
    double arglist[1];
    int ierflg = 0;
    gMinuit->mninit( ilun, olun, elun );
    // set starting values, step size and range for parameters
    double par_start[nparam] = {signalfrac};
    double par_step[nparam]  = {stepPar};
    double par_min[nparam]   = {minPar};
    double par_max[nparam]   = {maxPar};
    gMinuit->mnparm(ipar, "signal fraction", 
		    par_start[ipar], par_step[ipar], 
		    par_min[ipar],par_max[ipar],
		    ierflg);
    arglist[0] = 2;
    // ???
    gMinuit->mnexcm("set str",arglist,1,ierflg);
    

    // do minimisation
    chi = 0;
    arglist[0] = migrad_call;
    gMinuit->mnexcm( "call fcn", arglist, 1, ierflg );
    gMinuit->mnexcm( "migrad",   arglist, 0, ierflg );
    gMinuit->mnexcm( "minos",    arglist, 0, ierflg );
    //gMinuit->mnexcm( "stop",     arglis, 0, ierflg );
    
    // some output
    double parval;
    TString name;
    double pmin;
    double pmax;
    double errp;
    double errl;
    double errh;
    int ivar;
    double erro;
    double cglo;
    gMinuit->mnpout( ipar, name, parval, erro, pmin, pmax, ivar );
    gMinuit->mnerrs( ipar, errp, errl, errh, cglo );
    double fedm; 
    double errdef; 
    int npari; 
    int nparx; 
    int istat;
    gMinuit->mnstat(chi, fedm, errdef, npari, nparx, istat);

    std::cout << "-----------------------------------------------------------------"<< std::endl;
    std::cout << "| " << name << " (parameter " << ipar << "): " 
	              << parval << " -+ " << errp << std::endl;
    std::cout << "| w/ minos error neg= " << errl << " minos error pos= " << errh << std::endl;
    std::cout << "| erro: " << erro << " pmin: " << pmin << " pmax: " << pmax
	      << " ivar: "  << ivar << " cglo: " << cglo << std::endl;
    fitparams[i]=parval;
    // set error on params
    if ( errp!=0. ) fiterr[i]=errp;
    else  fiterr[i]= (((errh)>(fabs(errl)))?(errh):(fabs(errl)));    
    chi2[i] = chi;
    ndf[i] = ndof - nparam;
    double chi_grad = chi2[i] / double( ndf[i] ); 
    std::cout << "| ndf[i]: "     << ndf[i]   << std::endl;
    std::cout << "| chi2 / ndf: " << chi_grad << std::endl;
    std::cout << "-----------------------------------------------------------------"<< std::endl;
    //    std::cout << "fedm: " << fedm << std::endl;
    //    std::cout << "errdef: " << errdef << std::endl;
    //    std::cout << "npari: " << npari << std::endl;
    //    std::cout << "nparx: " << nparx << std::endl;
    //    std::cout << "istat: " << istat << std::endl;
    
    /*
    // output latex file
    ofstream latexfile;
    char val[50];
    latexfile.open(outputFileName+".tex", ofstream::out | ofstream::app);
    latexfile << "\\hline" << endl;
    sprintf(val,"%.2f \\pm %.2f",(*val_sgnFrac)(nrow,0),(*err_sgnFrac)(nrow,0));
    latexfile << "signal fraction: & $" << val << "$\\\\" << endl;
    sprintf(val,"%.2f / %.0f",(*ki2)(nrow,0),(*NDF)(nrow,0));
    latexfile << "$\\chi^2 / ndf$: & $" << val << "$\\\\" << endl;
    sprintf(val,"%.0f",(*nTot)(nrow,0));
    latexfile << "$N_{evt}$: & $" << val << "$\\\\" << endl;     
    latexfile << "\\hline" << endl;
    */

    double chi2_grad       = chi2[i]/ndf[i];
    double fittedSgn_value = fitparams[i];
    double fittedSgn_error = fiterr[i];
    double fittedSgnEntries_value = fitparams[i]*pseudoDataEntries;
    double fittedSgnEntries_error = fittedSgn_error*pseudoDataEntries;
    double fittedBkgEntries_value = (1.-fittedSgn_value)*pseudoDataEntries;
    double fittedBkgEntries_error =      fittedSgn_error*pseudoDataEntries;
    // Display version
    // plot fitted distributions
    fitHistoGeneration(backgroundHisto,signalHisto,fittedSgn_value,fittedSgn_error);
    backgroundHisto->SetNormFactor(fittedBkgEntries_value);
    signalHisto->SetNormFactor(fittedSgnEntries_value);

    drawCanvas(canvas[i],totHisto,pseudoData,backgroundHisto,signalHisto,
	       fittedBkgEntries_value,fittedBkgEntries_error,
	       fittedSgnEntries_value,fittedSgnEntries_error,
	       chi2[i],ndf[i],i);
    if ( i < 10 ) pseudoexpPostscript->NewPage();
    if ( i < 10 ) canvas[i]->Update();

    // plot pull distributions
    minChi2->Fill(chi2_grad);
    fittedSgnFracVSminChi2->Fill(chi2_grad,fittedSgn_value);
    fittedSgnFrac_value->Fill(fittedSgn_value);
    fittedSgnFrac_error->Fill(fittedSgn_error);

    fittedSgnEventsVSminChi2->Fill(chi2_grad,fittedSgnEntries_value);
    fittedSgnEvents_value->Fill(fittedSgnEntries_value);
    fittedSgnEvents_error->Fill(fittedSgnEntries_error);
    double expectedSgnFrac_value = sgn_in_pseudoDataEntries/
                                   pseudoDataEntries;
    double constrainedSgnFrac_error = uncert_bkgEntries;
    computePull(expectedSgnFrac_value,constrainedSgnFrac_error,fittedSgn_value,fittedSgn_error,fittedSgnFrac_significance,fittedSgnFrac_pull);
  }

  pseudoexpPostscript->Close();
  //  // invoke Postscript viewer
  //   gSystem->Exec("gs file.ps");

  double fitSgnMean = fittedSgnEvents_value->GetMean();
  double fitSgnErrr = fittedSgnEvents_value->GetMeanError();
  double minunc = 3*fitSgnErrr;
  int pseudoexp = 0;
  for ( int i = 0; i < nPseudoExp; ++i ) {
    double fitSgn = fitparams[i]*pseudoDataEntries;
    double unc = TMath::Abs(fitSgn-fitSgnMean);
    if (unc<minunc) {
      minunc=unc;
      pseudoexp=i;
      std::cout << "minunc: "    << minunc << std::endl;
      std::cout << "pseudoexp: " << pseudoexp << std::endl;
    }
  }
  gStyle->SetOptStat(0);
  canvas[pseudoexp]->Print(outputFileName+"_pseudoExp.eps");


  double x    = 0.;
  double ymin = 0.;
  double ymax = 0.;
  x = sgn_in_pseudoDataEntries/
    pseudoDataEntries;
  ymax = fittedSgnFrac_value->GetMaximum();
  TLine * expSgnFrac = new TLine(x,ymin,x,ymax);
  expSgnFrac->SetLineColor(kGreen);
  expSgnFrac->SetLineWidth(2);
  x = sgn_in_pseudoDataEntries;
  ymax = fittedSgnEvents_value->GetMaximum();
  TLine * expSgnEvents = new TLine(x,ymin,x,ymax);
  expSgnEvents->SetLineColor(kGreen);
  expSgnEvents->SetLineWidth(2);

  //  pseudoexpPostscript->NewPage();
  gStyle->SetOptStat(111110);
  canvas_fittedSgnFrac->UseCurrentStyle();
  canvas_fittedSgnFrac->cd(1);
  minChi2->Draw();
  canvas_fittedSgnFrac->cd(2);
  fittedSgnFracVSminChi2->Draw();
  canvas_fittedSgnFrac->cd(3);
  fittedSgnFrac_value->Draw();
  expSgnFrac->Draw("same");
  canvas_fittedSgnFrac->cd(4);
  fittedSgnFrac_error->Draw();
  canvas_fittedSgnFrac->Update();


  //  pseudoexpPostscript->NewPage();
  canvas_fittedSgn->UseCurrentStyle();
  canvas_fittedSgn->cd(1);
  minChi2->Draw();
  canvas_fittedSgn->cd(2);
  fittedSgnEventsVSminChi2->Draw();
  canvas_fittedSgn->cd(3);
  fittedSgnEvents_value->Draw();
  expSgnEvents->Draw("same");
  canvas_fittedSgn->cd(4);
  fittedSgnEvents_error->Draw();
  canvas_fittedSgn->Update();

  //  pseudoexpPostscript->NewPage();
  gStyle->SetOptFit(1111111);
  canvas_pull->UseCurrentStyle();
  canvas_pull->cd(1);
  fittedSgnFrac_value->Draw();
  expSgnFrac->Draw("same");
  canvas_pull->cd(3);
  fittedSgnEvents_value->Draw();
  expSgnEvents->Draw("same");
  canvas_pull->cd(2);
  fittedSgnFrac_significance->Draw();
  canvas_pull->cd(4);
  fittedSgnFrac_pull->Fit("gaus");
  fittedSgnFrac_pull->Draw();
  canvas_pull->Update();
  canvas_pull->Print(outputFileName+"_pull.eps");

  
  TF1 * gausPull = fittedSgnFrac_pull->GetFunction("gaus");
  double mean_value  = gausPull->GetParameter(1);
  double mean_error  = gausPull->GetParError(1);
  double sigma_value = gausPull->GetParameter(2);
  double sigma_error = gausPull->GetParError(2);
  double significance_mean  = fittedSgnFrac_significance->GetMean();
  double significance_error = fittedSgnFrac_significance->GetMeanError();

  // output latex file
  ofstream latexfile;
  char val[50];
  string variable(bkHisto);
  variable.substr(0,variable.find("_"));
  latexfile.open(outputFileName+".tex", ofstream::out | ofstream::app);
  latexfile << "\\hline" << endl;
  sprintf(val,"%s",variable.c_str());
  latexfile << "variable: & " << val << "$\\\\" << endl;
  latexfile << "\\hline" << endl;
  sprintf(val,"%.3f \\pm %.3f",mean_value,mean_error);
  latexfile << "pull mean: & $" << val << "$\\\\" << endl;
  sprintf(val,"%.3f / %.3f",sigma_value,sigma_error);
  latexfile << "pull sigma: & $" << val << "$\\\\" << endl;
  latexfile << "\\hline" << endl;
  sprintf(val,"%.3f / %.3f",significance_mean,significance_error);
  latexfile << "significance: & $" << val << "$\\\\" << endl;
  latexfile << "\\hline" << endl;

  latexfile.close();

  std::cout << "sgnEntries*sgnScale= " << sgnEntries << "*" << sgnScale << "= sgn: " << sgn << std::endl;
  std::cout << "bkgEntries*bkgScale= " << bkgEntries << "*" << bkgScale << "= bkg: " << bkg << std::endl;
  std::cout << "sgnfracvalueEstimated: " << sgnfracvalueEstimated 
	    << "sgnfracerrorEstimated: " << sgnfracerrorEstimated 
	    << " sgnevtsvalueEstimated: " << sgnevtsvalueEstimated
	    << " sgnevtserrorEstimated: " << sgnevtserrorEstimated 
	    << " significanceEstimated: " << significanceEstimated << std::endl;

  outputFile->cd();
  backgroundHisto->Write();
  signalHisto->Write();
  minChi2->Write();
  fittedSgnFracVSminChi2->Write();
  fittedSgnEventsVSminChi2->Write();
  fittedSgnFrac_value->Write();
  fittedSgnFrac_error->Write();
  fittedSgnEvents_value->Write();
  fittedSgnEvents_error->Write();
  fittedSgnFrac_significance->Write();
  fittedSgnFrac_pull->Write();

  // end of macro
  return;

}

// construct a smeared template from the original one
// smear within xsigma each bin
TH1D* smearedHistoConstruction( TH1D *histo, double xsigma ) {

  TString histoName = histo->GetName();
  // create smeared histo
  TH1D* smearedHisto = new TH1D(histoName+"_smeared","Smeared "+histoName,nbin_,xmin_,xmax_);
  smearedHisto->Sumw2();
  for ( int ibin=1; ibin<=nbin_; ++ibin ){
    double content = histo->GetBinContent(ibin);
    double error   = TMath::Sqrt(content);
    // fluctuate bin
    double newcontent = gRandom->Gaus(content,xsigma*error); 
    // content within x sigma
    smearedHisto->SetBinContent(ibin,newcontent);
    smearedHisto->SetBinError(ibin,TMath::Sqrt(newcontent));
  }
  
  return smearedHisto;
}

void pseudoDataGeneration(TH1D* bkgHisto, TH1D* sgnHisto) {
  // initialize random generator
  // seed number taken from clock (microsecond precision):
  struct timeval tp;
  gettimeofday(&tp, 0);
  gRandom = new TRandom3(tp.tv_usec+tp.tv_sec*1000000);

//  bkg_in_pseudoDataEntries = bkgEntries*bkgScale; 
//  sgn_in_pseudoDataEntries = sgnEntries*sgnScale;
//  int bkgEntries_gaus =  int(gRandom->Gaus(bkgEntries));
//  int sgnEntries_gaus =  int(gRandom->Gaus(sgnEntries));
  double bkgEntries_poisson = gRandom->Poisson(bkgEntries); // std::cout << "bkgEntries_poisson*bkgScale: " << bkgEntries_poisson*bkgScale << std::endl;
  double sgnEntries_poisson = gRandom->Poisson(sgnEntries); // std::cout << "sgnEntries_poisson*sgnScale: " << sgnEntries_poisson*sgnScale << std::endl;
  bkg_in_pseudoDataEntries = bkgEntries_poisson*bkgScale;
  sgn_in_pseudoDataEntries = sgnEntries_poisson*sgnScale;

  if ( flag_bkgSmear ) pseudoData->FillRandom(smearedHistoConstruction(bkgHisto, sigma_bkgSmear),Int_t(bkg_in_pseudoDataEntries));
  else pseudoData->FillRandom(bkgHisto,Int_t(bkg_in_pseudoDataEntries));
  if ( flag_sgnSmear ) pseudoData->FillRandom(smearedHistoConstruction(sgnHisto, sigma_sgnSmear),Int_t(sgn_in_pseudoDataEntries));
  else pseudoData->FillRandom(sgnHisto,Int_t(sgn_in_pseudoDataEntries));
  pseudoDataEntries=0;
  for ( int ibin=1; ibin<=nbin_; ibin++ ) {
    //  double contents = bkgHisto->GetBinContent(ibin)+sgnHisto->GetBinContent(ibin);
    //  pseudoData->SetBinContent(ibin,contents);
    double contents = pseudoData->GetBinContent(ibin);
    pseudoDataEntries+=contents;
  }
  //  pseudoDataEntries = int(pseudoData->GetEntries());

  std::cout << "pseudoDataEntries: " << pseudoDataEntries << std::endl;
  std::cout << "bkg_in_pseudoDataEntries: " << bkg_in_pseudoDataEntries << std::endl;
  std::cout << "sgn_in_pseudoDataEntries: " << sgn_in_pseudoDataEntries << std::endl;

  // calculation of error on pseudo-data
  // using an approximation for low statistics
  // -see N. Gehrels APJ, 303 
  // or http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?explain_errors
  for ( int ibin=1; ibin <= nbin_; ++ibin ) {
    double content = pseudoData->GetBinContent(ibin);
    double error = 0.0;
    if ( content < 10. ) error = (1 + TMath::Sqrt(content + 0.75));
    else error = TMath::Sqrt(content);
    pseudoData->SetBinError(ibin,error);
  }
}

void fitHistoGeneration(TH1D* bkgHisto, TH1D* sgnHisto, double par_val, double par_err) {
 double bkgFrac = (1.-par_val)/bkgEntries*pseudoDataEntries;
 //  std::cout << "bkgFrac: " << bkgFrac << std::endl; 
 double sgnFrac = par_val/sgnEntries*pseudoDataEntries;    
 //   std::cout << "sgnFrac: " << sgnFrac << std::endl;  
  totHisto->Add(bkgHisto,bkgFrac);
  totHisto->Add(sgnHisto,sgnFrac);
  // calculation of error on fit histo
  for ( int ibin=1; ibin <= nbin_; ++ibin ) {
    bkg_content = bkgHisto->GetBinContent(ibin);
    sgn_content = sgnHisto->GetBinContent(ibin);
    double bkg_error2 = 0.0;
    double sgn_error2 = 0.0;
    if ( bkg_content < 10. ) bkg_error2 = (1 + TMath::Sqrt(bkg_content + 0.75))*(1 + TMath::Sqrt(bkg_content + 0.75));
    else bkg_error2 = bkg_content;
    if ( sgn_content < 10. ) sgn_error2 = (1 + TMath::Sqrt(sgn_content + 0.75))*(1 + TMath::Sqrt(sgn_content + 0.75));
    else sgn_error2 = sgn_content;
    double error = TMath::Sqrt(sgn_content*sgn_content*par_err*par_err+
			       par_val*par_val*sgn_error2+
			       bkg_content*bkg_content*par_err*par_err+
			       (1.-par_val)*(1.-par_val)*bkg_error2);
    totHisto->SetBinError(ibin,error);
  }
}

void drawCanvas(TCanvas* c, TH1D* fithisto, TH1D* data, TH1D* bkg, TH1D* sgn, 
		double fittedBkgvalue, double fittedBkgerror, 
		double fittedSgnvalue, double fittedSgnerror, 
		double k2, int nf, int ipseudoexp) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);

  fithisto->SetLineColor(7); 
  fithisto->SetFillColor(7); 
  fithisto->SetLineWidth(1); 
  
  data->GetXaxis()->SetTitle("pseudo data");
  data->GetXaxis()->SetTitleSize(0.05);
  data->SetTitle("pseudo data fit");
  data->SetStats(1111111);
  data->SetMarkerStyle(21);
  data->SetMarkerSize(0.8);
  data->SetMarkerColor(1); 
  
  bkg->SetFillColor(0);
  bkg->SetLineWidth(2);
  bkg->SetLineColor(kRed);
  bkg->SetLineStyle(1); 
  
  sgn->SetFillColor(0);
  sgn->SetLineWidth(2);
  sgn->SetLineColor(kBlue);   
  sgn->SetLineStyle(2); 
  
  c->UseCurrentStyle();
  c->cd();
  c->GetPad(0)->SetLogy();
  fithisto->Draw("e2same");
  data->Draw("epsame");
  bkg->Draw("samehisto");
  sgn->Draw("samehisto");
  
  char c2[50];
  char nev[50];
  // draw the legend
  TLegend *legend=new TLegend(0.67,0.76,1.0,1.0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);  
  legend->SetTextFont(72);
  legend->SetTextSize(0.027);
  legend->SetFillColor(0);
  sprintf(nev,"data events: %d",int(pseudoDataEntries));
  legend->AddEntry(data,nev,"lpe");
  sprintf(nev,"#color[2]{bkg: %.2f #pm %.2f}",fittedBkgvalue,fittedBkgerror);
  legend->AddEntry(bkg,nev,"l");
  sprintf(nev,"#color[4]{sgn: %.2f #pm %.2f}",fittedSgnvalue,fittedSgnerror);
  legend->AddEntry(sgn,nev,"l");
  legend->AddEntry(fithisto,"#color[7]{Global Fit}","l");
  sprintf(c2,"#chi^{2}/ndf: %.2f / %d",k2,nf);
  legend->AddEntry(fithisto,c2,"");
  legend->Draw();
  /*
    stringstream cName << "fitCanvasPseudoExp_" << ipseudoexp << ".eps";
    c->Print(TString(cName.str()));
  */


}

void computePull(double expectedParvalue, double constrainedParerror, double fittedParvalue,double fittedParerror, TH1D* significance, TH1D* pull) {
  // significance computation
 double fittedParsignificance = 0.;
 if (fittedParerror != 0.) 
   fittedParsignificance = fittedParvalue/fittedParerror;
 // pull computation
 double fittedParpull = 0.;
 double fittedParmean = fittedParvalue-expectedParvalue;
 // compute fittedParsigmaPull as TMath::Sqrt(fittedParerror*fittedParerror - fittedParvalue*fittedParvalue*constrainedParerror*constrainedParerror); 
 // since it is necessary to teke into account the constraint on the background knowledge in the chi2 function
// double fittedParsigmaPull = TMath::Sqrt(fittedParvalue*fittedParvalue*constrainedParerror*constrainedParerror - 
//					   fittedParerror*fittedParerror                                           ); 
 double fittedParsigma = 0.;
 if ( constrainedParerror != 0. ) fittedParsigma = TMath::Sqrt(constrainedParerror*constrainedParerror - 
							      fittedParerror*fittedParerror              ); 
 else fittedParsigma = fittedParerror;
 if (fittedParsigma != 0.) 
   fittedParpull = fittedParmean/fittedParsigma;
 std::cout << "expectedParvalue: "      << expectedParvalue      << std::endl;
 std::cout << "constrainedParerror: "   << constrainedParerror   << std::endl;
 std::cout << "fittedParvalue: "        << fittedParvalue        << std::endl;
 std::cout << "fittedParerror: "        << fittedParerror        << std::endl;
 std::cout << "fittedParsignificance: " << fittedParsignificance << std::endl;
 std::cout << "fittedParmean: "         << fittedParmean         << std::endl;
 std::cout << "fittedParsigma: "        << fittedParsigma        << std::endl;
 std::cout << "fittedParpull: "         << fittedParpull         << std::endl;
 std::cout << "-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8-8" << std::endl;
 if (fittedParerror != 0) significance->Fill(fittedParsignificance);
 if (fittedParerror != 0) pull->Fill(fittedParpull);
}

// Sum of background and peak function
//extern 
double func(double sgn, double bkg, double *par) {
 double value = par[0]*sgn + (1-par[0])*bkg;
 return value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate chisquare
void chi2function(int &npar, double *grad, double &fval, double *par, int flag) {

  if (flag == 1 ) std::cout << "flag=1" << std::endl;
  //  int pseudoDataEntries = pseudoData->GetEntries();
  chi = 0;
  ndof=0;
  double delta = 0.;      
  for ( int ibin=1; ibin <= nbin_; ++ibin ) {
    sgn_content  = signalHisto->GetBinContent(ibin);                     // std::cout << "sgn_content: "  << sgn_content  << std::endl;
    bkg_content  = backgroundHisto->GetBinContent(ibin);              // std::cout << "bkg_content: "  << bkg_content  << std::endl;
    data_content = pseudoData->GetBinContent(ibin);              // std::cout << "data_content: " << data_content << std::endl;
    //    data_error   = TMath::Sqrt(pseudoData->GetBinContent(ibin)); // std::cout << "data_error: "   << data_error   << std::endl;
    data_error   = pseudoData->GetBinError(ibin); // std::cout << "data_error: "   << data_error   << std::endl;
    if (data_error!=0. && (sgn_content!=0. || bkg_content!=0.) ) {
      delta = (data_content - func(sgn_content/sgnEntries,bkg_content/bkgEntries,par)*pseudoDataEntries)/data_error;
      ++ndof;
    } else { 
//      std::cout << "DATA ERROR = 0!!!";
//      std::cout << "ibin: " << ibin 
//		<< " data_content: " << data_content 
//		<< " sgn_content: " << sgn_content 
//		<< " bkg_content: " << bkg_content << std::endl;
    }
    chi += delta*delta;	
  }
  double constraint = 0.;
  // constraint on background knowledge
  if ( uncert_bkgEntries != 0. ) {
    constraint = (bkg_in_pseudoDataEntries - 
		  (1.-par[0])*pseudoDataEntries)/(uncert_bkgEntries*bkg_in_pseudoDataEntries);
    //  std::cout << "constraint: " << constraint*constraint << std::endl;
  }
  chi += constraint*constraint;

  fval = chi;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate likelihood
void likelihoodfunction(int &npar, double *grad, double &fval, double *par, int flag) {

  like      = -9.;
  ndof      = 0;
  double delta = 0.;      
  for ( int ibin = 1; ibin <= nbin_; ++ibin ) { // loop over bins
    sgn_content  = signalHisto->GetBinContent(ibin);                     // std::cout << "sgn_content: "  << sgn_content  << std::endl;
    bkg_content  = backgroundHisto->GetBinContent(ibin);              // std::cout << "bkg_content: "  << bkg_content  << std::endl;
    data_content = pseudoData->GetBinContent(ibin);              // std::cout << "data_content: " << data_content << std::endl;
    // for pseudo experiments 
    // we must not include the error on the templates 
    // used to fit the 'data' distribution ! 
    data_error   = pseudoData->GetBinError(ibin); // std::cout << "data_error: "   << data_error   << std::endl;
    if (sgn_content!=0. || bkg_content!=0.) {
      delta = func(sgn_content/sgnEntries,bkg_content/bkgEntries,par);
      ++ndof;
      if ( delta > 0. )	like -= data_content*log(delta)-delta;
    }
  } // end loop over bins

  like *= 2;

  fval = like;
  // cout << "fval= " << fval << endl;
}
