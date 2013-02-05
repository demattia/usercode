#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooThresholdCategory.h"
#include "RooBinningCategory.h"
#include "RooWorkspace.h"


#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TPaveLabel.h"

using namespace RooFit;

void readData(RooWorkspace&, TString);
void buildPdf(RooWorkspace&);

void fitYield(RooWorkspace&, TString, RooAbsData* d=0);
void plotMass(TString, RooWorkspace&, RooAbsData* d=0, RooAbsPdf* p=0);

void xsectTot(RooWorkspace&);
void xsectDif(RooWorkspace&, bool);

void getInitialParValues(RooWorkspace&);
void setDefaultFitParValues(RooWorkspace&,RooAbsData* d=0);

vector<double> initParVal_;
vector<string> initParName_;

#define mmin_ 8.0 
#define mmax_ 12.0
const TString figs_("figs/"); /* store location for produced file output */
const TString res_("results"); /* root folder in file for storing results */
const TString plt_("plots"); /* root folder in file for storing plots */
const TString dirname_("upsilonYield"); /* tree location in input file */
const TString ext_(".gif"); /* save plots format */
const double lumi_ = 1870; /*sample lumi in 1/nb*/

void oniafitter(TString finput = "upsilonYieldWeighted.root", TString foutput = "xsection.root" ) {
  
  /// roofit workspace to manage fitting info
  RooWorkspace* ws = new RooWorkspace("ws","upsilon mass");

  /// create fitting model
  buildPdf(*ws);

  /// read the data
  readData(*ws,finput);
  plotMass("data",*ws,ws->data("data"));

  /// cache the initial parameter
  getInitialParValues(*ws);

  /// open output root file
  TFile file(foutput,"update");
  gDirectory->mkdir(res_);
  gDirectory->mkdir(plt_);
  gDirectory->Cd(res_);

  /* PLEASE UNCOMMENT BELOW ONLY THOSE FUNCTION CALLS
     YOU WANT TO BE EXECUTED EACH TIME THE MACRO IS RUN */

  /// fit raw yield
  fitYield(*ws,"total");
  plotMass("mass_raw",*ws);

  /// fit corrected yield 
  //((RooDataSet*)ws->data("data"))->setWeightVar(*ws->var("weight"));
  //fitYield(*ws,"total");
  //plotMass("mass_wei",*ws);

  /// total cross section
  //xsectTot(*ws);

  /// cross section vs pT
  //xsectDif(*ws,true);

  ///save workspace
  //ws->Print();
  //ws->Write();
  //ws->writeToFile("ws_fitter.root") ;
  //gDirectory->Add(ws) ;
  file.Close();
}

/* total cross section
   retrieve and normalize global fit results
*/
void xsectTot(RooWorkspace& w) {

  cout << "computing total cross section...\n" << std::flush;
  gDirectory->Cd("../"+res_);

  ///raw yield
  RooFitResult* fitres_r  = (RooFitResult *) gROOT->FindObject("fit_result_raw_total");
  double nsigVal_r        = ((RooRealVar*)fitres_r->floatParsFinal().find("nsig1"))->getVal();
  double nsigErr_r        = ((RooRealVar*)fitres_r->floatParsFinal().find("nsig1"))->getError();

  ///corrected yield
  RooFitResult* fitres_w  = (RooFitResult *) gROOT->FindObject("fit_result_corrected_total");
  double nsigVal_w        = ((RooRealVar*)fitres_w->floatParsFinal().find("nsig1"))->getVal();
  double nsigErr_w        = ((RooRealVar*)fitres_w->floatParsFinal().find("nsig1"))->getError();

  /// normalize
  double xsec  = nsigVal_w / lumi_;
  double xsecE = nsigErr_w / lumi_;

  /// check & print
  double mean_weight = ((RooDataSet*)w.data("data"))->mean(*w.var("weight"));
  cout << "raw  yield " << nsigVal_r << "+-" << nsigErr_r << "\n";
  cout << "cross check raw * weight(" <<  mean_weight << ") = " 
       << nsigVal_r*mean_weight << "+-" << nsigErr_r*mean_weight << "\n";
  cout << "weighted yield " << nsigVal_w << "+-" << nsigErr_w << "\n";
  //cout << "cross section:" << xsection << "+-" << xsection_err << " nb\n";

  /// save value as root object
  TGraphAsymmErrors xst;
  xst.SetPoint(0,0, xsec);
  xst.SetPointError(0,0.1,0.1, xsecE, xsecE);
  gDirectory->Cd("../"+plt_);
  xst.SetName("xsection_total");
  xst.Write();

}

/* differential cross section, sigma vs pT
-split the dataset in pT regions
-fit the sub-samples, extract the yield, and normalize   
*/
void xsectDif(RooWorkspace& w, bool dofit) {

  RooAbsPdf*  pdf   = w.pdf("pdf");
  RooDataSet* data  = (RooDataSet*)w.data("data");
  RooRealVar* upsPt = w.var("upsPt");

  /// assign upsilon-pT range categories
  const int nptbins = 6;
  double UpsPtBinEdges[nptbins+1] = {0,2.5,4.5,7,10,15,20};
  double UpsPtBinCenter[nptbins], UpsPtBinEdgeLo[nptbins], UpsPtBinEdgeHi[nptbins];
  RooThresholdCategory ptRegion("ptRegion","region of pt",*upsPt) ;
  for(int i=0; i<nptbins; i++) {
    TString reg = TString::Format("PtBin%d",i+1);    
    ptRegion.addThreshold(UpsPtBinEdges[i+1],reg) ;    
  }
  data->addColumn(ptRegion);

  /// process each pt subsample
  RooDataSet* dataPt;
  for(int i=0; i<nptbins; i++) {
    cout << "processing subsample " << i << "\n" << std::flush;
    dataPt = (RooDataSet*) data->reduce(TString::Format("ptRegion==ptRegion::PtBin%d",i+1));
    dataPt->setWeightVar(*w.var("weight"));
    UpsPtBinCenter[i] =  dataPt->mean(*w.var("upsPt")); 
    UpsPtBinEdgeLo[i] =  UpsPtBinCenter[i] - UpsPtBinEdges[i];
    UpsPtBinEdgeHi[i] = -UpsPtBinCenter[i] + UpsPtBinEdges[i+1];
    //cout << "NTOTAL sample " << i << " = " << dataPt->sumEntries() << "\n";
    //if(i!=0) continue; /* re-fit sub-sample */
    if(dofit) {
      fitYield(w, TString::Format("pt%d",i), dataPt);
      plotMass(TString::Format("massfit_pt%d",i),w, dataPt,pdf);  
    }
  }


  /// store yield
  double UpsYieldPt  [nptbins];
  double UpsYieldPt_e[nptbins];
  gDirectory->Cd("../"+res_);
  RooFitResult* fitresPt;
  for(int i=0; i<nptbins; i++) {
    fitresPt = (RooFitResult *) gROOT->FindObject(TString::Format("fit_result_corrected_pt%d",i));
    UpsYieldPt  [i] = ((RooRealVar*)fitresPt->floatParsFinal().find("nsig1"))->getVal();
    UpsYieldPt_e[i] = ((RooRealVar*)fitresPt->floatParsFinal().find("nsig1"))->getError();
  }

  /// normalize yield
  for(int i=0; i<nptbins; i++) {
    double binw =  UpsPtBinEdges[i+1]-UpsPtBinEdges[i];
    UpsYieldPt  [i] *= 1./lumi_/binw;
    UpsYieldPt_e[i] *= 1./lumi_/binw;
  }

  /// check
  for(int i=0; i<nptbins; i++) {
    cout << "UpsPtBinCenter [" << i << "]=" << UpsPtBinCenter[i] 
	 << " UpsYieldPt:" << UpsYieldPt[i]
	 << "+-" << UpsYieldPt_e[i]
	 << "\n";
  }

  /// plot cross section
  TGraphAsymmErrors *fitgr = new 
    TGraphAsymmErrors(nptbins,   UpsPtBinCenter,   UpsYieldPt,  
		      UpsPtBinEdgeLo, UpsPtBinEdgeHi, 
		      UpsYieldPt_e,UpsYieldPt_e);
  gDirectory->Cd("../"+plt_);
  TCanvas c; c.cd();
  double x[2] = {0, 18};
  double y[2] = {0.,10};
  TGraph frame(2,x,y);
  TString ytitle("#frac{d#sigma}{dp_{T}} . BR(#Upsilon#rightarrow#mu#mu) [nb/GeV]");
  frame.SetTitle( "" );
  frame.GetXaxis()->SetTitle("p_{T} (#Upsilon) [GeV]");
  frame.GetYaxis()->SetTitle(ytitle);
  frame.Draw("AP");
  fitgr->Draw("P");
  //gPad->SetLogy();
  fitgr->SetName("xsection_pt");
  fitgr->Write();
  c.SaveAs(figs_+"xsecdiff"+ext_);
}


/* plot the mass distribution and fitted model
 */
void plotMass(TString hname, RooWorkspace& w, RooAbsData *data, RooAbsPdf* pdf) {

  gROOT->SetStyle("Plain");

  //float yield[4];
  bool dataonly = (data && !pdf);
  if(!data)
    data = (RooDataSet*)w.data("data");
  if(!pdf)
    pdf = (RooAbsPdf*)w.pdf("pdf");
  
  RooPlot* frame = w.var("invariantMass")->frame() ;
  data->plotOn(frame);
  if(!dataonly) {
    pdf->paramOn(frame,Layout(0.6));
    pdf->plotOn(frame,Components("bkg"),
		LineStyle(1),LineWidth(2),LineColor(16)) ; 
    pdf->plotOn(frame);
  }
  data->plotOn(frame);
  frame->SetTitle( "" );
  frame->GetXaxis()->SetTitle("#mu#mu invariant mass [GeV]");
  TCanvas c; c.cd();
  //gPad->SetLogy();
  frame->Draw();
  c.SaveAs(figs_+hname+ext_);
  gDirectory->Cd("../"+plt_);
  frame->SetName(hname);
  frame->Write(hname);
}

/* fit the mass distribution
- save fit results to file
 */
void fitYield(RooWorkspace& w, TString name, RooAbsData* data) {
  cout << "fitting the upsilon mass...\n" << std::flush;

  if(!data)
    data = (RooDataSet*)w.data("data");

  /// reset the fit parameters, retune the yields
  setDefaultFitParValues(w,data);

  RooFitResult* fitres = w.pdf("pdf")->fitTo(*data, Save(), Extended(kTRUE)); //SumW2Error(data->isWeighted()));//,Range(mmin_,mmax_));
			  
  TString fres_n("fit_result_");
  fres_n += data->isWeighted()?"corrected":"raw";
  fres_n += ("_"+name);
  fitres->SetName(fres_n);
  gDirectory->Cd("../"+res_);
  fitres->Write();
  fitres->Print();
  cout << "/tsaved results in " << fres_n << "\n" << std::flush;
  //ws->pdf("pdf")->fitTo(ws->data("data")) ;//, Save(), Extended(kTRUE));
}

/* read the data from ttree
 */
void readData(RooWorkspace& w, TString fname) {
  TFile f(fname,"read");
  TString dirname("upsilonYield"); /* tree location in input file */
  gDirectory->Cd(fname+":/"+dirname);
  TTree* theTree     = (TTree*)gROOT->FindObject("upsilonYield");
  RooRealVar* mass   = w.var("invariantMass");
  RooRealVar* upsPt  = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,50,"GeV"); 
  RooRealVar* weight = new RooRealVar("weight",  "weight"  ,0,100);

  RooDataSet* data0   = new RooDataSet("data","data",theTree,RooArgSet(*mass,*upsPt,*weight)); 
				       
  RooDataSet* data = ( RooDataSet*) data0;//->reduce(EventRange(0,1000),Cut("invariantMass<9.8"));
  w.import(*data);
  //data->Print();
  f.Close();
}

/* define the fit model
+ signal
for each upsilon 1S/2S/3S: double gaussian, or gaussian + 'crystal ball' 
at each peak, the two gaussians have common center
default fitting parameters:
- mass of upsilon 1S floats
- mass shift between three peaks are fixed to pdg (or common scale factor allowed to float)
- the three peaks have a common shape 
- relative rates of three peaks float
+ background
- second order polynominal
*/
void buildPdf(RooWorkspace& w) {

  RooRealVar* mass  = new RooRealVar("invariantMass","#mu#mu mass",mmin_,mmax_,"GeV"); 

  const double M1S = 9.46;   //upsilon 1S pgd mass value
  const double M2S = 10.02;  //upsilon 2S pgd mass value
  const double M3S = 10.35;  //upsilon 3S pgd mass value

  RooRealVar *mean    = new RooRealVar("mass_mean","#Upsilon mean",M1S,M1S-0.5,M1S+0.5,"GeV");
  RooRealVar *shift21 = new RooRealVar("shift2","mass diff #Upsilon(1,2S)",M2S-M1S);
  RooRealVar *shift31 = new RooRealVar("shift3","mass diff #Upsilon(1,3S)",M3S-M1S);

  RooRealVar *mscale  = new RooRealVar("mscale","mass scale factor",1.,0.7,1.3);
  mscale->setConstant(kTRUE); /* the def. parameter value is fixed in the fit */

  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",
  					    RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0+@1*@2",
					    RooArgList(*mean,*mscale,*shift21));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0+@1*@2",
					    RooArgList(*mean,*mscale,*shift31));

  RooRealVar *sigma1 = new RooRealVar("sigma1","Sigma_1",0.06,0.01,0.30);
  RooRealVar *sigma2 = new RooRealVar("sigma2","Sigma_2",0.10,0.01,0.30);

  /// to describe final state radiation tail on the left of the peaks
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",1.5,0.2,4);
  RooRealVar *npow   = new RooRealVar("npow","power order",2);
  npow ->setConstant(kTRUE);
  alpha->setConstant(kTRUE);

  /// relative fraction of the two peak components 
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.3,0.,1.);

  /// Upsilon 1S
  RooGaussian *gauss1S1 = new RooGaussian("gauss1S1","1S Gaussian_1",*mass,*mean1S,*sigma1);
  //RooGaussian *gauss1S2 = new RooGaussian("gauss1S2","1S Gaussian_2",*mass,*mean1S,*sigma2);
  RooCBShape  *gauss1S2 = new RooCBShape ("gauss1S2", "FSR cb 1s", 
					  *mass,*mean1S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig1S      = new RooAddPdf  ("sig1S","1S mass pdf",
					  RooArgList(*gauss1S1,*gauss1S2),*sigmaFraction);

  /// Upsilon 2S
  RooGaussian *gauss2S1 = new RooGaussian("gauss2S1","2S Gaussian_1",*mass,*mean2S,*sigma1);
  //RooGaussian *gauss2S2 = new RooGaussian("gauss2S2","2S Gaussian_2",*mass,*mean2S,*sigma2);
  RooCBShape  *gauss2S2 = new RooCBShape ("gauss2S2", "FSR cb 2s", 
					  *mass,*mean2S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig2S      = new RooAddPdf  ("sig2S","2S mass pdf",
					  RooArgList(*gauss2S1,*gauss2S2),*sigmaFraction);

  /// Upsilon 3S
  RooGaussian *gauss3S1 = new RooGaussian("gauss3S1","3S Gaussian_1",*mass,*mean3S,*sigma1);
  RooGaussian *gauss3S2 = new RooGaussian("gauss3S2","3S Gaussian_2",*mass,*mean3S,*sigma2);
  //RooCBShape  *gauss3S2 = new RooCBShape ("gauss3S2", "FSR cb 3s", 
  //					 *mass,*mean2S,*sigma2,*alpha,*npow); 
  RooAddPdf *sig3S      = new RooAddPdf  ("sig3S","3S mass pdf",
					  RooArgList(*gauss3S1,*gauss3S2),*sigmaFraction);

  /// Background
  RooRealVar *bkg_a1  = new RooRealVar("bkg_a1", "background a1", 0, -1, 1);
  RooRealVar *bkg_a2  = new RooRealVar("bkg_a2", "background a2", 0, -1, 1);
  RooAbsPdf  *bkgPdf  = new RooChebychev("bkg","linear background",
					 *mass, RooArgList(*bkg_a1,*bkg_a2));

  /// Combined pdf
  int nt = 30000;
  RooRealVar *nsig1 = new RooRealVar("nsig1","signal 1s yield",nt*0.25*0.6,0,10*nt); 
  RooRealVar *nsig2 = new RooRealVar("nsig2","signal 2s yield",nt*0.25*0.3,0,10*nt); 
  RooRealVar *nsig3 = new RooRealVar("nsig3","signal 3s yield",nt*0.25*0.1,0,10*nt); 
  RooRealVar *nbkgd = new RooRealVar("nbkgd","brackground yield",nt*0.75,0,10*nt); 
  RooAbsPdf  *pdf   = new RooAddPdf ("pdf","total signal+background pdf", 
				    RooArgList(*sig1S,*sig2S,*sig3S,*bkgPdf),
				    RooArgList(*nsig1,*nsig2,*nsig3,*nbkgd));
 
  w.import(*pdf);

}


/* cache the default fit parameters
 */
void getInitialParValues(RooWorkspace& w) {
  RooArgSet* pars = (RooArgSet*)w.pdf("pdf")->getParameters(*w.data("data"));//->selectByAttrib("Constant",kFALSE);
  TIterator* coefIter = pars->createIterator() ;
  RooRealVar* coef ;
  while((coef = (RooRealVar*)coefIter->Next())) {
    cout << "caching parameter " << coef->GetName() << "\t" << flush;
    coef->Print();
    initParVal_.push_back(coef->getVal());
    initParName_.push_back(coef->GetName());
  }  
}

/* reset initial parameters of the fit pdf
- to their default values
- if the dataset is specified, retune the yield to counting-based estimates
 */
void setDefaultFitParValues(RooWorkspace& w,RooAbsData* data) {

  RooArgSet* pars = (RooArgSet*)w.pdf("pdf")->getParameters(*w.data("data"));//->selectByAttrib("Constant",kFALSE);
  RooRealVar* par = 0;
  for(int i=0; i<pars->getSize();i++) {
    par = (RooRealVar*)pars->find(initParName_[i].data());
    cout << "reseting parameter " << par->GetName() << " " << par->getVal() << " ->\t" << flush;
    par->setVal(initParVal_[i]);
    par->Print();
  }
  if(data) {
    double nsig = data->sumEntries("invariantMass>9.2 & invariantMass<9.7");
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig1"))->setVal(nsig*0.6);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig2"))->setVal(nsig*0.3);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nsig3"))->setVal(nsig*0.1);
    ((RooRealVar*)w.pdf("pdf")->getParameters(data)->find("nbkgd"))->setVal(nsig*3);
  }
}
