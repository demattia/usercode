// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>

#include <TPad.h>
#include <TAxis.h>
#include "RooHist.h"
#include "TPaveLabel.h"

#include <TCanvas.h>
#include <TLatex.h>
//#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooCategory.h"
#include "RooWorkspace.h"

using namespace RooFit;
using namespace std;


void defineBackground(RooWorkspace *ws, int& index)
{
  cout << "bkg definition" << endl;

  ws->factory("Exponential::expFunct(BuMass,tau[-4.,-10.,10])");
  ws->factory("EXPR::err('(TMath::Erfc((BuMass-erfMean)/erfWidth))',BuMass,erfMean[5.15,5.1,5.2],erfWidth[0.08,0.001,100])");
//  ws->factory("SUM::bkgPDF(erfFrac[0.3,0.,1]*err,expFunct)");
    ws->factory("Chebychev::che(BuMass,{m[.3,-1,1],c[.5,-1,1]})");
//   ws->factory("SUM::bkgPDF(coeffMix1[0.1,0.,1.]*expFunct,che)");
  ws->factory("SUM::bkgPDF(expFunct)");
//  ws->factory("SUM::bkgPDF(erfFrac[0.3,0.,1]*err,expFunct)");

  return;
}

void defineSignal(RooWorkspace *ws,int& index)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  if (index==0) ws->factory("Gaussian::sigG1(BuMass,meanSig1[5.27,5.2,5.33],sigmaSig1[0.04,0.025,0.05])");
  else ws->factory("Gaussian::sigG1(BuMass,meanSig1[5.27,5.2,5.33],sigmaSig1[0.05,0.025,0.10])");

  ws->factory("Gaussian::sigG2(BuMass,meanSig1,sigmaSig2[0.03,0.,0.2])");

  //Gaussian with same mean as signalG1
  //ws->factory("Gaussian::signalG2OneMean(DBuMass,meanSig1,sigmaSig2)");

  //Crystall Ball
  if (index==0) ws->factory("CBShape::sigCB(BuMass,meanSig1,sigmaSig1,alpha[2.9,0.1,3],enne[1.5,1.5,1.5])");
  else ws->factory("CBShape::sigCB(BuMass,meanSig1,sigmaSig1,alpha[1.5,.1,3],enne[1.5,1.5,1.5])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::twoGaus(coeffGauss[0.5,0.,1.]*sigG1,sigG2)");

  //Sum of Gaussians with same mean
  //ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,sigG2)");

//  ws->factory("Landau::lan(BuMass,lmean[5.35,5.3,-5.9], lsigma[0.03,5,-20])");
//  ws->factory("SUM::sigPDF(coeffLan[0.05,0,-0.1]*lan,twoGaus)");

  ws->factory("Landau::lan(BuMass,lmean[5.357,5.357,5.357], lsigma[0.04885,0.04885,0.04885])");
//  ws->factory("SUM::sigPDF(coeffLan[0.049,0.049,0.049]*lan,twoGaus)");
  ws->factory("SUM::sigPDF(twoGaus)");

  return;
}

void prefitSideband(RooWorkspace *ws, const int barrel)
{

  //cout << "XXXXXXXXXXXXXXX-> index " << barrel << endl;
  RooDataSet *tmpdata;
  RooDataHist *tmpdataBin; 


  if (barrel==0) tmpdata = (RooDataSet*)ws->data("dataBarrel");
  else tmpdata= (RooDataSet*)ws->data("dataEndCap");

  if (barrel==0) tmpdataBin=tmpdata->binnedClone("dataBinBar");
  else tmpdataBin=tmpdata->binnedClone("dataBinEnd");

  RooPlot *mframe = ws->var("BuMass")->frame();

  tmpdata->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(70));
  ws->pdf("bkgPDF")->plotOn(mframe,LineColor(kRed));
  ws->pdf("bkgPDF")->fitTo(*tmpdataBin,Range("left"),SumW2Error(kTRUE));
  ///ws->pdf("bkgPDF")->plotOn(mframe,LineColor(kBlue));
  ws->pdf("bkgPDF")->plotOn(mframe,Components("lan"),LineColor(kGreen));
  //ws->pdf("bkgPDF")->plotOn(mframe,Components("err"),LineColor(kGreen));
  //ws->pdf("bkgPDF")->plotOn(mframe,Components("expFunct"),LineColor(k));
  //ws->var("c1")->setConstant(kTRUE);
  //ws->var("c")->setConstant(kTRUE);
  //ws->var("m")->setConstant(kTRUE);
  //ws->var("coefExp2")->setConstant(kTRUE);

  TCanvas c1;
  c1.cd();  mframe->Draw();

  string out;
  if (barrel==0) out ="Prefit_Barrel.pdf";
  else out ="Prefit_EndCap.pdf";
  c1.SaveAs(out.c_str());

  /*if (barrel==0) ws->var("alpha")->setVal(2.9);
  else ws->var("alpha")->setVal(1.5);
  ws->var("alpha")->setConstant(kTRUE);
  if (barrel==0) ws->var("c1")->setConstant(kTRUE);
  if (barrel==0) ws->var("c2")->setConstant(kTRUE);*/


  /*ws->var("enne")->setVal(1.5);
    ws->var("enne")->setConstant(kTRUE);*/


  
  //ws->var("meanSig1")->setVal(5.27);
  //ws->var("meanSig1")->setConstant(kTRUE);
  
  if (barrel==0) ws->var("sigmaSig1")->setVal(0.034);
  else ws->var("sigmaSig1")->setVal(0.055);
  ws->var("sigmaSig1")->setConstant(kTRUE);

   
  return;
}

void setRanges(RooWorkspace *ws, int &index)
{
  const float DMassMin = 4.5;
  const float DMassMax = 5.55;

  ws->var("BuMass")->setRange("all",DMassMin,DMassMax);
  if (index!=0) ws->var("BuMass")->setRange("left",DMassMin,5.15);
  else ws->var("BuMass")->setRange("left",DMassMin,5.15);
  ws->var("BuMass")->setRange("right",5.5,DMassMax);

  return;
}

void drawResults(RooWorkspace *ws, string& sname, int& index, int& npar)
{

  RooDataSet *data;

  if (index==0) data = (RooDataSet*)ws->data("dataBarrel");
  else data= (RooDataSet*)ws->data("dataEndCap");

  RooAbsPdf *totPDF = ws->pdf("totPDF");
  RooAbsPdf *BkgPDF = ws->pdf("bkgPDF");

  char reducestr[200];

  RooPlot *mframe = ws->var("BuMass")->frame();

  sprintf(reducestr,"Bu mass fit");

  mframe->SetTitle(reducestr);

  data->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(70));
  totPDF->plotOn(mframe,Components("bkgPDF"),LineColor(kRed),LineStyle(kDashed)/*,Normalization(1.0,RooAbsReal::RelativeExpected)*/);
  totPDF->plotOn(mframe,Components("err"),LineColor(kBlue),LineStyle(kDashed)/*,Normalization(1.0,RooAbsReal::RelativeExpected)*/);
  totPDF->plotOn(mframe,Components("lan"),LineColor(kGreen));
  totPDF->paramOn(mframe,Layout(0.2,0.5,0.5),Format("NELU",AutoPrecision(1)),Label("")); 
  totPDF->plotOn(mframe,LineColor(kBlack)/*,Normalization(1.0,RooAbsReal::RelativeExpected)*/);
/*
  RooPlot* rframe = ws->var("BuMass")->frame();
  RooHist* hresid = mframe->pullHist();
  hresid->SetLineColor(2);
  rframe->addPlotable(hresid,"P") ;
  //rframe->addObject(hresid);
  int nresid=6;
  rframe->SetTitle("");
  rframe->GetXaxis()->SetTitle("");
  rframe->GetYaxis()->SetLabelSize(0.12);
  rframe->GetYaxis()->SetNdivisions(2*nresid+1);
  rframe->SetAxisRange(-1*nresid,nresid,"Y");
  //rframe->SetMarkerStyle(7);
  //rframe->Draw();

  TCanvas c1;
  c1.cd();  

  TPaveLabel *t2 = new TPaveLabel(0.01,0.1,0.07,0.9, "Residual pull", "brNDC"); 
  t2->SetFillColor(0);
  t2->SetTextAngle(90);
  t2->SetTextSize(0.2);
  t2->SetBorderSize(0);

  TPad *p1 = new TPad("p1","p1",0.1,0.25,0.9,0.901);// ,0,0,0); 
  //p1->SetBottomMargin(0.);
  //p1->SetBorderMode(0); 
  p1->Draw(); 
  //p1->cd();
  // mframe->Draw();
   
  TPad *p2 = new TPad("p2","p2",0.1,0.1,0.9,0.24); 
  //p2->SetTopMargin(0.);    
  //p2->SetBorderMode(0); 
  //p2->SetTicks(1,0); 
  p2->Draw(); 
  //p2->cd();
  //rframe->Draw();

  p1->cd();
  mframe->Draw();
  p2->cd();
  rframe->Draw();
  t2->Draw();
*/
  
  //mframe->Draw();
  //rframe->Draw();

  double chi=mframe->chiSquare(npar);

  TCanvas c1;
  c1.cd();  mframe->Draw();

  char sig[80], ms[80],SovB[80];
  double errsig = ws->var("NSig")->getError();
  double Nsig   = ws->var("NSig")->getVal();
  sprintf(sig,"N_{Sig}= %5.0f #pm %3.0f", Nsig,errsig);

  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(13);
  t->SetTextFont(63);
  t->SetTextSizePixels(20);
  t->DrawLatex(0.55,0.85,sig);
  t->SetTextSizePixels(20);
  double Nbkg = ws->var("NBkg")->getVal();
  double coeffExp=0;
  //coeffExp= ws->var("coefExp")->getVal();
  double sigma = 1000*ws->var("sigmaSig1")->getVal();
  double esigma = 1000*ws->var("sigmaSig1")->getError();
  double m0 = ws->var("meanSig1")->getVal();
  double em0 =0;// ws->var("meanSig1")->getError();
  sprintf(ms,"m_{0} = %6.4f #pm %5.4f GeV/c^{2}", m0,em0);
  t->DrawLatex(0.55,0.80,ms);
  sprintf(ms,"#sigma = %4.1f #pm %3.1f MeV/c^{2}",sigma,esigma);
  t->DrawLatex(0.55,0.75,ms);
  double a=ws->var("alpha")->getVal();
  double n=ws->var("enne")->getVal();


  std::cout << " sigma= " << sigma << " mean= " << m0 << " coeff exp= " << coeffExp << " alpha= " << a << " n= " << n << std::endl;
  std::cout << "Nsig= " << Nsig << " Nbkg= " << Nbkg << std::endl;


  //ws->var("DMass")->setRange("integral",m0-2*(sigma/1000),m0+2*(sigma/1000));

  ws->var("BuMass")->setRange("integral",5.2,5.35);
  
  RooAbsReal* N= ws->pdf("bkgPDF")->createIntegral(RooArgSet(*ws->var("BuMass")),NormSet(RooArgSet(*ws->var("BuMass"))),Range("integral"));
  RooAbsReal* S= ws->pdf("sigCB")->createIntegral(RooArgSet(*ws->var("BuMass")),NormSet(RooArgSet(*ws->var("BuMass"))),Range("integral"));

  cout << "B= " << N->getVal() << " S= " << S->getVal() << endl;
  double SN=(Nsig*S->getVal())/(Nbkg*N->getVal());

  sprintf(SovB,"S/B= %5.2f",SN);
  t->DrawLatex(0.55,0.7,SovB);

  sprintf(SovB,"#Chi^{2}/d.o.f = %5.2f",chi/npar);
  t->DrawLatex(0.55,0.65,SovB);

  string out;
  if (index==0) out ="MassFitBu_Barrel_"+sname+".png";
  else out ="MassFitBu_EndCap_"+sname+".png";
  c1.SaveAs(out.c_str());

  if (index==0) out ="MassFitBu_Barrel_"+sname+".root";
  else out ="MassFitBu_EndCap_"+sname+".root";
  c1.SaveAs(out.c_str());

  return;

}

string cleanName(string& fname){

  string fnameOld=string(fname);
  string fnameNew;

  size_t pos=0;

  pos = fnameOld.find("/");
  fnameNew=fnameOld.substr(pos+1);

  while(pos != string::npos){
    pos = fnameNew.find("/");
    fnameNew=fnameNew.substr(pos+1);
  }
  
  pos=fnameNew.find(".root");
  fnameNew=fnameNew.substr(0,pos);
  //cout << fnameNew << endl;
  return fnameNew;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol)
{
  Nsig   = ws->var("NSig")->getVal();
  errSig = ws->var("NSig")->getError();
  //const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  //const double sigmaSig2 = ws->var("sigmaSig2")->getVal();

  resol = sigmaSig1;
  //resol = (coeffGauss*coeffGauss*sigmaSig1 + (1-coeffGauss)*(1-coeffGauss)*sigmaSig2) / (coeffGauss*coeffGauss + (1-coeffGauss)*(1-coeffGauss)); 

  return;
}

int main(int argc, char* argv[])
{
 
  //gROOT->ProcessLine(".L setTDRStyle_modified.C");
  //gROOT->ProcessLine("setTDRStyle()");

  char *filename;
  bool sidebandPrefit = false;
  
  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];
    
    switch(pchar[0]){
      
    case '-':{
      
      switch(pchar[1]){
	
      case 'f':{
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
	
      case 's':{
        sidebandPrefit = true;
        cout << "Sideband pre-fitting activated" << endl;
        break;
      }
      }
    }
    }
  }

 

  /*string strname=string(filename);
  size_t p;
  
  if (strname.find("Ds") != string::npos){
    p=strname.find("Ds");
    strname.replace(p,2,"D");
  }
  cout << filename << "  e dopo " << strname << endl;*/


  TFile fIn(filename);
  fIn.cd();

  string strname=string(filename);
  strname=cleanName(strname);

  for (int i=0; i<2; i++){

    //if (i==0) continue;
    RooWorkspace *ws = new RooWorkspace("ws");

    if (i==0) cout << "Barrel" << endl;
    else cout << "EndCap" << endl;

    RooDataSet *reddata=new RooDataSet(); 
    if (i==0) reddata= (RooDataSet*)fIn.Get("dataBarrel");
    else reddata= (RooDataSet*)fIn.Get("dataEndCap");
    reddata->Print(); 

//    reddata->reduce(Cut("BuMass>5.0 && BuMass<5.55"));

    ws->import(*reddata);

    cout<<"hello1"<<endl;

    ws->var("BuMass")->setBins(70);
    cout<<"hello2"<<endl;

    RooDataHist *reddataBin;
    /*if (i==0) reddataBin= new RooDataHist("dataBinBar","dataBinBar",RooArgSet(*(ws->var("BuMass"))),*reddata);
      else reddataBin= new RooDataHist("dataBinEnd","dataBinEnd",RooArgSet(*(ws->var("BuMass"))),*reddata);*/
    if (i==0) reddataBin=reddata->binnedClone("dataBinBar");
    else reddataBin=reddata->binnedClone("dataBinEnd");

    //ws->import(*reddataBin);
    setRanges(ws,i);

    //DEFINE SIGNAL AND BACKGROUND
    defineSignal(ws, i);
    defineBackground(ws, i);

    // Total PDF (signal CB+Gauss)
    if (i==0) ws->factory("SUM::totPDF(NSig[500000.,300000.,800000.]*sigPDF,NBkg[5000000.,1000000.,10000000.]*bkgPDF)");
    else ws->factory("SUM::totPDF(NSig[150000.,90000.,200000.]*sigPDF,NBkg[1500000.,1000000.,10000000.]*bkgPDF)");
    // ws->factory("SUM::totPDF(NSig[5000.,10.,100000.]*sigCB,NBkg[5000.,10.,100000.]*expFunct)");

    if (sidebandPrefit) prefitSideband(ws,i);
    cout<<"hello3"<<endl;

//    ws->var("BuMass")->setRange("all",5.,5.5);
    //ws->var("coefExp")->setVal(-2.02);
    //ws->var("coefExp")->setConstant(kTRUE);
    ws->var("BuMass")->setRange("fit",5.0,5.55);

    RooFitResult *rfr=ws->pdf("totPDF")->fitTo(*reddataBin,Extended(1),Range("fit"),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
    rfr->Print();
    cout<<"hello6"<<endl;
    int npar=rfr->floatParsFinal().getSize();
    cout<<"hello4"<<endl;

    drawResults(ws,strname,i,npar);
    cout<<"hello5"<<endl;

    double NSig, errSig,resol;
    printResults(ws,NSig,errSig,resol);

    cout << endl << "D yield:"<< endl;
    cout << "Fit : " << NSig << " +/- " << errSig << endl;
    cout << "Resolution : " << resol*1000. << endl; 
    
    
    delete ws;
    //delete reddataBin;
    
  }
  return 1;
}
