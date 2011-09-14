/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #201
// 
// Composite p.d.f with signal and background component
//
// pdf = f_bkg * bkg(x,a0,a1) + (1-fbkg) * (f_sig1 * sig1(x,m,s1 + (1-f_sig1) * sig2(x,m,s2)))
//
//
// 07/2008 - Wouter Verkerke 
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;
#include "TH1F.h"
#include <vector>

void Convolution_variationWithTau()
{
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  // Declare observable x
  RooRealVar x("x","x",0,3);

  double sigmaValue = 0.04;
  double tauValue = 0.02;
  std::vector<double> maxValues;

  int bins = 100;
  TH1F * positionOfTheMaximum = new TH1F("positionOfTheMaximum", "position of the maximum", bins, tauValue, tauValue+1);

  for( int i=0; i<bins; ++i ) {

    // sigmaValue += float(i)/float(bins);
    tauValue += 1./float(bins);
    // std::cout << "float("<<i<<")/float("<<bins<<") = " << float(i)/float(bins) << std::endl;

    // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their paramaters
    RooRealVar mean("mean","mean of gaussians",0);
    RooRealVar sigma("sigma","width of gaussians",sigmaValue);
    RooGaussModel resolModel("sig1","Signal component 1",x,mean,sigma);  
  
    RooRealVar tau("tau", "lifetime", tauValue);
    RooDecay model("sigNP", "Signal component NP", x, tau, resolModel, RooDecay::SingleSided);

    // S a m p l e ,   f i t   a n d   p l o t   m o d e l
    // ---------------------------------------------------

    RooAbsReal * deriv = model.derivative(x, 1);
    Double_t xmax = deriv->findRoot(x,0,3,0);

    std::cout << "for ctau = " << tauValue << ", xmax = " << xmax << std::endl;

    maxValues.push_back(xmax);

    // Plot data and PDF overlaid
    if( i == 1 ) {
      RooPlot* xframe = x.frame(Title("Example of composite pdf=(sig1+sig2)+bkg")) ;
      model.plotOn(xframe) ;
      // Draw the frame on the canvas
      new TCanvas("rf201_composite","rf201_composite",600,600) ;
      gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw();
    }
  }

  std::vector<double>::const_iterator it = maxValues.begin();
  for( int i=1; it != maxValues.end(); ++it, ++i ) {
    positionOfTheMaximum->SetBinContent(i, *it);
  }

  TCanvas * canvas = new TCanvas;
  positionOfTheMaximum->Draw();
  positionOfTheMaximum->GetXaxis()->SetTitle("c#tau");
  positionOfTheMaximum->GetYaxis()->SetTitle("maximum position on x");
  canvas->SaveAs("canvas.pdf");

}
