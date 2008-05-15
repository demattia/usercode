#ifndef EFFICIENCYCANVAS_H
#define EFFICIENCYCANVAS_H

#include <vector>
#include <string>

#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"

using namespace std;

/**
 * Template class. The first argument is the type of the MultiHisto class, the second is the type of the histograms <br>.
 * Takes histograms from two different MultiHisto objects and draws them stacked and with a legend. <br>
 * By default they are not normalized and the option used for the stack is "nostack". <br>
 * The number of histograms per canvas can be specified as the fifth parameter to the constructor. By default 4 is used.
 * A subgroup of the total multiHisto histogram vector can be selected by giving the numbers firstHisto and lastHisto as
 * the sixth and sevent parameters to the constructor. The default are 0 for both, where 0 for the lastHisto means all the
 * the end of the vector. So by default all the vector is considered. <br>
 * The canvases are created and written in the Write method.
 */

template<class T1, class T2>
class EfficiencyCanvas {
 public:

  /**
   * Gets the same arguments as the TCanvas constructor, plus extra arguments to specify the number of histograms per
   * canvas and indeces to delimit a subgroup of histograms. <br>
   * The last paramenter determines if the histograms are to be normalized. Deafault is false.
   */
  EfficiencyCanvas( const TString& CANVASTITLE,
		    const TString& CANVASNAME,
		    const int      WIDTH          = 1000,
		    const int      HEIGHT         = 800,
		    const int      HISTOINDEX,
		    const bool     NORMALIZE      = false );

  ~EfficiencyCanvas();

  /// Gets two MultiHisto object pointers (can be MultiTH1F or MultiTProfile) 
  void Fill( const T1 * MULTIHISTO1,
	     const T1 * MULTIHISTO2 );

  /// Gets the option for the stack. Default is "nostack".
  void Write( const TString & OPTION = "nostack" );

 private:
  TString canvasName_;
  TString canvasTitle_;
  int width_;
  int height_;
  int histoIndex_;
  bool normalize_;
  vector<T2*> vecMultiHisto1_;
  vector<T2*> vecMultiHisto2_;
};

// Methods
// -------

template<class T1, class T2>
EfficiencyCanvas<T1, T2>::EfficiencyCanvas( const TString& CANVASNAME,
					    const TString& CANVASTITLE,
					    const int      WIDTH,
					    const int      HEIGHT,
					    const int      HISTOINDEX,
					    const bool     NORMALIZE)
{
  canvasName_     = CANVASNAME;
  canvasTitle_    = CANVASTITLE;
  height_         = HEIGHT;
  width_          = WIDTH;
  histoIndex_     = HISTOINDEX;
  normalize_      = NORMALIZE;
}

template<class T1, class T2>
EfficiencyCanvas<T1, T2>::~EfficiencyCanvas()
{
}

template<class T1, class T2>
void EfficiencyCanvas<T1, T2>::Fill( const T1 * MULTIHISTO1,
				     const T1 * MULTIHISTO2 )
{
  vecMultiHisto1_ = MULTIHISTO1->multiHistos();
  vecMultiHisto2_ = MULTIHISTO2->multiHistos();
}

template<class T1, class T2>
void EfficiencyCanvas<T1, T2>::Write( const TString & OPTION )
{

  // Create the new canvases
    TCanvas * canvas = new TCanvas( canvasName_, canvasTitle_, width_, height_ );
    canvas->Divide(2, 2);

  // Loop on the histograms
  //  typename vector<T2*>::const_iterator histo1_it = vecMultiHisto1_.begin();
  //  typename vector<T2*>::const_iterator histo2_it = vecMultiHisto2_.begin();

  // Loop only on selected histograms

    T2* histo1_ptr = vecMultiHisto1_[histoIndex_];
    T2* histo2_ptr = vecMultiHisto2_[histoIndex_];

    // Generate the name for the stack
    string histoName( histo1_ptr->GetName() );
    // Take the index of the last character in the string matching the underscore.
    // In MultiHisto the index is added as _num. This gives us the part we need to
    // create the index to append to the stack name and title.
    TString histoIndex = histoName.substr( histoName.find_last_of('_') );
    TString eff = "efficiency " + canvasName_;

    THStack * stackVariable   = new THStack( canvasName_, canvasTitle_ );
    THStack * stackEfficiency = new THStack( eff, eff );
    TLegend * legendVariable   = new TLegend( 0.55, 0.65, 0.76, 0.82 );
    TLegend * legendEfficiency   = new TLegend( 0.55, 0.65, 0.76, 0.82 );
    TLegend * legendSBeff  = new TLegend( 0.55, 0.65, 0.76, 0.82 );
    TLegend * legendS2Beff = new TLegend( 0.55, 0.45, 0.90, 0.82 );

    int bin_ = histo1_ptr->GetNbinsX();
    float xmin_ = histo1_ptr->GetXaxis()->GetXmin();
    float xmax_ = histo1_ptr->GetXaxis()->GetXmax();
    TH1F* eff_histo1_ptr = new TH1F(histo1_ptr->GetName(),histo1_ptr->GetTitle(),bin_,xmin_,xmax_);
    TH1F* eff_histo2_ptr = new TH1F(histo2_ptr->GetName(),histo2_ptr->GetTitle(),bin_,xmin_,xmax_);
    TH1F* ratioSBeff     = new TH1F("ratioS_B","signal/background distribution",bin_,xmin_,xmax_);
    TH1F* ratioS2Beff    = new TH1F("ratioS_sqrtB","signal^2/background distribution",bin_,xmin_,xmax_);


    double entries1 = histo1_ptr->Integral(1,bin_);
    double entries2 = histo2_ptr->Integral(1,bin_);
    for(int ibin=1;ibin<=bin_;++ibin) {
      double ientries1 = histo1_ptr->Integral(ibin,bin_); 
      double ientries2 = histo2_ptr->Integral(ibin,bin_); 
      double sgnEff = ientries1/entries1;
      double bkgEff = ientries2/entries2;
      double sgnEff2 = sgnEff*sgnEff;
      double bkgEff2 = bkgEff*bkgEff;
      eff_histo1_ptr->SetBinContent(ibin,sgnEff);
      eff_histo2_ptr->SetBinContent(ibin,bkgEff);
      // binomial statistic errors
      double iSgnError = 1.;
      double iBkgError = 1.;
      if (entries1!=0) iSgnError = sqrt(sgnEff*(1-sgnEff)/entries1);
      if (entries2!=0) iBkgError = sqrt(bkgEff*(1-bkgEff)/entries2);
//      // poisson statistic errors
//      if (ientries1!=0) iSgnError = sgnEff*sqrt(1/ientries1+1/entries1);
//      if (ientrie2!=0) iBkgError = bkgEff*sqrt(1/ientries2+1/entries2);
//      // root errors
//      float iSgnError = eff_histo1_ptr->GetBinError(ibin);
//      float iBkgError = eff_histo2_ptr->GetBinError(ibin);
      eff_histo1_ptr->SetBinError(ibin,iSgnError);
      eff_histo2_ptr->SetBinError(ibin,iBkgError);
      double iSgnError2 = iSgnError*iSgnError;
      double iBkgError2 = iBkgError*iBkgError;
      double ierror = 1.;
      if(bkgEff!=0) {
	ratioSBeff->SetBinContent(ibin,sgnEff/bkgEff);
	ierror = sqrt(1/bkgEff2*iSgnError2+sgnEff2/(bkgEff*bkgEff2)*iBkgError2);
      } else ratioSBeff->SetBinContent(ibin,1);
      ratioSBeff->SetBinError(ibin,ierror);
      if(bkgEff!=0) {
	ratioS2Beff->SetBinContent(ibin,sgnEff*sgnEff/bkgEff);
	ierror =  sqrt(4*sgnEff2/bkgEff2*iSgnError2+sgnEff2*sgnEff2/(bkgEff2*bkgEff2)*iBkgError2);
      } else ratioS2Beff->SetBinContent(ibin,1);      
      ratioS2Beff->SetBinError(ibin,ierror);
    }


    // Normalize only if required
    if ( normalize_ ) {
      Double_t integral1 = histo1_ptr->GetEntries();
      Double_t integral2 = histo2_ptr->GetEntries();
      if ( integral1 != 0 && integral2 != 0 ) {
        histo1_ptr->Scale(1./integral1);
        histo2_ptr->Scale(1./integral2);
      }
    }

    float Q    = ratioS2Beff->GetMaximum();
    int im     = ratioS2Beff->FindBin(Q);
    float errQ =  ratioS2Beff->GetBinError(im);
    float m    = im*(20./bin_);
    cout << "Q value: " << Q << " pm " << errQ << " @ m = " << m << endl;
    stringstream Qvalue;
    stringstream mvalue;
    Qvalue << "Q value: " << Q << " pm " << errQ;
    mvalue << "@ m: " << m;

    histo1_ptr->SetLineColor(kRed);
    stackVariable->Add(histo1_ptr);
    stackVariable->Add(histo2_ptr);
    legendVariable->SetFillColor(0);
    legendVariable->AddEntry(histo1_ptr, histo1_ptr->GetName(), "l");
    legendVariable->AddEntry(histo2_ptr, histo2_ptr->GetName(), "l");
    canvas->cd(1);
    stackVariable->Draw( OPTION );
    legendVariable->Draw();

    eff_histo1_ptr->SetLineColor(kRed);
    stackEfficiency->Add(eff_histo1_ptr);
    stackEfficiency->Add(eff_histo2_ptr);
    legendEfficiency->SetFillColor(0);
    legendEfficiency->AddEntry(eff_histo1_ptr, eff_histo1_ptr->GetName(), "pe");
    legendEfficiency->AddEntry(eff_histo2_ptr, eff_histo2_ptr->GetName(), "pe");
    canvas->cd(2);
    stackEfficiency->Draw( OPTION );
    legendEfficiency->Draw();

    legendSBeff->SetFillColor(0);
    legendSBeff->AddEntry(ratioSBeff,ratioSBeff->GetName(), "pe");
    canvas->cd(3);
    ratioSBeff->Draw();
    legendSBeff->Draw();

    legendS2Beff->SetFillColor(0);
    legendS2Beff->AddEntry(ratioS2Beff,ratioS2Beff->GetName(), "pe");
    legendS2Beff->AddEntry((Qvalue.str()).c_str(),(Qvalue.str()).c_str(), "");
    legendS2Beff->AddEntry((mvalue.str()).c_str(),(mvalue.str()).c_str(), "");
    canvas->cd(4);
    ratioS2Beff->Draw();
    legendS2Beff->Draw();

    // Write the canvas
    canvas->Write();

}

#endif //EFFICIENCYCANVAS_H
