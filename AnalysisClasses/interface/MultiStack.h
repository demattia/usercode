#ifndef MULTISTACK_H
#define MULTISTACK_H

#include <vector>
//#include <sstream>

#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"
#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"

using namespace std;

template<class T1, class T2>
class MultiStack {
 public:

  /// Gets the same arguments as the TCanvas constructor
  MultiStack( const TString& CANVASTITLE,
	      const TString& CANVASNAME,
	      const int      WIDTH      = 1000,
	      const int      HEIGHT     = 800 );

  ~MultiStack();

  /// Gets two MultiHisto objects pointers (can be MultiTH1F or MultiTProfile) 
  void Fill( const T1 * MULTIHISTO1,
	     const T1 * MULTIHISTO2 );

  void Write();

 private:
  TString canvasName_;
  TString canvasTitle_;
  int width_;
  int height_;
  vector<T2*> vecMultiHisto1_;
  vector<T2*> vecMultiHisto2_;
  TCanvas * canvas_;
};

// Methods
// -------

template<class T1, class T2>
MultiStack::MultiStack( const TString& CANVASNAME,
			const TString& CANVASTITLE,
			const int      WIDTH      = 1000,
			const int      HEIGHT     = 800 )
{
  canvasName_  = CANVASNAME;
  canvasTitle_ = CANVASTITLE;
  height_ = HEIGHT;
  width_ = WIDTH;
}

template<class T1, class T2>
MultiStack::~MultiStack()
{
}

template<class T1, class T2>
void MultiStack::Fill( const T1 * MULTIHISTO1,
		       const T1 * MULTIHISTO2 )
{
  vecMultiHisto1_  = MULTIHISTO1->multiHistos();
  vecMultiHisto2_  = MULTIHISTO2->multiHistos();
  TCanvas* canvas_ = new TCanvas( canvasName_, canvasTitle_, width_, height_ );
  canvas_->Divide(2,2);
}

template<class T1, class T2>
void MultiStack::Write()
{
  // Generate the name for the stack

  // Loop on the histograms
  int i = 0;
  vector<T2*>::const_iterator histo2_it = vecMultiHisto2_.begin();
  for ( vector<T2*>::const_iterator histo1_it = vecMultiHisto1_.begin();
	histo1_it != vecMultiHisto1_.end(); ++histo1_it, ++histo2_it, ++i ) {

    THStack * stack = new THStack( CANVASNAME, CANVASTITLE );
    TLegend * legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );

    // Normalize only if required

    Double_t integral1 = (*histo1_it)->GetEntries();
    Double_t integral2 = (*histo2_it)->GetEntries();
    if ( integral1 != 0 && integral2 != 0 ) {
      (*histo1_it)->Scale(1./integral1);
      (*histo2_it)->Scale(1./integral2);
    }
    (*histo_it)->SetLineColor(kRed);
    stack->Add(*histo1_it);
    stack->Add(*histo2_it);
    legend->AddEntry(*histo1_it, (*histo1_it)->GetName(), "l");
    legend->AddEntry(*histo2_it, (*histo2_it)->GetName(), "l");
    //  if(padCounter>1 && padCounter<6 ){
    canvas_->cd(i);
    stack->Draw("nostack");
    legend->Draw();
    // }
  } // end loop on histogram vectors
  canvas_->Write();

}

#endif //MULTISTACK_H
