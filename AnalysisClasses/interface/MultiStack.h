#ifndef MULTISTACK_H
#define MULTISTACK_H

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
class MultiStack {
 public:

  /**
   * Gets the same arguments as the TCanvas constructor, plus extra arguments to specify the number of histograms per
   * canvas and indeces to delimit a subgroup of histograms. <br>
   * The last paramenter determines if the histograms are to be normalized. Deafault is false.
   */
  MultiStack( const TString& CANVASTITLE,
	      const TString& CANVASNAME,
	      const int      WIDTH          = 1000,
	      const int      HEIGHT         = 800,
              const int      HISTOPERCANVAS = 4,
              const int      FIRSTHISTO     = 1,
              const int      LASTHISTO      = 0,
              const bool     NORMALIZE      = false );

  ~MultiStack();

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
  bool normalize_;
  vector<T2*> vecMultiHisto1_;
  vector<T2*> vecMultiHisto2_;
  vector<TCanvas*> vecCanvasPtr_;
  int histoPerCanvas_;
  int firstHisto_;
  int lastHisto_;
  int canvasNum_;
};

// Methods
// -------

template<class T1, class T2>
MultiStack<T1, T2>::MultiStack( const TString& CANVASNAME,
                                const TString& CANVASTITLE,
                                const int      WIDTH,
                                const int      HEIGHT,
                                const int      HISTOPERCANVAS,
                                const int      FIRSTHISTO,
                                const int      LASTHISTO,
                                const bool     NORMALIZE)
{
  canvasName_     = CANVASNAME;
  canvasTitle_    = CANVASTITLE;
  height_         = HEIGHT;
  width_          = WIDTH;
  normalize_      = NORMALIZE;
  histoPerCanvas_ = HISTOPERCANVAS;
  firstHisto_     = FIRSTHISTO-1;
  lastHisto_      = LASTHISTO-1;
  canvasNum_      = 0;
}

template<class T1, class T2>
MultiStack<T1, T2>::~MultiStack()
{
}

template<class T1, class T2>
void MultiStack<T1, T2>::Fill( const T1 * MULTIHISTO1,
                               const T1 * MULTIHISTO2 )
{
  vecMultiHisto1_ = MULTIHISTO1->multiHistos();
  vecMultiHisto2_ = MULTIHISTO2->multiHistos();
}

template<class T1, class T2>
void MultiStack<T1, T2>::Write( const TString & OPTION )
{
  // If lastHisto_ == 0 loop on all histograms
  if ( lastHisto_ == 0 ) vecMultiHisto1_.size();

  // Create the new canvases
  canvasNum_ = (lastHisto_-firstHisto_)/histoPerCanvas_;
  // If the number of histogram is less than the number of histogram per canvas, generate one canvas anyway.
  if ( canvasNum_ == 0 ) canvasNum_ = 1;
  // If the residual of the integer division is != 0 create an additional canvas.
  if ( (lastHisto_-firstHisto_)%histoPerCanvas_ != 0 ) ++canvasNum_;  

  for ( int iCanvas=0; iCanvas<canvasNum_; ++iCanvas ) {
    stringstream canvasIndex;
    canvasIndex << iCanvas;
    TCanvas * canvasPtr = new TCanvas( canvasName_ + "_" + canvasIndex.str(), canvasTitle_ + "_" + canvasIndex.str(), width_, height_ );
    canvasPtr->Divide(2, int((histoPerCanvas_+1)/2) );
    vecCanvasPtr_.push_back( canvasPtr );
  }

  // Loop on the histograms
  //  typename vector<T2*>::const_iterator histo1_it = vecMultiHisto1_.begin();
  //  typename vector<T2*>::const_iterator histo2_it = vecMultiHisto2_.begin();

  // Canvas index
  int canvas = 0;

  // Loop only on selected histograms
  int histoIndex=firstHisto_;
  for ( int i=0; i<(lastHisto_-firstHisto_); ++i, ++histoIndex ) {
    //, ++histo1_it, ++histo2_it ) {

    T2* histo1_ptr = vecMultiHisto1_[histoIndex];
    T2* histo2_ptr = vecMultiHisto2_[histoIndex];

    // Generate the name for the stack
    string histoName( histo1_ptr->GetName() );
    // Take the index of the last character in the string matching the underscore.
    // In MultiHisto the index is added as _num. This gives us the part we need to
    // create the index to append to the stack name and title.
    TString histoIndex = histoName.substr( histoName.find_last_of('_') );

    THStack * stack = new THStack( canvasName_ + histoIndex, canvasTitle_ + histoIndex );
    TLegend * legend = new TLegend( 0.55, 0.65, 0.76, 0.82 );

    // Normalize only if required
    if ( normalize_ ) {
      Double_t integral1 = histo1_ptr->GetEntries();
      Double_t integral2 = histo2_ptr->GetEntries();
      if ( integral1 != 0 && integral2 != 0 ) {
        histo1_ptr->Scale(1./integral1);
        histo2_ptr->Scale(1./integral2);
      }
    }

    histo1_ptr->SetLineColor(kRed);
    stack->Add(histo1_ptr);
    stack->Add(histo2_ptr);
    legend->AddEntry(histo1_ptr, histo1_ptr->GetName(), "l");
    legend->AddEntry(histo2_ptr, histo2_ptr->GetName(), "l");
    vecCanvasPtr_[canvas]->cd(i+1-(canvas*histoPerCanvas_));
    stack->Draw( OPTION );
    legend->Draw();
    // Increase the canvas index if needed
    if ( i==(histoPerCanvas_-1) ) {
      ++canvas;
    }

  } // end loop on histogram vectors

  // Write all the canvases
  for ( vector<TCanvas*>::const_iterator canvas_it = vecCanvasPtr_.begin();
        canvas_it != vecCanvasPtr_.end(); ++canvas_it ) {
    (*canvas_it)->Write();
  }
}

#endif //MULTISTACK_H
