#ifndef THSTACKLEGEND_HH
#define THSTACKLEGEND_HH

/**
 * Simple THStack container that automatically generates a legend and modifies the
 * color of the different histograms to make them recognizable.
 * The THStack works with any TH1 inheriting histograms, so it works on TH1F and
 * TProfile as well.
 */

#include "THStack.h"
#include "TH1.h"
//#include "TH1F.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"

#include <iostream>

class THStackLegend {
 public:
  THStackLegend( TString TITLE, 
		 double RIGHT  = 0.55, 
		 double BOTTOM = 0.82, 
		 double LEFT   = 0.76, 
		 double TOP    = 0.65
		 ) {
    canvas_ = new TCanvas( TITLE + "_canvas", TITLE, 1000, 800 );
    stack_  = new THStack( TITLE + "_stack",  TITLE );
    legend_ = new TLegend( RIGHT,BOTTOM,LEFT,TOP );
    counter_ = 1;
    drawn    = false;
  }
  ~THStackLegend() {
    delete stack_;
    delete canvas_;
  }

  TH1* Add(const TH1  * HISTO, 
	   const char * LEGEND    = "", 
	   const bool & NORMALIZE = false, 
	   const char * FILL      = "l", 
	   const bool & ERRORS    = false,
	   const int    COLOR     = 0,
	   const int    MARKER    = 0
	   );
  TProfile* Add(const TProfile * HISTO, const char* LEGEND = "", const bool & NORMALIZE = false, const char* FILL = "l" );

  TAxis * GetXaxis() const {
    return ( stack_->GetXaxis() );
  }

  TAxis * GetYaxis() const {
    return ( stack_->GetHistogram()->GetYaxis() );
  }

  void Draw(TString OPTION = "");
  void Write(TString OPTION = "");

 private:
  THStack * stack_;
  TCanvas * canvas_;
  TLegend * legend_;
  int  counter_;
  bool drawn;
};

#endif
