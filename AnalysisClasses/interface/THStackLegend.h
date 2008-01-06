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
  THStackLegend( TString TITLE ) {
    canvas_ = new TCanvas( TITLE + "_canvas", TITLE, 1000, 800 );
    stack_ = new THStack( TITLE + "_stack", TITLE );
    legend_= new TLegend( 0.55, 0.65, 0.76, 0.82 );
    counter_ = 1;
    drawn = false;
  }
  ~THStackLegend() {
    delete stack_;
    delete canvas_;
  }
//  void Add(TH1F * HISTO, char* LEGEND = "", char* FILL = "l");
  TH1* Add(const TH1 * HISTO, const char* LEGEND = "", const bool & NORMALIZE = false, const char* FILL = "l");
  TProfile* Add(const TProfile * HISTO, const char* LEGEND = "", const bool & NORMALIZE = false, const char* FILL = "l");

  TAxis * GetXaxis() const {
//    return ( stack_->GetHistogram()->GetXaxis() );
    return ( stack_->GetXaxis() );
  }

  TAxis * GetYaxis() const {
    return ( stack_->GetHistogram()->GetYaxis() );
  }

  void Draw(TString OPTION = "");
  void Write(TString OPTION = "") const;
 private:
//  TH1* Add(const TH1 * HISTO, const char* LEGEND = "", const bool & NORMALIZE = false, const char* FILL = "l");
  THStack * stack_;
  TCanvas * canvas_;
  TLegend * legend_;
  int counter_;
  bool drawn;
};

#endif
