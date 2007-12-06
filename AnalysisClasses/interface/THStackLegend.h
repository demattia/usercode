#ifndef THSTACKLEGEND_HH
#define THSTACKLEGEND_HH

#include "THStack.h"
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
  }
  ~THStackLegend() {
    delete stack_;
    delete canvas_;
  }
//  void Add(TH1F * HISTO, char* LEGEND = "", char* FILL = "l");
  TH1F* Add(const TH1F * HISTO, const char* LEGEND = "", const bool & NORMALIZE = false, const char* FILL = "l");
  void Draw(TString OPTION = "") const;
  void Write(TString OPTION = "") const;
 private:
  THStack * stack_;
  TCanvas * canvas_;
  TLegend * legend_;
  int counter_;
};

#endif
