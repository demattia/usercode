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
#include "AnalysisExamples//AnalysisClasses/interface/ListFashionAttributedHisto.h"
//#include "TH1F.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"

#include <iostream>

template <class T>
class THStackLegend {
  
 public:
  THStackLegend( TString TITLE, 
		 double RIGHT  = 0.8, 
		 double BOTTOM = 0.68, 
		 double LEFT   = 0.98, 
		 double TOP    = 0.97
		 );

  ~THStackLegend();

  T* THStackLegend::Add(const T    * HISTO, 
			const char * LEGEND    = "", 
			const bool & NORMALIZE = false, 
			const char * FILL      = "l", 
			const bool & ERRORS    = false,
			const bool & DRESSED   = false,
			const int    COLOR     = 0,
			const int    MARKER    = 0
			);
  
  TAxis * GetXaxis() const;
  TAxis * GetYaxis() const;

  void Draw(TString OPTION = "");
  void Write(TString OPTION = "");

 private:
  THStack * stack_;
  TCanvas * canvas_;
  TLegend * legend_;
  int  counter_;
  bool drawn;
};

// methods
// -------
template<class T>
THStackLegend<T>::THStackLegend( TString TITLE, 
				 double RIGHT,
				 double BOTTOM,
				 double LEFT,
				 double TOP
				 ) {
  canvas_ = new TCanvas( TITLE + "_canvas", TITLE, 1000, 800 );
  stack_  = new THStack( TITLE + "_stack",  TITLE );
  legend_ = new TLegend( RIGHT,BOTTOM,LEFT,TOP );
  legend_->SetFillColor(0);
  legend_->SetBorderSize(0);  
  legend_->SetTextFont(72);
  legend_->SetTextSize(0.035);
  counter_ = 1;
  drawn    = false;
}

template<class T>
THStackLegend<T>::~THStackLegend() {
  delete stack_;
  delete canvas_;
}

template<class T>
T* THStackLegend<T>::Add(const T    * HISTO, 
			 const char * LEGEND,
			 const bool & NORMALIZE,
			 const char * FILL,
			 const bool & ERRORS,
			 const bool & DRESSED,
			 const int    COLOR,
			 const int    MARKER
			 ) {
  T * HISTO_ = (T*)HISTO->Clone();
  // Do not save this histogram in the current directory
  HISTO_->SetDirectory(0);
  
  if ( ERRORS ) {
    HISTO_->Sumw2();
    HISTO_->SetMarkerSize(0.7);
  }
  
  if (NORMALIZE) {
    Double_t integral = HISTO_->Integral();
    if ( integral != 0. ) {
      HISTO_->Scale(1./integral);
    }
    else {
      std::cout << "integral = " << integral << std::endl;
    }
  }
  stack_->Add(HISTO_);
  if ( FILL == "l" || FILL == "p" ) {
    if ( !DRESSED ) {
      if ( COLOR ) HISTO_->SetLineColor(COLOR);
      else HISTO_->SetLineColor(counter_);
      if ( MARKER ) {
	HISTO_->SetMarkerStyle(MARKER);
	if ( COLOR ) HISTO_->SetMarkerColor(COLOR);
	else HISTO_->SetMarkerColor(counter_);
      } else {
	HISTO_->SetMarkerStyle(counter_+20);
	if (COLOR) HISTO_->SetMarkerColor(COLOR);
	else HISTO_->SetMarkerColor(counter_);
      }
    }
    HISTO_->SetLineWidth(2);
  }
  else if ( FILL == "f" ) {
    if ( !DRESSED ) {
      if ( COLOR ) HISTO_->SetFillColor(COLOR);
      else HISTO_->SetFillColor(counter_+1);
    }
  }
  if ( LEGEND == "" ) {
    legend_->AddEntry(HISTO_, HISTO_->GetName(), FILL);
  }
  else {
    legend_->AddEntry(HISTO_, LEGEND, FILL);
  }
  ++counter_;
  
  return HISTO_;
}

template<class T>
TAxis * THStackLegend<T>::GetXaxis() const {
  return ( stack_->GetXaxis() );
}

template<class T>
TAxis * THStackLegend<T>::GetYaxis() const {
  return ( stack_->GetHistogram()->GetYaxis() );
}

template<class T>
void THStackLegend<T>::Draw(TString OPTION) {
  if ( !drawn ) {
    canvas_->cd();
    stack_->Draw(OPTION);
    legend_->Draw();
  }
  drawn = true;
}

template<class T>
void THStackLegend<T>::Write(TString OPTION) {
  if ( !drawn ) {
    canvas_->cd();
    stack_->Draw(OPTION);
    legend_->Draw();
  }
  canvas_->Write();
  drawn = true;
}

#endif
