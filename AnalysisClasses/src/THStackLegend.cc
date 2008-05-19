

#ifndef THSTACKLEGEND_CC
#define THSTACKLEGEND_CC

#include "AnalysisExamples/AnalysisClasses/interface/THStackLegend.h"

TH1* THStackLegend::Add(const TH1  * HISTO, 
			const char * LEGEND, 
			const bool & NORMALIZE, 
			const char * FILL, 
			const bool & ERRORS,
			const int    COLOR,
			const int    MARKER
			) {
  TH1 * HISTO_ = (TH1*)HISTO->Clone();
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
    if ( COLOR ) HISTO_->SetLineColor(COLOR);
    else HISTO_->SetLineColor(counter_);
    if ( MARKER ) {
      HISTO_->SetMarkerStyle(MARKER);
      if (COLOR) HISTO_->SetMarkerColor(COLOR);
      else HISTO_->SetMarkerColor(counter_);
    } else {
      HISTO_->SetMarkerStyle(counter_+20);
      if (COLOR) HISTO_->SetMarkerColor(COLOR);
      else HISTO_->SetMarkerColor(counter_);
    }
    HISTO_->SetLineWidth(2);
  }
  else if ( FILL == "f" ) {
    if (COLOR) HISTO_->SetFillColor(COLOR);
    else HISTO_->SetFillColor(counter_+1);
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

TProfile* THStackLegend::Add(const TProfile * HISTO, const char* LEGEND, const bool & NORMALIZE, const char * FILL) {
  TProfile * HISTO_ = new TProfile(*HISTO);
  // Do not save this histogram in the current directory
  HISTO_->SetDirectory(0);

  if (NORMALIZE) {
    Double_t integral = HISTO_->Integral();
    if ( integral != 0 ) {
      HISTO_->Scale(1./integral);
    }
    else {
      std::cout << "integral = " << integral << std::endl;
    }
  }
  stack_->Add(HISTO_);
  if ( FILL == "l" || FILL == "p" ) {
    HISTO_->SetLineColor(counter_);
    HISTO_->SetMarkerStyle(counter_+20);
    HISTO_->SetMarkerColor(counter_);
    HISTO_->SetLineWidth(2);
  }
  else if ( FILL == "f" ) {
    //   if (counter_ == 4) counter_=counter_+2;
    HISTO_->SetFillColor(counter_+1);
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

void THStackLegend::Draw(TString OPTION) {
  if ( !drawn ) {
    canvas_->cd();
    stack_->Draw(OPTION);
    legend_->Draw();
  }
  drawn = true;
}

void THStackLegend::Write(TString OPTION) {
  if ( !drawn ) {
    canvas_->cd();
    stack_->Draw(OPTION);
    legend_->Draw();
  }
  canvas_->Write();
  drawn = true;
}

#endif
