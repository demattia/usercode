#ifndef THSTACKLEGEND_CC
#define THSTACKLEGEND_CC

#include "AnalysisExamples/AnalysisClasses/interface/THStackLegend.h"

// void THStackLegend::Add(TH1F * HISTO, char* LEGEND, char* FILL) {
//   stack_->Add(HISTO);
//   if ( FILL == "l") {
//     HISTO->SetLineColor(counter_);
//   }
//   else if ( FILL == "f" ) {
// //    if (counter_ == 4) counter_=counter_+2;
//     HISTO->SetFillColor(counter_+1);
//   }
//   if ( LEGEND == "" ) {
//     legend_->AddEntry(HISTO, HISTO->GetName(), FILL);
//   }
//   else {
//     legend_->AddEntry(HISTO, LEGEND, FILL);
//   }
//   ++counter_;
// }

// This makes an internal copy of the histograms
// this way we avoid multiple normalizations
TH1F* THStackLegend::Add(const TH1F * HISTO, const char* LEGEND, const bool & NORMALIZE, const char * FILL) {
  TH1F * HISTO_ = new TH1F(*HISTO);
  //  TH1F * HISTO_ = (TH1F*)HISTO->Clone();

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

void THStackLegend::Draw(TString OPTION) const {
  canvas_->cd();
  stack_->Draw(OPTION);
  legend_->Draw();
}

void THStackLegend::Write(TString OPTION) const {
  canvas_->cd();
  stack_->Draw(OPTION);
  legend_->Draw();
  canvas_->Write();
}

#endif
