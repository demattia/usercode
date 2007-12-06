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

TH1F* THStackLegend::Add(const TH1F * HISTO, const char* LEGEND, const bool & NORMALIZE, const char * FILL) {
  TH1F * HISTO_ = new TH1F(*HISTO);
  // Do not save this histogram in the current directory
  HISTO_->SetDirectory(0);
  //  TH1F * HISTO_ = (TH1F*)HISTO->Clone();

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

//   TH1* histo_ptr = HISTO_;
//   TH1F * returnHisto_ptr = dynamic_cast<TH1F*>( Add( histo_ptr, LEGEND, NORMALIZE, FILL ) );
//   return returnHisto_ptr;
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

//   TH1* histo_ptr = HISTO_;
//   TProfile * returnHisto_ptr = dynamic_cast<TProfile*>( Add( histo_ptr, LEGEND, NORMALIZE, FILL ) );
//   return returnHisto_ptr;

}

// /** This makes an internal copy of the histograms.
//  * This way we avoid multiple normalizations.
//  * The histogram is created with new, but it is never deleted.
//  * It should be managed by root, which is the case if it is put in a file.
//  */
// TH1* THStackLegend::Add(const TH1 * HISTO, const char* LEGEND, const bool & NORMALIZE, const char * FILL) {

// }

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
