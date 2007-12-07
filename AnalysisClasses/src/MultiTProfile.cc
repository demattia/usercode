
#ifndef MULTITPROFILE_CC
#define MULTITPROFILE_CC

#include "AnalysisExamples/AnalysisClasses/interface/MultiTProfile.h"

MultiTProfile::MultiTProfile ( const char* NAME, const char* TITLE,
                       const int & BINVALX, const double & FIRSTVALX, const double & LASTVALX,
                       const double & FIRSTVALY, const double & LASTVALY,
                       const int & BINPAR, const double & FIRSTPAR, const double & LASTPAR,
                       TFile* OUTFILE ) {

  // White background for the canvases
  gROOT->SetStyle("Plain");

  binpar_ = BINPAR;
  firstpar_ = FIRSTPAR;
  lastpar_ = LASTPAR;
  increment_ = ( lastpar_ - firstpar_ )/double(binpar_);
  TString Name( NAME );
  TString Title( TITLE );
  TString Mean( "Mean_" );
  TString Stack( "Stack_" );
  TString SparseStack( "SparseStack_" );
  outfile_ = OUTFILE;

  Directory_ = outfile_->mkdir( Name );

  Directory_->cd();

  // THStackLegend histogram
  StackLegend_ = new THStackLegend( Stack + Name );
  // Sparse THStackLegend histogram
  SparseStackLegend_ = new THStackLegend( SparseStack + Name );

  // Mean histogram
  //   HistoMean_ = new TH1F( Mean + Name, Mean + Title, BINPAR, FIRSTPAR, LASTPAR );

  // Multiple histograms
  double Multi = 0.;
  for ( int num = 0; num < BINPAR; ++num ) {
    Multi = FIRSTPAR + increment_*double(num);
    snum_ << Multi;
    vec_MultiHisto_.push_back( new TProfile( Name + "_" + snum_.str(), Title + "_" + snum_.str(), BINVALX, FIRSTVALX, LASTVALX, FIRSTVALY, LASTVALY ) );
    // Empty the ostringstream
    snum_.str("");
  }
  // Go back to the base dir
  outfile_->cd();
}

/// Fills the histogram corresponding to the index passed as second parameter with the value passed as first parameter
void MultiTProfile::Fill( const double & VALX, const double & VALY, const int & PAR ) {
  // Check that the index is in the allowed range
  if ( PAR >= 0 && PAR < binpar_ ) {
    vec_MultiHisto_[PAR]->Fill( VALX, VALY );
  }
}

/// Writes the histograms to file
void MultiTProfile::Write() {

  // Put all the duplicate histograms in the subdir

  int sparse = binpar_/5;
  // Center labels for the mean histogram and put one on each bin
  //   HistoMean_->GetXaxis()->CenterLabels();
  //   HistoMean_->GetXaxis()->SetNdivisions(HistoMean_->GetSize()-2, false);
  //   HistoMean_->GetXaxis()->SetLabelSize(0.02);
  double Multi = 0;
  for( int num = 0; num < binpar_; ++num ) {
    Multi = firstpar_ + increment_*double(num);
    snum_ << Multi;

    TProfile* MultiHisto_ptr = vec_MultiHisto_[num];

    // Fill the mean histogram
    //     HistoMean_->SetBinContent( num+1, MultiHisto_ptr->GetMean() );
    //     HistoMean_->SetBinError( num+1, MultiHisto_ptr->GetMeanError() );

    // Fill the THStackLegend histogram
    StackLegend_->Add( MultiHisto_ptr, snum_.str().c_str(), false );

    // When the bin number is a multiple of one-tenth of the total number of bins
    if ( float(num+1)/float(sparse) == (num+1)/sparse ) {
      SparseStackLegend_->Add( MultiHisto_ptr, snum_.str().c_str(), false );
    }

    // Empty the ostringstream
    snum_.str("");
  }

  // Put the stack histograms in the correct directory
  Directory_->cd();

  StackLegend_->Write("nostack");
  SparseStackLegend_->Write("nostack");
  // Go back to the base dir
  outfile_->cd();
}

std::vector<TProfile*> MultiTProfile::multiProfiles() const {
  
  return vec_MultiHisto_;
}

#endif // MULTITPROFILE_CC
