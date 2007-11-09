#ifndef MULTITH1F_CC
#define MULTITH1F_CC

#include "AnalysisExamples/AnalysisClasses/interface/MultiTH1F.h"

MultiTH1F::MultiTH1F ( const char* NAME, const char* TITLE,
                       const int & BINVAL, const double & FIRSTVAL, const double & LASTVAL,
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
  HistoMean_ = new TH1F( Mean + Name, Mean + Title, BINPAR, FIRSTPAR, LASTPAR );

  // Multiple histograms
  double Multi = 0.;
  for ( int num = 0; num < BINPAR; ++num ) {
    Multi = FIRSTPAR + increment_*double(num);
    snum_ << Multi;
    vec_MultiHisto_.push_back( new TH1F( Name + snum_.str(), Title + snum_.str(), BINVAL, FIRSTVAL, LASTVAL ) );
    // Empty the ostringstream
    snum_.str("");
  }
  // Go back to the base dir
  outfile_->cd();
}

/// Fills the histogram corresponding to the index passed as second parameter with the value passed as first parameter
void MultiTH1F::Fill( const double & VAL, const int & PAR ) {
  // Check that the index is in the allowed range
  if ( PAR >= 0 && PAR < binpar_ ) {
    vec_MultiHisto_[PAR]->Fill( VAL );
  }
}

/// Writes the histograms to file
void MultiTH1F::Write() {

  // Put all the duplicate histograms in the subdir

  int sparse = binpar_/5;
  // Center labels for the mean histogram and put one on each bin
  HistoMean_->GetXaxis()->CenterLabels();
  HistoMean_->GetXaxis()->SetNdivisions(HistoMean_->GetSize()-2, false);
  HistoMean_->GetXaxis()->SetLabelSize(0.02);
  double Multi = 0;
  for( int num = 0; num < binpar_; ++num ) {
    Multi = firstpar_ + increment_*double(num);
    snum_ << Multi;

    TH1F* MultiHisto_ptr = vec_MultiHisto_[num];

    // Fill the mean histogram
    HistoMean_->SetBinContent( num+1, MultiHisto_ptr->GetMean() );
    HistoMean_->SetBinError( num+1, MultiHisto_ptr->GetMeanError() );

    // Fill the THStackLegend histogram
    StackLegend_->Add( MultiHisto_ptr, snum_.str().c_str(), true );

    // When the bin number is a multiple of one-tenth of the total number of bins
    if ( float(num+1)/float(sparse) == (num+1)/sparse ) {
      SparseStackLegend_->Add( MultiHisto_ptr, snum_.str().c_str(), true );
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

#endif // MULTITH1F_CC
