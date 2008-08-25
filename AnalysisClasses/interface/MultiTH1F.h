#ifndef MULTITH1F_H
#define MULTITH1F_H

#include <vector>
#include <sstream>

#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"

#include "AnalysisExamples/AnalysisClasses/interface/THStackLegend.h"

/**
 * This class can be used to manage loops in which to produce
 * many histograms for the same quantities while varying some parameter.
 * It has the same methods of the TH1F but it puts everything inside a
 * directory with the name passed as histogram name.
 *
 * The constructor is called like a TH2F but requiring also a TFile*:
 * MultiTH1F multi( "name", "title", bins val, first val, last val, bin par, first par, last par, TFile* outfile ).
 *
 * Where the "val" parameters refer to the variable into consideretion ( for
 * example eta, or pt).
 * The "par" parameters refer to the parameter which is being varied (for
 * example dz, when considering distribution of the pt of the primary vertex
 * as a function of the minimum dz used to reconstruct it).
 * "bin val" is used as the size of the internal vector of TH1F.
 * The "first par" and "last par" values are used as the first and last values for the mean histogram.
 * The "first val" and "last val" values are used as the first and last values for the multiple histograms.
 * It also requires a TFile*, for the file where to create the directory and save the histograms.
 *
 * ATTENTION
 * If using a variable not limited (or semilimited), the "first val" and "last val" values
 * should be given big, because underflow and overflow bins are not considered in mean calculation,
 * and the mean histogram would result in wrong values.
 *
 * The Write method evaluates the mean histogram, which contains mean values from all the multiple
 * histograms. It has the name = Mean_"name".
 * It also creates a THStackLegend with all the multiple histograms superimposed, rescaled and not
 * stacked. It has the name = Stack_"name".
 *
 * Author M. De Mattia
 *
 * Updated to use the templated THStackLegend - 28/5/2008
 */

class MultiTH1F {
 public:

  MultiTH1F ( const char* NAME, const char* TITLE,
              const int & BINVAL, const double & FIRSTVAL, const double & LASTVAL,
              const int & BINPAR, const double & FIRSTPAR, const double & LASTPAR,
              TFile* OUTFILE );

  /// Fills the histogram corresponding to the index passed as second parameter 
  /// with the value passed as first parameter
  void Fill( const double & VAL, const int & PAR );

  /// Writes the histograms to file
  void Write();

  /// Returns vector of internal histograms
  std::vector<TH1F*> multiHistos() const;

 private:
  TDirectory * Directory_;
  TH1F * HistoMean_;
  THStackLegend<TH1F> * StackLegend_;
  THStackLegend<TH1F> * SparseStackLegend_;
  std::vector<TH1F*> vec_MultiHisto_;
  std::ostringstream snum_;
  double firstpar_;
  double lastpar_;
  int binpar_;
  double increment_;
  TFile* outfile_;
};

#endif //MULTITH1F_H
