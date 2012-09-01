#ifndef PileupReweighting_h
#define PileupReweighting_h

#include "TH1.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <TCanvas.h>
#include <iostream>


class PileupReweighting {
  public:
    PileupReweighting() { };
    PileupReweighting( TString dataFileName );
    ~PileupReweighting();
    double weight( float npv ) ;

  private:
    TString dataFileName_;
    TFile * dataFile_;

    TH1F * MC_distr_; // Normalised MC distribution
    TH1F * Data_distr_; // Normalised Data distribution

    TH1F * weights_; // Weights

};

#endif
