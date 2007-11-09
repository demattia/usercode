#ifndef OFFLINEJET_H
#define OFFLINEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseAll.h"

#include <cmath>
#include <vector>

/**
 *
 * Used for offline jets. Includes b-tagging informations.
 * Inherits from BaseJet.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class OfflineJet : public BaseJet {
  public:
    OfflineJet( const double & ET, const double & ETA, const double & PHI, const double & UNCORRET, const double & EMENERGYFRACTION ) : BaseJet( ET, ETA, PHI ) {
      uncorrEt_ = UNCORRET;
      emEnergyFraction_ = EMENERGYFRACTION;
    }
    // Default constructor, only needed for classes.h
    OfflineJet() : BaseJet( 0., 0., 0. ) {
      uncorrEt_ = 0.;
      emEnergyFraction_ = 0.;
    }
    double uncorrEt() const { return uncorrEt_; }
    double emEnergyFraction() const { return emEnergyFraction_; }
    void setUncorrEt( const double & UNCORRET ) { uncorrEt_ = UNCORRET; }
    void setEmEnergyFraction( const double & EMENERGYFRACTION ) { emEnergyFraction_ = EMENERGYFRACTION; }
  protected:
    double uncorrEt_;
    double emEnergyFraction_;
  };

  typedef std::vector<OfflineJet> OfflineJetCollection;

}

#endif // OFFLINEJET_H
