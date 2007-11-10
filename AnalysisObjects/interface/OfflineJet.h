#ifndef OFFLINEJET_H
#define OFFLINEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"

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
    OfflineJet( const double & ET, const double & ETA, const double & PHI, const double & UNCORRET, const double & EMENERGYFRACTION,
                const float & DISCRIMINATORHIGHEFF, const float & DISCRIMINATORHIGHPUR, const int TKNUM,
                const double & TKSUMPT, const double & JETMASS ) : BaseJet( ET, ETA, PHI ) {
      uncorrEt_ = UNCORRET;
      emEnergyFraction_ = EMENERGYFRACTION;
      discriminatorHighEff_ = DISCRIMINATORHIGHEFF;
      discriminatorHighPur_ = DISCRIMINATORHIGHPUR;
      tkNum_ = TKNUM;
      tkSumPt_ = TKSUMPT;
      jetMass_ = JETMASS;
    }
    // Default constructor, only needed for classes.h
    OfflineJet() : BaseJet( 0., 0., 0. ) {
      uncorrEt_ = 0.;
      emEnergyFraction_ = 0.;
      discriminatorHighEff_ = 0.;
      discriminatorHighPur_ = 0.;
      tkNum_ = 0;
      tkSumPt_ = 0.;
      jetMass_ = 0.;
    }
    double uncorrEt() const { return uncorrEt_; }
    double emEnergyFraction() const { return emEnergyFraction_; }
    float discriminatorHighEff() const { return discriminatorHighEff_; }
    float discriminatorHighPur() const { return discriminatorHighPur_; }
    int tkNum() const { return tkNum_; }
    double tkSumPt() const { return tkSumPt_; }
    double jetMass() const { return jetMass_; }
    void setUncorrEt( const double & UNCORRET ) { uncorrEt_ = UNCORRET; }
    void setEmEnergyFraction( const double & EMENERGYFRACTION ) { emEnergyFraction_ = EMENERGYFRACTION; }
    void setDiscriminatorHighEff( const float & DISCRIMINATORHIGHEFF ) { discriminatorHighEff_ = DISCRIMINATORHIGHEFF; }
    void setDiscriminatorHighPur( const float & DISCRIMINATORHIGHPUR ) { discriminatorHighPur_ = DISCRIMINATORHIGHPUR; }
    void setTkNum( const int TKNUM ) { tkNum_ = TKNUM; }
    void setTkSumPt( const double & TKSUMPT ) { tkSumPt_ = TKSUMPT; }
    void setJetMass( const double & JETMASS ) { jetMass_ = JETMASS; }
  protected:
    double uncorrEt_;
    double emEnergyFraction_;
    float discriminatorHighEff_;
    float discriminatorHighPur_;
    int tkNum_;
    double tkSumPt_;
    double jetMass_;
  };

  typedef std::vector<OfflineJet> OfflineJetCollection;

}

#endif // OFFLINEJET_H
