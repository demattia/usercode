#ifndef GLOBALMUON_H
#define GLOBALMUON_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

#include <cmath>
#include <vector>

/**
 *
 * Stores a ParamGlobalMuon. Includes em, had and ho energy and
 * sumPt from tracks for isolation.
 *
 * Inherits from BaseParticle.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class GlobalMuon : public BaseParticle {
  public:
    GlobalMuon( const double & PT, const double & ETA, const double & PHI,
                const double & CALOEM, const double & CALOHAD, const double & CALOHO,
                const double & ISO03SUMPT, const double & NORMCHI2, const int HITSNUM ) : BaseParticle( PT, ETA, PHI ) {
      caloEm_ = CALOEM;
      caloHad_ = CALOHAD;
      caloHo_ = CALOHO;
      iso03SumPt_ = ISO03SUMPT;
      normChi2_ = NORMCHI2;
      hitsNum_ = HITSNUM;
    }
    /// Default constructor, only needed for classes.h
    GlobalMuon() : BaseParticle( 0., 0., 0. ) {
      caloEm_ = 0.;
      caloHad_ = 0.;
      caloHo_ = 0.;
      iso03SumPt_ = 0.;
      normChi2_ = 0.;
      hitsNum_ = 0;
    }
    double caloEm() const { return caloEm_; }
    double caloHad() const { return caloHad_; }
    double caloHo() const { return caloHo_; }
    double iso03SumPt() const { return iso03SumPt_; }
    double normChi2() const { return normChi2_; }
    int hitsNum() const { return hitsNum_; }
    void setCaloEm( const double & CALOEM ) { caloEm_ = CALOEM; }
    void setCaloHad( const double & CALOHAD ) { caloHad_ = CALOHAD; }
    void setCaloHo( const double & CALOHO ) { caloHo_ = CALOHO; }
    void setIso03SumPt( const double & ISO03SUMPT ) { iso03SumPt_ = ISO03SUMPT; }
    void setNormChi2( const double & NORMCHI2 ) { normChi2_ = NORMCHI2; }
    void setHitsNum( const int HITSNUM ) { hitsNum_ = HITSNUM; }
  protected:
    double caloEm_;
    double caloHad_;
    double caloHo_;
    double iso03SumPt_;
    double normChi2_;
    int hitsNum_;
  };

  typedef std::vector<GlobalMuon> GlobalMuonCollection;

}

#endif // GLOBALMUON_H
