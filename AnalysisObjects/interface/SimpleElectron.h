#ifndef SIMPLEELECTRON_H
#define SIMPLEELECTRON_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"

#include <cmath>
#include <vector>

/**
 *
 * Stores an offline electron. Inlcudes Et, had/em fraction and
 * isolation.
 *
 * Inherits from BaseParticle.
 *
 * Author M. De Mattia - 8/11/2007
 *
 */

namespace anaobj {

  class SimpleElectron : public BaseParticle {
  public:
    SimpleElectron( const double                  & PT, 
		    const double                  & ETA, 
		    const double                  & PHI,
		    const int                       CHARGE,
		    const math::XYZTLorentzVector & P4,    
		    const math::XYZPoint          & VERTEX,    
                    const double                  & ET, 
		    const double                  & HADOVEREM, 
		    const double                  & ISOVAL,
		    const double                  & IMPACTPARAMXY,
		    const double                  & ERRORIMPACTPARAMXY) : BaseParticle( PT, ETA, PHI ) {
      charge_             = CHARGE;
      p4_                 = P4;
      vertex_             = VERTEX;
      et_                 = ET;
      hadOverEm_          = HADOVEREM;
      isoVal_             = ISOVAL;
      impactParamXY_      = IMPACTPARAMXY;
      errorImpactParamXY_ = ERRORIMPACTPARAMXY;
    }
    /// Default constructor, only needed for classes.h
    SimpleElectron() : BaseParticle( 0., 0., 0. ) {
      charge_             = 0;
      p4_                 = math::XYZTLorentzVector(0.,0.,0.,0.);
      vertex_             = math::XYZPoint(0.,0.,0.);
      et_                 = 0.;
      hadOverEm_          = 0.;
      isoVal_             = 0.;
      impactParamXY_      = 0.;
      errorImpactParamXY_ = 0.;
    }
    int           charge()             const { return charge_;             }
    math::XYZTLorentzVector p4()       const { return p4_;                 }
    math::XYZPoint vertex()            const { return vertex_;             }
    double        et()                 const { return et_;                 }
    double        hadOverEm()          const { return hadOverEm_;          }
    double        isoVal()             const { return isoVal_;             }
    double        impactParamXY()      const { return impactParamXY_;      }
    double        errorImpactParamXY() const { return errorImpactParamXY_; }

    void setCharge(             const int                       CHARGE             ) { charge_             = CHARGE;             }
    void setP4(                 const math::XYZTLorentzVector & P4                 ) { p4_                 = P4;                 }
    void setVertex(             const math::XYZPoint          & VERTEX             ) { vertex_             = VERTEX;             }
    void setEt(                 const double                  & ET                 ) { et_                 = ET;                 }
    void setHadOverEm(          const double                  & HADOVEREM          ) { hadOverEm_          = HADOVEREM;          }
    void setIsoVal(             const double                  & ISOVAL             ) { isoVal_             = ISOVAL;             }
    void setImpactParamXY(      const double                  & IMPACTPARAMXY      ) { impactParamXY_      = IMPACTPARAMXY;      }
    void setErrorImpactParamXY( const double                  & ERRORIMPACTPARAMXY ) { errorImpactParamXY_ = ERRORIMPACTPARAMXY; }
  protected:
    int                     charge_;
    math::XYZTLorentzVector p4_;
    math::XYZPoint          vertex_;
    double                  et_;
    double                  hadOverEm_;
    double                  isoVal_;
    double                  impactParamXY_;
    double                  errorImpactParamXY_;
  };
  
  typedef std::vector<SimpleElectron> SimpleElectronCollection;
  
}

#endif // GLOBALMUON_H
