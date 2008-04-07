#ifndef GLOBALMUON_H
#define GLOBALMUON_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Candidate/interface/Particle.h"
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
    GlobalMuon( const double                  & PT, 
		const double                  & ETA, 
		const double                  & PHI,
		const int                       CHARGE,
		const math::XYZTLorentzVector & P4,    
		const math::XYZPoint          & VERTEX,
                const double                  & CALOEM, 
		const double                  & CALOHAD, 
		const double                  & CALOHO,
                const float                   & ISO03SUMPT, 
		const float                   & ISO03EMET, 
		const float                   & ISO03HADET, 
		const float                   & ISO03HOET, 
		const int                     & ISO03NJETS, 
		const int                     & ISO03NTRACKS, 
		const double                  & NORMCHI2, 
		const int                       HITSNUM,
		const double                  & IMPACTPARAMXY, 
		const double                  & ERRORIMPACTPARAMXY ) : BaseParticle( PT, ETA, PHI ) {
      charge_             = CHARGE;
      p4_                 = P4;
      vertex_             = VERTEX;
      caloEm_             = CALOEM;
      caloHad_            = CALOHAD;
      caloHo_             = CALOHO;
      iso03SumPt_         = ISO03SUMPT;
      iso03emEt_          = ISO03EMET;
      iso03hadEt_         = ISO03HADET;
      iso03hoEt_          = ISO03HOET;
      iso03nJets_         = ISO03NJETS;
      iso03nTracks_       = ISO03NTRACKS;
      normChi2_           = NORMCHI2;
      hitsNum_            = HITSNUM;
      impactParamXY_      = IMPACTPARAMXY;
      errorImpactParamXY_ = ERRORIMPACTPARAMXY;
    }
    /// Default constructor, only needed for classes.h
    GlobalMuon() : BaseParticle( 0., 0., 0. ) {
      charge_             = 0;
      p4_                 = math::XYZTLorentzVector(0.,0.,0.,0.);
      vertex_             = math::XYZPoint(0.,0.,0.);
      caloEm_             = 0.;
      caloHad_            = 0.;
      caloHo_             = 0.;
      iso03SumPt_         = 0.;
      iso03emEt_          = 0.;
      iso03hadEt_         = 0.;
      iso03hoEt_          = 0.;
      iso03nJets_         = 0;
      iso03nTracks_       = 0;
      normChi2_           = 0.;
      hitsNum_            = 0;
      impactParamXY_      = 0.;
      errorImpactParamXY_ = 0.;
    }
    int           charge()             const { return charge_;             }
    math::XYZTLorentzVector p4()       const { return p4_;                 }
    math::XYZPoint vertex()            const { return vertex_;             }
    double        caloEm()             const { return caloEm_;             }
    double        caloHad()            const { return caloHad_;            }
    double        caloHo()             const { return caloHo_;             }
    float         iso03SumPt()         const { return iso03SumPt_;         }
    float         iso03emEt()          const { return iso03emEt_;          }
    float         iso03hadEt()         const { return iso03hadEt_;         }
    float         iso03hoEt()          const { return iso03hoEt_;          }
    int           iso03nJets()         const { return iso03nJets_;         }
    int           iso03nTracks()       const { return iso03nTracks_;       }
    double        normChi2()           const { return normChi2_;           }
    int           hitsNum()            const { return hitsNum_;            }
    double        impactParamXY()      const { return impactParamXY_;      }
    double        errorImpactParamXY() const { return errorImpactParamXY_; }

    void setCharge(             const int                       CHARGE             ) { charge_             = CHARGE;             }
    void setP4(                 const math::XYZTLorentzVector & P4                 ) { p4_                 = P4;                 }
    void setVertex(             const math::XYZPoint          & VERTEX             ) { vertex_             = VERTEX;             }
    void setCaloEm(             const double                  & CALOEM             ) { caloEm_             = CALOEM;             }
    void setCaloHad(            const double                  & CALOHAD            ) { caloHad_            = CALOHAD;            }
    void setCaloHo(             const double                  & CALOHO             ) { caloHo_             = CALOHO;             }
    void setIso03SumPt(         const float                   & ISO03SUMPT         ) { iso03SumPt_         = ISO03SUMPT;         }
    void setIso03emEt(          const float                   & ISO03EMET          ) { iso03emEt_          = ISO03EMET;          }
    void setIso03hadEt(         const float                   & ISO03HADET         ) { iso03hadEt_         = ISO03HADET;         }
    void setIso03hoEt(          const float                   & ISO03HOET          ) { iso03hoEt_          = ISO03HOET;          }
    void setIso03nJets(         const int                     & ISO03NJETS         ) { iso03nJets_         = ISO03NJETS;         }
    void setIso03nTracks(       const int                     & ISO03NTRACKS       ) { iso03nTracks_       = ISO03NTRACKS;       }
    void setNormChi2(           const double                  & NORMCHI2           ) { normChi2_           = NORMCHI2;           }
    void setHitsNum(            const int                       HITSNUM            ) { hitsNum_            = HITSNUM;            }
    void setImpactParamXY(      const double                  & IMPACTPARAMXY      ) { impactParamXY_      = IMPACTPARAMXY;      }
    void setErrorImpactParamXY( const double                  & ERRORIMPACTPARAMXY ) { errorImpactParamXY_ = ERRORIMPACTPARAMXY; }

  protected:
    int                     charge_;
    math::XYZTLorentzVector p4_;
    math::XYZPoint          vertex_;
    double                  caloEm_;
    double                  caloHad_;
    double                  caloHo_;
    float                   iso03SumPt_;
    float                   iso03emEt_;
    float                   iso03hadEt_;
    float                   iso03hoEt_;
    int                     iso03nJets_;
    int                     iso03nTracks_;
    double                  normChi2_;
    int                     hitsNum_;
    double                  impactParamXY_;
    double                  errorImpactParamXY_;
  };

  typedef std::vector<GlobalMuon> GlobalMuonCollection;

}

#endif // GLOBALMUON_H
