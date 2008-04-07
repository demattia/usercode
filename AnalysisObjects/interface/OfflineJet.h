#ifndef OFFLINEJET_H
#define OFFLINEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Candidate/interface/Particle.h"
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
    OfflineJet( const double         & ET, 
		const double         & ETA, 
		const double         & PHI, 
		const double         & UNCORRET, 
		const double         & EMENERGYFRACTION,
		const math::XYZPoint & VERTEX,
                const float          & DISCRIMINATORHIGHEFF, 
		const float          & DISCRIMINATORHIGHPUR, 
		const double & JETMASS,
                const int      TKNUMS1, 
		const double & TKSUMPTS1, 
		const double & TAGTKMASSS1,
                const int      TKNUMS2, 
		const double & TKSUMPTS2, 
		const double & TAGTKMASSS2,
                const int      TKNUMS3, 
		const double & TKSUMPTS3, 
		const double & TAGTKMASSS3 ) : BaseJet( ET, ETA, PHI ) {
      uncorrEt_             = UNCORRET;
      emEnergyFraction_     = EMENERGYFRACTION;
      vertex_               = VERTEX;
      discriminatorHighEff_ = DISCRIMINATORHIGHEFF;
      discriminatorHighPur_ = DISCRIMINATORHIGHPUR;
      jetMass_              = JETMASS;

      tkNumS1_     = TKNUMS1;
      tkSumPtS1_   = TKSUMPTS1;
      tagTkMassS1_ = TAGTKMASSS1;

      tkNumS2_     = TKNUMS2;
      tkSumPtS2_   = TKSUMPTS2;
      tagTkMassS2_ = TAGTKMASSS2;

      tkNumS3_     = TKNUMS3;
      tkSumPtS3_   = TKSUMPTS3;
      tagTkMassS3_ = TAGTKMASSS3;
    }
    // Default constructor, only needed for classes.h
    OfflineJet() : BaseJet( 0., 0., 0. ) {
      uncorrEt_             = 0.;
      emEnergyFraction_     = 0.;
      vertex_               = math::XYZPoint(0.,0.,0.);
      discriminatorHighEff_ = 0.;
      discriminatorHighPur_ = 0.;
      jetMass_              = 0.;

      tkNumS1_     = 0;
      tkSumPtS1_   = 0.;
      tagTkMassS1_ = 0.;

      tkNumS2_     = 0;
      tkSumPtS2_   = 0.;
      tagTkMassS2_ = 0.;

      tkNumS3_     = 0;
      tkSumPtS3_   = 0.;
      tagTkMassS3_ = 0.;
    }
    double         uncorrEt()             const { return uncorrEt_;             }
    double         emEnergyFraction()     const { return emEnergyFraction_;     }
    math::XYZPoint vertex()               const { return vertex_;               }
    float          discriminatorHighEff() const { return discriminatorHighEff_; }
    float          discriminatorHighPur() const { return discriminatorHighPur_; }
    double         jetMass()              const { return jetMass_;              }

    int    tkNumS1()     const { return tkNumS1_;     }
    double tkSumPtS1()   const { return tkSumPtS1_;   }
    double tagTkMassS1() const { return tagTkMassS1_; }

    int    tkNumS2()     const { return tkNumS2_;     }
    double tkSumPtS2()   const { return tkSumPtS2_;   }
    double tagTkMassS2() const { return tagTkMassS2_; }

    int    tkNumS3()     const { return tkNumS3_;     }
    double tkSumPtS3()   const { return tkSumPtS3_;   }
    double tagTkMassS3() const { return tagTkMassS3_; }

    void setUncorrEt(             const double         & UNCORRET             ) { uncorrEt_             = UNCORRET;             }
    void setEmEnergyFraction(     const double         & EMENERGYFRACTION     ) { emEnergyFraction_     = EMENERGYFRACTION;     }
    void setVertex(               const math::XYZPoint & VERTEX               ) { vertex_               = VERTEX;               }
    void setDiscriminatorHighEff( const float          & DISCRIMINATORHIGHEFF ) { discriminatorHighEff_ = DISCRIMINATORHIGHEFF; }
    void setDiscriminatorHighPur( const float          & DISCRIMINATORHIGHPUR ) { discriminatorHighPur_ = DISCRIMINATORHIGHPUR; }
    void setJetMass(              const double         & JETMASS              ) { jetMass_              = JETMASS;              }

    void setTkNumS1(     const int      TKNUMS1     ) { tkNumS1_     = TKNUMS1;     }
    void setTkSumPtS1(   const double & TKSUMPTS1   ) { tkSumPtS1_   = TKSUMPTS1;   }
    void setTagTkMassS1( const double & TAGTKMASSS1 ) { tagTkMassS1_ = TAGTKMASSS1; }

    void setTkNumS2(     const int      TKNUMS2     ) { tkNumS2_     = TKNUMS2;     }
    void setTkSumPtS2(   const double & TKSUMPTS2   ) { tkSumPtS2_   = TKSUMPTS2;   }
    void setTagTkMassS2( const double & TAGTKMASSS2 ) { tagTkMassS2_ = TAGTKMASSS2; }

    void setTkNumS3(     const int      TKNUMS3     ) { tkNumS3_     = TKNUMS3;     }
    void setTkSumPtS3(   const double & TKSUMPTS3   ) { tkSumPtS3_   = TKSUMPTS3;   }
    void setTagTkMassS3( const double & TAGTKMASSS3 ) { tagTkMassS3_ = TAGTKMASSS3; }

  protected:
    double         uncorrEt_;
    double         emEnergyFraction_;
    math::XYZPoint vertex_;
    float          discriminatorHighEff_;
    float          discriminatorHighPur_;
    double         jetMass_;

    int    tkNumS1_;
    double tkSumPtS1_;
    double tagTkMassS1_;

    int    tkNumS2_;
    double tkSumPtS2_;
    double tagTkMassS2_;

    int    tkNumS3_;
    double tkSumPtS3_;
    double tagTkMassS3_;
  };

  typedef std::vector<OfflineJet> OfflineJetCollection;

}

#endif // OFFLINEJET_H
