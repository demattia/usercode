#ifndef OFFLINEJET_H
#define OFFLINEJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <cmath>
#include <vector>

#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"

namespace anaobj {

  class SimpleTrack;
  
  /**
   *
   * Used for offline jets. Includes b-tagging informations.
   * Inherits from BaseJet.
   *
   * Author M. De Mattia - 8/11/2007
   *
   */

  class OfflineJet : public BaseJet {
  public:
    OfflineJet( const double         & ET, 
		const double         & ETA, 
		const double         & PHI, 
		const double         & UNCORRET, 
		const double         & EMENERGYFRACTION,
		const math::XYZTLorentzVector & P4,    
		const math::XYZPoint & VERTEX,
                const float          & DISCRIMINATORHIGHEFF, 
		const float          & DISCRIMINATORHIGHPUR, 
		const double         & JETMASS,
                const double         & BTAGTKINVMASS ) : BaseJet( ET, ETA, PHI ) {
      uncorrEt_             = UNCORRET;
      emEnergyFraction_     = EMENERGYFRACTION;
      p4_                   = P4;
      vertex_               = VERTEX;
      discriminatorHighEff_ = DISCRIMINATORHIGHEFF;
      discriminatorHighPur_ = DISCRIMINATORHIGHPUR;
      jetMass_              = JETMASS;
      bTagTkInvMass_        = BTAGTKINVMASS;
    }
    // Default constructor, only needed for classes.h
    OfflineJet() : BaseJet( 0., 0., 0. ) {
      uncorrEt_             = 0.;
      emEnergyFraction_     = 0.;
      p4_                   = math::XYZTLorentzVector(0.,0.,0.,0.);
      vertex_               = math::XYZPoint(0.,0.,0.);
      discriminatorHighEff_ = 0.;
      discriminatorHighPur_ = 0.;
      jetMass_              = 0.;
      bTagTkInvMass_        = 0.;
    }
    double         uncorrEt()             const { return uncorrEt_;             }
    double         emEnergyFraction()     const { return emEnergyFraction_;     }
    math::XYZTLorentzVector p4()          const { return p4_;                   }
    math::XYZPoint vertex()               const { return vertex_;               }
    float          discriminatorHighEff() const { return discriminatorHighEff_; }
    float          discriminatorHighPur() const { return discriminatorHighPur_; }
    double         jetMass()              const { return jetMass_;              }
    double         bTagTkInvMass()        const { return bTagTkInvMass_;        }
    const edm::RefVector<std::vector<SimpleTrack> > tkRefVec() const {
      return vecTrackRef_;
    }

    void setUncorrEt(             const double         & UNCORRET             ) { uncorrEt_             = UNCORRET;             }
    void setEmEnergyFraction(     const double         & EMENERGYFRACTION     ) { emEnergyFraction_     = EMENERGYFRACTION;     }
    void setP4(                   const math::XYZTLorentzVector & P4          ) { p4_                   = P4;                   }    
    void setVertex(               const math::XYZPoint & VERTEX               ) { vertex_               = VERTEX;               }
    void setDiscriminatorHighEff( const float          & DISCRIMINATORHIGHEFF ) { discriminatorHighEff_ = DISCRIMINATORHIGHEFF; }
    void setDiscriminatorHighPur( const float          & DISCRIMINATORHIGHPUR ) { discriminatorHighPur_ = DISCRIMINATORHIGHPUR; }
    void setJetMass(              const double         & JETMASS              ) { jetMass_              = JETMASS;              }
    void setbTagTkInvMass(        const double         & BTAGTKINVMASS        ) { bTagTkInvMass_        = BTAGTKINVMASS;        }
    void addTkRef( edm::Ref<std::vector<SimpleTrack> > TRACKREF ) {
      vecTrackRef_.push_back(TRACKREF);
    }

  protected:
    double         uncorrEt_;
    double         emEnergyFraction_;
    math::XYZTLorentzVector p4_;    
    math::XYZPoint vertex_;
    float          discriminatorHighEff_;
    float          discriminatorHighPur_;
    double         jetMass_;
    double         bTagTkInvMass_;
    edm::RefVector<std::vector<SimpleTrack> > vecTrackRef_;
  };

  typedef std::vector<OfflineJet> OfflineJetCollection;
  typedef edm::Ref<OfflineJetCollection> OfflineJetRef;
  typedef edm::RefProd<OfflineJetCollection> OfflineJetRefProd;
  typedef edm::RefVector<OfflineJetCollection> OfflineJetRefVector;
}

#endif // OFFLINEJET_H
