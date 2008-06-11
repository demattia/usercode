#ifndef SIMPLETRACK_H
#define SIMPLETRACK_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseParticle.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <cmath>
#include <vector>

#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"

namespace anaobj {

  class OfflineJet;

  /**
   *
   * Used for tracks with vertex.
   * Inherits from BaseParticle.
   *
   * Author M. De Mattia - 10/4/2008
   *
   * Method zError()
   * Modified by Roberto Casagrande - 15/04/2008
   *
   * Further modifications -  11/06/2008
   * Added the 2D impact parameter significance
   *
   */
  
  class SimpleTrack : public BaseParticle {
  public:
    SimpleTrack( const double & PT, const double & ETA, const double & PHI, const double & Z, const double & ZERROR, const double & IP2DSIGNIFICANCE ) : BaseParticle( PT, ETA, PHI ) {
      z_ = Z;
      zerror_ = ZERROR;
      ip2Dsignificance_ = IP2DSIGNIFICANCE;
    }
    // Default constructor, only needed for classes.h
    SimpleTrack() : BaseParticle( 0., 0., 0. ) {
      z_ = 0.;
      zerror_ = 0.;
      ip2Dsignificance_ = 0.;
    }
    double z() const { return z_; }
    double zError() const { return zerror_; }
    double ip2Dsignificance() const { return ip2Dsignificance_; }
    const edm::RefVector<std::vector<OfflineJet> > jetRefVec() const {
      return vecJetRef_;
    }
    void setZ( const double & Z ) { z_ = Z; }
    void setZError( const double & ZERROR ) { zerror_ = ZERROR; }
    void ip2Dsignificance( const double IP2DSIGNIFICANCE ) { ip2Dsignificance_ = IP2DSIGNIFICANCE; }
    void addJetRef( edm::Ref<std::vector<OfflineJet> > JETREF ) {
      vecJetRef_.push_back(JETREF);
    }

  protected:
    double z_;
    double zerror_;
    double ip2Dsignificance_;
    edm::RefVector<std::vector<OfflineJet> > vecJetRef_;
  };

  typedef std::vector<SimpleTrack> SimpleTrackCollection;
  typedef edm::Ref<SimpleTrackCollection> SimpleTrackRef;
  typedef edm::RefProd<SimpleTrackCollection> SimpleTrackRefProd;
  typedef edm::RefVector<SimpleTrackCollection> SimpleTrackRefVector;
}

#endif // SIMPLETRACK_H
