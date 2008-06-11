#ifndef BASEVERTEX_H
#define BASEVERTEX_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/RefVector.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Simple vertex class storing the z and its corresponding error for a vertex.
   * Used to store the recoVertex collection by the OfflineProducer.
   *
   * Author M. De Mattia - 11/6/2008
   *
   */
  
  class BaseVertex {
  public:
    BaseVertex( const double & Z, const double & ZERROR ) {
      z_ = Z;
      zError_ = ZERROR;
    }
    // Default constructor, only needed for classes.h
    BaseVertex() {
      z_ = 0.;
      zError_ = 0.;
    }
    double z() const { return z_; }
    double zError() const { return zError_; }

    void setZ( const double & Z ) { z_ = Z; }
    void setZError( const double & ZERROR ) { zError_ = ZERROR; }

  protected:
    double z_;
    double zError_;
  };

  typedef std::vector<BaseVertex> BaseVertexCollection;
  typedef edm::Ref<BaseVertexCollection> BaseVertexRef;
}

#endif // BASEVERTEX_H
