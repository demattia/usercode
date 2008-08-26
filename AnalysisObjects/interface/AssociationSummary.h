#ifndef ASSOCIATIONSUMMARY_H
#define ASSOCIATIONSUMMARY_H

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Used for jet association to vertex
   * 
   *
   * Author Roberto Casagrande - 06/05/2008
   *
   * 
   * 
   *
   */
  
  class AssociationSummary {
  public:
    AssociationSummary( const int & IDJET, const int & IDVTX, const double & BEST, const int & IDSECONDBEST, const double & SECONDBEST, const bool & ASSOCIATED ) {

      idjet_ = IDJET;
      idvtx_ = IDVTX;
      best_ = BEST;
      idsecondbest_ = IDSECONDBEST;
      secondbest_ = SECONDBEST;
      associated_ = ASSOCIATED;      
    }
    // Default constructor, only needed for classes.h
    AssociationSummary(){
      
      idjet_ = 0;
      idvtx_ = 0; 
      best_ = 0.;
      idsecondbest_ = 0;
      secondbest_ = 0.;
      associated_ = false;
      
    }

    int idJet() const { return idjet_; }
    int idVtx() const { return idvtx_; }
    double Best() const { return best_; }
    int idSecondBest() const { return idsecondbest_; }
    double SecondBest() const { return secondbest_; }
    bool Associated() const { return associated_; }
    
    void setIdJet( const int & IDJET ) { idjet_ = IDJET; }
    void setIdVtx( const int & IDVTX ) { idvtx_ = IDVTX; }
    void setBest  ( const double &  BEST )  { best_ = BEST;}
    void setIdSecondBest  ( const int & IDSECONDBEST ) { idsecondbest_ = IDSECONDBEST;}
    void setSecondBest ( const double & SECONDBEST ) { secondbest_ = SECONDBEST;}
    void setAssociated ( const bool &  ASSOCIATED ) { associated_ = ASSOCIATED;}


	
  protected:
 
    int idjet_;
    int idvtx_;
    double best_;
    int idsecondbest_; 
    double secondbest_;
    bool associated_;

  };
  
  typedef std::vector<AssociationSummary> AssociationSummaryCollection;
  
}

#endif 
