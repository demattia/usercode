#ifndef SIMPLECALOJET_H
#define SIMPLECALOJET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"

#include <cmath>
#include <vector>

namespace anaobj {

  /**
   *
   * Used for calojet with vertex.
   * Building associator algorithm 
   * Inherits from BaseJet.
   *
   * Author Roberto Casagrande - 28/4/2008
   *
   */
  
  class SimpleCaloJet : public BaseJet {
  public:
    SimpleCaloJet( const double & ET, const double & ETA, const double & PHI, const double & Z, const double & ZERROR , const int & IDJET,
		   const int IDBESTVTX, const int IDSECONDBESTVTX, const double & DZBESTVTX, const double &  DZSECONDBESTVTX,
		   const double & PROBBESTVTX, const double & PROBSECONDBESTVTX, const OfflineJet * OFFLINEJETPTR = 0 ) : BaseJet( ET, ETA, PHI ) {

      z_ = Z;
      zerror_ = ZERROR;
      idjet_ = IDJET;
      idBestVtx_ = IDBESTVTX;
      idSecondBestVtx_ = IDSECONDBESTVTX;
      dzBestVtx_ = DZBESTVTX; 
      dzSecondBestVtx_ = DZSECONDBESTVTX;
      probBestVtx_ = PROBBESTVTX; 
      probSecondBestVtx_ = PROBSECONDBESTVTX; 
      offlineJetPtr_ = OFFLINEJETPTR;

    }
      // Default constructor, only needed for classes.h
      SimpleCaloJet() : BaseJet( 0., 0., 0. ) {

	z_ = 0.;
	zerror_ = 0.;
	idjet_ = 0;
	idBestVtx_ = 0;
	idSecondBestVtx_ = 0;
	dzBestVtx_ = 0.; 
	dzSecondBestVtx_ = 0.;
	probBestVtx_ = 0.; 
	probSecondBestVtx_ = 0.; 
	offlineJetPtr_ = 0;

      }
	void setZ( const double & Z ) { z_ = Z; }
	void setZError( const double & ZERROR ) { zerror_ = ZERROR; }
	void setIdJet(const int & IDJET) {idjet_ = IDJET; }
 	void setIdBestVtx( const int IDBESTVTX) { idBestVtx_ = IDBESTVTX; }
	void setIdSecondBestVtx( const int IDSECONDBESTVTX){idSecondBestVtx_ = IDSECONDBESTVTX;}
     	void setDzBestVtx( const double & DZBESTVTX) {dzBestVtx_ = DZBESTVTX; }
	void setDzSecondBestVtx( const double &DZSECONDBESTVTX) {dzSecondBestVtx_ = DZSECONDBESTVTX;}
	void setProbBestVtx( const double & PROBBESTVTX) { probBestVtx_ = PROBBESTVTX; }
	void setProbSecondBestVtx ( const double & PROBSECONDBESTVTX ){ probSecondBestVtx_ = PROBSECONDBESTVTX; }
	void setOfflineJetPtr( const OfflineJet * OFFLINEJETPTR ) { offlineJetPtr_ = OFFLINEJETPTR; }

	double z() const { return z_; }
	double zError() const { return zerror_; }
	int idJet() const {return idjet_; }
 	int idBestVtx() const { return idBestVtx_; }
	int idSecondBestVtx() const { return idSecondBestVtx_; }
	double dzBestVtx() const { return dzBestVtx_;}
	double dzSecondBestVtx() const { return	dzSecondBestVtx_;}
	double probBestVtx() const { return probBestVtx_;}
	double probSecondBestVtx() const { return probSecondBestVtx_;}
	const OfflineJet * offlineJetPtr() const { return offlineJetPtr_; }

	
  protected:
	double z_;
	double zerror_;
	int idjet_;
       	int idBestVtx_, idSecondBestVtx_;
	double dzBestVtx_, dzSecondBestVtx_;
	double probBestVtx_, probSecondBestVtx_; 
	const OfflineJet * offlineJetPtr_;
  };
  
  typedef std::vector<SimpleCaloJet> SimpleCaloJetCollection;
  
}

#endif // SIMPLECALOJET_H
