#ifndef D2ASSOCIATOR_H
#define D2ASSOCIATOR_H

#include <vector>
#include <algorithm>
#include "AnalysisExamples/L1PixelAnalyzer/interface/DzAssociator.h"

using namespace std;


//funzione che restituisce una coppia di coppie e vuole in input una 
//SimpleTrackCollection e 2 double minimizza la differenza tra il D2
//e tutti gli elementi della collezionerestituisce la minore e la 
//seconda minore differenza, e un contatore che identifica l'elemento 
//del vettore cui e' associata quella differenza
//-------------------------------------------------------------------
pair<pair<int, double>, pair<int, double> > d2Associator( const anaobj::SimpleTrackCollection & vecSimTk, const double & eta, const double phi ) {
  int simTkCounter = 0;   
    vector< pair< int, double> > d2Vec;
      
      anaobj::SimpleTrackCollection::const_iterator vecSimTkIt = vecSimTk.begin(); 
	for(; vecSimTkIt !=  vecSimTk.end(); ++vecSimTkIt, ++simTkCounter ) { 
	    double d2 = 0.;
	    //Delta2 simTrack - reco selected track
	    d2 = pow( M_PI - fabs( fabs( vecSimTkIt->phi() - phi) - M_PI), 2) + pow( vecSimTkIt->eta() - eta,2);
	    d2Vec.push_back( make_pair( simTkCounter, d2 ) );
	  
	  }
  
  sort( d2Vec.begin(), d2Vec.end(), pairSecondSort );
	  
	  //se il vettore ha almeno 2 coppie ritorno le prime 2 coppie
	    if( d2Vec.size() > 1 ) {
	      return(make_pair( *(d2Vec.begin()), *(d2Vec.begin()+1) ) );
		}
	  //se non e' vuoto ritorno la prima e la seconda con valori fake
	    else if (d2Vec.size() > 0 ) {
	      return(make_pair( *(d2Vec.begin()), make_pair(-1,1000.) ) );
		}
	  //altrimenti ritorno valori fake
	    else {
	      return(make_pair( make_pair(-1,1000.), make_pair(-1,1000.) ) );
		}	  
	  
}


#endif
