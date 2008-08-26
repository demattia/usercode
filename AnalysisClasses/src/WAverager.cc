#ifndef WAVERAGER_CC
#define WAVERAGER_CC

#include "AnalysisExamples/AnalysisClasses/interface/WAverager.h"

#include <vector>
#include <algorithm>

#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"

using namespace std;

/**
 * Used to make a weighted average for a simpletrackcollection z value
 * 
 *
 *  Author Roberto Casagrande - 15/04/2008
 */


pair< int, pair< double, double > > wAverager(const anaobj::SimpleTrackCollection & vecTk){
  
  double sumTk  = 0.; //somma delle coordinate z delle tracce
  double sumWeight = 0.; //somma pesi
  double wAvg = 0.; //media pesata
  double wAvgError = 0.; //errore media pesata
  int nTk = 0; //contatore # valori
  
  anaobj::SimpleTrackCollection::const_iterator vecTkIt = vecTk.begin();
  for( ; vecTkIt  != vecTk.end(); ++vecTkIt){ 
    //sommo i valori
    sumTk += vecTkIt->z()/pow(vecTkIt->zError(),2);
    //sommatoria dei pesi
    sumWeight += 1/pow(vecTkIt->zError(),2);
    //conto il # di valori  
    nTk++;
  }//end loop 
  
  //se il # di valori sommati e' diverso da 0 faccio la media pesata
  if (nTk != 0){ 
    wAvg = sumTk /sumWeight; // Z Jet
    wAvgError = sqrt( 1/sumWeight ); // Z Jet error
    
    return (make_pair ( nTk , make_pair( wAvg , wAvgError) ) ); 
  }  
  //  else {
    return (make_pair ( 0 , make_pair( 1000. , 1000. ) ) );
    //  }
  
}  


#endif
