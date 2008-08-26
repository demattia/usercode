#ifndef DZASSOCIATOR_CC
#define DZASSOCIATOR_CC

#include "AnalysisExamples/AnalysisClasses/interface/DzAssociator.h"

#include <vector>
#include <map>
#include <algorithm>

// For the SimVertex
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

/**
 * Function DzAssociator
 * Overloaded function, accept as input vector of template, or SimVertex and a double 
 * Minimize the differenze between the value and the z of vertices vector.
 * Define the sort for pair respect to the second element and the sort for the map repect the 
 * fabs of the key
 *
 * Author M.De Mattia - 16/04/2008
 *
 */

using namespace std;

//funzione che restituisce una coppia di coppie e vuole in input un vettore di XYZTLorentzVector o un vettore di SimVertex  e un double
//minimizza la differenza tra il double e tutti gli elementi del vettore
//restituisce la minore e la seconda minore differenza, e un contatore che identifica l'elemento del vettore cui e' associata quella differenza
//----------------------------------------------------------------------------------------------------------------------------------------------


// For the simulated vertexes collection
pair<pair<int, double>, pair<int, double> > dzAssociator( const vector<SimVertex> & vecSimVertex, const double & zMax ) {
  
  map<double, int, mapFabsSort> mapDz;
  
  //itero sul vettore e riempio un vettore di coppie ( contatore, valore ) con tutte le differenze 
  int vCounter = 0;
  vector<SimVertex>::const_iterator verticesIt = vecSimVertex.begin();
  for(; verticesIt != vecSimVertex.end() ; verticesIt++, vCounter++) {
    
    mapDz.insert( make_pair( zMax - (verticesIt->position()).z(), vCounter )  );
    
  }
  
  //  cout<<"map size= "<<mapDz.size()<<endl;
  //  cout<<"map begin= "<<mapDz.begin()->first<<endl;
  //se il vettore ha almeno 2 coppie ritorno le prime 2 coppie
  if( mapDz.size() > 1 ) {
    
    map<double, int>::const_iterator mapDzIt = mapDz.begin();
    mapDzIt++;
    
    return(make_pair( make_pair( mapDz.begin()->second, mapDz.begin()->first ), make_pair( mapDzIt->second, mapDzIt->first ) ) );
  }
  //se non e' vuoto ritorno la prima e la seconda con valori fake
  else if (  mapDz.size() > 0 ) {
    return(make_pair(  make_pair( mapDz.begin()->second, mapDz.begin()->first ), make_pair(-1,1000.) ) );
  }
  //altrimenti ritorno valori fake
  else {
    return(make_pair( make_pair(-1,1000.), make_pair(-1,1000.) ) );
  }
  
}


#endif
