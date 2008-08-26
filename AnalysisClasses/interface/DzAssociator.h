#ifndef DZASSOCIATOR_H
#define DZASSOCIATOR_H

#include <vector>
#include <map>
#include <algorithm>

#include "DataFormats/Math/interface/LorentzVector.h"
// For the SimVertex
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "AnalysisExamples/AnalysisClasses/interface/SortingDef.h"

using namespace std;

/**
 * Function DzAssociator
 * Overloaded function, accept as input vector of XYZTLorentzVector, or SimVertex and a double 
 * Minimize the differenze between the value and the z of vertices vector.
 * Define the sort for pair respect to the second element and the sort for the map repect the 
 * fabs of the key
 *
 * Author M.De Mattia - 16/04/2008
 *
 */

//funzione che restituisce una coppia di coppie e vuole in input un vettore di XYZTLorentzVector o un vettore di SimVertex  e un double
//minimizza la differenza tra il double e tutti gli elementi del vettore
//restituisce la minore e la seconda minore differenza, e un contatore che identifica l'elemento del vettore cui e' associata quella differenza
//----------------------------------------------------------------------------------------------------------------------------------------------

namespace {
  pair<pair<int, double>, pair<int, double> > selectFirstSecond( vector<pair<int, double> > & dzVec ) {
    
    //lo ordino dal piu' piccolo al piu' grande con il criterio dell'operatore minore definito prima 
    sort( dzVec.begin(), dzVec.end(), pairSecondSort );
    
    //se il vettore ha almeno 2 coppie ritorno le prime 2 coppie
    if( dzVec.size() > 1 ) {
      return(make_pair( *(dzVec.begin()), *(dzVec.begin()+1) ) );
    }
    //se non e' vuoto ritorno la prima e la seconda con valori fake
    else if (dzVec.size() > 0 ) {
      return(make_pair( *(dzVec.begin()), make_pair(-1,1000.) ) );
    }
    //altrimenti ritorno valori fake
    else {
      return(make_pair( make_pair(-1,1000.), make_pair(-1,1000.) ) );
    }
  }
}

template <class T> pair<pair<int, double>, pair<int, double> > dzAssociator( const vector<T> & vec_Vertices, const double & AvgW_Tk_dz ) {

  vector<pair<int, double> > dzVec;
  //itero sul vettore e riempio un vettore di coppie ( contatore, valore ) con tutte le differenze 
  int Vmin_Counter = 0;
  typename vector<T>::const_iterator Vertices_it_1 = vec_Vertices.begin();
  for(; Vertices_it_1 != vec_Vertices.end() ; Vertices_it_1++, Vmin_Counter++) {
    dzVec.push_back( make_pair( Vmin_Counter, AvgW_Tk_dz - Vertices_it_1->z() ) );
  }
  return selectFirstSecond( dzVec );
}


// For the simulated verteces collection
pair<pair<int, double>, pair<int, double> > dzAssociator( const vector<SimVertex> & vecSimVertex, const double & zMax );

#endif
