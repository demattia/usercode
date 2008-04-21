#ifndef DZASSOCIATOR_H
#define DZASSOCIATOR_H

#include <vector>
#include <map>
#include <algorithm>


using namespace std;
//definizione dell'operatore minore per le coppie per il sort 
bool pairSecondSort(const pair<int, double> & a, const pair<int, double> & b) {
  return fabs(a.second) < fabs(b.second);
}

struct mapFabsSort {
  bool operator() (const double & a, const double & b) {
    return fabs(a) < fabs(b);
  }
};

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

//funzione che restituisce una coppia di coppie e vuole in input un vettore di XYZTLorentzVector o un vettore di SimVertex  e un double
//minimizza la differenza tra il double e tutti gli elementi del vettore
//restituisce la minore e la seconda minore differenza, e un contatore che identifica l'elemento del vettore cui e' associata quella differenza
//----------------------------------------------------------------------------------------------------------------------------------------------


pair<pair<int, double>, pair<int, double> > dzAssociator( const vector<math::XYZTLorentzVector> & vec_Vertices, const double & AvgW_Tk_dz ) {
  
  vector<pair<int, double> > dzVec;
  //itero sul vettore e riempio un vettore di coppie ( contatore, valore ) con tutte le differenze 
  int Vmin_Counter = 0;
  vector<math::XYZTLorentzVector>::const_iterator Vertices_it_1 = vec_Vertices.begin();
  for(; Vertices_it_1 != vec_Vertices.end() ; Vertices_it_1++, Vmin_Counter++) {
    dzVec.push_back( make_pair( Vmin_Counter, AvgW_Tk_dz - Vertices_it_1->z() ) );
  }
  return selectFirstSecond( dzVec );

}


// For the simulated verteces collection
pair<pair<int, double>, pair<int, double> > dzAssociator( const vector<SimVertex> & vecSimVertex, const double & zMax ) {
  
  map<double, int, mapFabsSort> mapDz;
  
  //itero sul vettore e riempio un vettore di coppie ( contatore, valore ) con tutte le differenze 
  int vCounter = 0;
  vector<SimVertex>::const_iterator verticesIt = vecSimVertex.begin();
  for(; verticesIt != vecSimVertex.end() ; verticesIt++, vCounter++) {
    
    mapDz.insert( make_pair( zMax - (verticesIt->position()).z(), vCounter )  );
    
  }
  
  cout<<"map size= "<<mapDz.size()<<endl;
  cout<<"map begin= "<<mapDz.begin()->first<<endl;
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
