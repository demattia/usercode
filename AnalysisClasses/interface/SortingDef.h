#ifndef SORTINGDEF_H
#define SORTINGDEF_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "AnalysisExamples/AnalysisObjects/interface/SimpleTrack.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleCaloJet.h"

using namespace std;
using namespace anaobj;

//definizione dell'operatore minore per le coppie per il sort 
bool pairDecrescentSecondSort(const pair<int, double> & a, const pair<int, double> & b);


//definizione dell'operatore minore per le coppie per il sort  (DzAssociator)
bool pairSecondSort(const pair<int, double> & a, const pair<int, double> & b);

//definizione dell'operatore minore per il sort di una SimpleTrackCollection (RJAlg)
bool simpleTkCollSort(const SimpleTrack & a, const SimpleTrack & b) ;

struct mapFabsSort {   // (DzAssociator)
  bool operator() (const double & a, const double & b) {
    return fabs(a) < fabs(b);
  }
};
//defionizione di sort di un vettore di simplecalojet ordinati per valori decrescenti
// del rapporto delle probabilita' di assegnazione ai 2 vertici piu' vicini
bool simpleCaloJetCollSort(const SimpleCaloJet & a, const SimpleCaloJet & b) ;

//definizione di sort per un vettore di simple calojet ordinati per valori decrescenti di Et
bool EtSort(const SimpleCaloJet & a, const SimpleCaloJet & b); 

#endif
