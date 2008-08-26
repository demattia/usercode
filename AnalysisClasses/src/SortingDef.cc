#ifndef SORTINGDEF_CC
#define SORTINGDEF_CC

#include "AnalysisExamples/AnalysisClasses/interface/SortingDef.h"

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

using namespace std;

//definizione dell'operatore minore per le coppie per il sort 
bool pairDecrescentSecondSort(const pair<int, double> & a, const pair<int, double> & b) {
  return fabs(a.second) > fabs(b.second);
}

//definizione dell'operatore minore per le coppie per il sort (DzAssociator)
bool pairSecondSort(const pair<int, double> & a, const pair<int, double> & b) {
  return fabs(a.second) < fabs(b.second);
}

//definizione dell'operatore minore per il sort di una SimpleTrackCollection (RJAlg)
bool simpleTkCollSort(const SimpleTrack & a, const SimpleTrack & b) {
  return fabs(a.z()) < fabs(b.z());
}
//defionizione di sort di un vettore di simplecalojet ordinati per valori crescenti
// del rapporto delle probabilita' di assegnazione ai 2 vertici piu' vicini
bool simpleCaloJetCollSort(const SimpleCaloJet & a, const SimpleCaloJet & b) {
  return a.probBestVtx()/a.probSecondBestVtx() < b.probBestVtx()/b.probSecondBestVtx();
}
//definizione di sort per un vettore di simple calojet ordinati per valori crescenti di Et
bool EtSort(const SimpleCaloJet & a, const SimpleCaloJet & b) {
  return a.et() < b.et() ;
}


#endif
